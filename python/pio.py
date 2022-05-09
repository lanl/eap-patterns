#!/usr/bin/env python3
"""
========================================================================================
 (C) (or copyright) 2022. Triad National Security, LLC. All rights reserved.

 This program was produced under U.S. Government contract 89233218CNA000001 for Los
 Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC
 for the U.S. Department of Energy/National Nuclear Security Administration. All rights
 in the program are reserved by Triad National Security, LLC, and the U.S. Department
 of Energy/National Nuclear Security Administration. The Government is granted for
 itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
 license in this material to reproduce, prepare derivative works, distribute copies to
 the public, perform publicly and display publicly, and to permit others to do so.
========================================================================================

To get latest official version of this script:
  git clone git@gitlab.com:lanl/piotools.git

A simple class to read and manipulate PIO files.
Demo at bottom expands fractional volume and fractional energy.

Author: Sriram Swaminarayan (sriram@lanl.gov)
Date: December 02, 2021
Version: 1.01

Change Log:
Version 1.01: 2021-12-02: Introduced sizeMatArrays to recognize fractional arrays
                          Set buffer size to 8GB for copyToOffset()
Version 1.00: 2021-03-01: Initial commit
"""
from __future__ import print_function
import numpy as np
import struct

try:
    import gzip
except:
    # ignore failure, we will do simple open
    pass


class pio:
    """

    This class will read PIO files and allow minimal manipulation.
    It is expected that other scripts will build on the PIO class
    to do heavy-duty lifting.

    Author: Sriram Swaminarayan (sriram@lanl.gov)
    Date: December 02, 2021
    Version: 1.01

    """

    def __init__(self, theFile, verbose=0):
        """
        Initializes a PIO class and returns a PIO object.
        Argument is the PIO file to be read.
        """
        self.verbose = verbose

        self.offset = 0

        self.gzip = 0
        try:
            self.fp = gzip.open(theFile, mode="rb")
            self.gzip = 1
            s = self.str()

            self.seek(0)
        except:
            if theFile.endswith(".gz"):
                raise ValueError(".gz files not supported: gzip module not found")

            if self.gzip:
                self.fp.close()
            self.gzip = 0
            self.fp = open(theFile, mode="rb")

        # read signature
        s = self.str(8)

        if self.verbose:
            print(f"  File type is: {s}")
        if not s.lower() == b"pio_file":
            raise ValueError("Invalid file type")

        # Ensure second value is two
        d = self.doubles()
        if d != 2.0:
            raise ValueError("Second value is not 2")

        # read version
        self.version = self.ints()
        if self.verbose:
            print(f"version={self.version}")

        # read element lengths
        self.lName = self.ints()
        self.lHeader = self.ints()
        self.lIndex = self.ints()

        # read data/time
        self.date = self.str(16)
        if self.verbose:
            print(self.date)

        # read number of variables and index location
        self.n = self.ints()
        self.position = self.ints()

        # read file signature
        self.signature = self.ints()

        # read the variable index
        self.names = {}
        self.xnames = []
        if self.verbose:
            print("position=", self.position)
        self.seek(self.position)
        for i in range(int(self.n)):
            hdf = self.readArrayHeader()
            idx = hdf["name"].strip() + b"_%d" % hdf["index"]
            self.names[idx.decode()] = hdf
            self.xnames.append(hdf)

        # get the number of dimensions based on cell data
        amhc_i = self.readArray("amhc_i_0").astype(np.int64)
        self.nummat = amhc_i[5]
        self.ndim = amhc_i[42]
        self.numcell = amhc_i[54] 
        amhc_i = None
        
        self.outOffset = -1

        # initialize the CSR data
        self.csrN = 0
        self.csrLen = 0
        self.csrVol = None
        self.csrInvVol = None
        self.csrID = None
        self.csrIdx = None

    def updateCsrIndices(self, csr_counts, csr_ids, csr_vols, shift=1):

        if self.csrIdx is not None:
            # Already inited, do nothing
            return

        print("updating")
        if not csr_counts.endswith("_0"):
            csr_counts += "_0"

        if not csr_ids.endswith("_0"):
            csr_ids += "_0"

        if not csr_vols.endswith("_0"):
            csr_vols += "_0"

        # read volume and invert
        self.csrVol = self.readArray(csr_vols)
        self.csrInvVol = 1.0 / self.csrVol

        # read, convert, and shift CSR IDs
        self.csrID = self.readArray(csr_ids).astype(np.int64)
        self.csrN = max(self.csrID)
        self.csrID -= shift

        # read counts and convert to indices
        self.csrIdx = self.readArray(csr_counts).astype(np.int64)  # read array
        self.csrLen = sum(self.csrIdx)
        self.csrIdx = np.append(self.csrIdx, 0)  # extend array
        iSum = 0  # indexing starts at 0
        for iCell in range(self.numcell):
            iNext = self.csrIdx[iCell]
            self.csrIdx[iCell] = iSum
            iSum += iNext
        self.csrIdx[self.numcell] = iSum

    def expandCsrVariable(self, name, scale=False):
        """
        Returns the expanded "index" version of the variable
        """
        if self.csrN == 0:
            return None

        if not name.endswith("_0"):
            name += "_0"
        data = self.readArray(name)
        print(f'read {name}={data}')        

        newArray = np.zeros((self.csrN, self.numcell))
        for iCell in range(self.numcell):
            iStart = self.csrIdx[iCell]
            iEnd = self.csrIdx[iCell + 1]
            for idx in range(iStart, iEnd):
                iCsr = self.csrID[idx]
                if scale:
                    newArray[iCsr][iCell] = data[idx] * self.csrInvVol[iCell]
                else:
                    newArray[iCsr][iCell] = data[idx]

        return newArray

    def writeHeader(self, fp, n=None, pos=None):
        """
        Writes PIO header to the file pointer sent.
        """
        if n is None:
            n = self.n

        if pos is None:
            pos = self.position
            
        self.outOffset = 0
        fp.write(b"pio_file")
        np.array(
            [2.0, self.version, self.lName, self.lHeader, self.lIndex],
            dtype="double",
        ).tofile(fp)
        fp.write(self.date)
        np.array([n, pos, self.signature], dtype="double").tofile(fp)
        np.zeros((self.lHeader - 11), dtype="double").tofile(fp)
        self.outOffset = self.lHeader

    def writeWithNewCellArray(self, outName, newName, newData):
        """
        Adds the new cell array with name `newName` and data
        `newData` to the PIO file and writes it to file named
        `outname`.
        """

        # now add in a new array
        self.oldPosition = self.position

        self.addCellArray(newName)

        with open(outName, "wb") as ofp:
            # write the header
            self.writeHeader(ofp)

            # copy rest of file till old Index offset
            self.copyToOffset(ofp, self.lHeader, self.oldPosition)

            # write new data to file
            newData.tofile(ofp)
            self.outOffset += len(newData)

            # write Index
            self.writeIndex(ofp)

            ofp.close()
        self.outOffset = -1

    def writeWithExpandedCsrArray(self, outName, myVars=None, outIndices=None):
        """
        Adds all chunk_ variables to the PIO file and
        writes it to file named `outname`.
        """

        # save old offset
        self.oldPosition = self.position

        # Update material indices
        self.updateCsrIndices("chunk_nummat", "chunk_mat", "vcell", 1)

        # Add variables to the list
        if myVars is None:
            myVars = []
            for name in self.names:
                if self.names[name]["length"] == self.csrLen:
                    if name.endswith("_0"):
                        name = name[:-2]
                    myVars.append(name)

        if outIndices is None:
            outIndices = [x + 1 for x in range(self.csrN)]

        for name in myVars:
            # add in an entry per CSR Index
            for iCsr in outIndices:
                self.addCellArray(f"{name}-{iCsr}")

        # open file and write data
        with open(outName, "wb") as ofp:
            # write the header
            self.writeHeader(ofp)

            # copy rest of file till old Index offset
            self.copyToOffset(ofp, self.lHeader, self.oldPosition)

            # write new data to file
            for name in myVars:
                theVar = p.expandCsrVariable(name, True)
                print(name, theVar.shape)
                for iCsr in outIndices:
                    newData = theVar[iCsr - 1, :]
                    print(
                        "  Writing:",
                        name + f"-{iCsr}",
                        iCsr,
                        newData.shape,
                        newData.dtype,
                    )
                    newData.tofile(ofp)
                    self.outOffset += self.numcell

            # write Index
            self.writeIndex(ofp)

            ofp.close()
        self.outOffset = -1

    def writeIndex(self, outfp):
        """
        Writes the index of variables to outfp.
        This is the trailer to the PIO file.
        """
        for x in self.xnames:
            outfp.write(self.xnames[x]["bytes"])
            self.outOffset += self.lIndex

    def copyToOffset(self, outfp, offsetStart, offsetEnd):
        """
        Copies data verbatim from self.fp to outfp from
        offsetStart to offsetEnd in self.fp.

        The copying is done in self.numcell double sized
        blocks unless cells are less than 1M, in which
        case it is done in blocks of size 1M doubles.
        """
        self.seek(offsetStart)
        sz = offsetEnd - offsetStart

        # Read / write at most 8GB at once
        bufsize = 1024 ** 3

        written = 0
        left = sz
        while left:
            if left < bufsize:
                bufsize = left
            buf = self.doubles(bufsize)
            buf.tofile(outfp)
            written += bufsize
            left -= bufsize
            if self.verbose:
                print(f"{0.95*(100.*written)/sz:.2f}%  written ")
        self.outOffset += sz

    def addCellArray(self, name, idx=0, copyFrom="pres_0"):
        """
        Appends metadata for a new variable to the list of
        variables (self.names, self.xnames).  The base data
        is copied from the variable 'copyFrom' which defaults
        to pres_0.

        No checks are done on copyFrom, and swift death can
        result if you are not careful.
        """
        cch = self.copyArrayHeader(self.names[copyFrom])
        fmt = "%%-%ds" % self.lName
        cch["name"] = bytes(name, "utf8")
        longName = bytes(fmt % name, "utf8")
        cch["offset"] = self.lIndex
        b = cch["bytes"]
        o = self.lName
        cch["bytes"] = longName + bytes(
            np.array([0, self.numcell, self.position, 0], dtype="double")
        )
        if self.verbose:
            print(cch["bytes"], len(cch["bytes"]))
        self.position += self.numcell
        self.n += 1
        self.names[f"{name}_{idx}"] = cch
        self.xnames.append(cch)
        return

    def readArray(self, name):
        """
        Reads given array from the file and
        returns it as an array of doubles.

        Returns None if array is not found.
        """
        if name not in self.names:
            return None
        hdr = self.names[name]
        self.seek(hdr["offset"])
        data = self.doubles(hdr["length"], force=True)
        return data

    def readArrayInt(self, name):
        """
        Reads given array from the file and
        returns it as an array of doubles.

        Returns None if array is not found.
        """
        if name not in self.names:
            return None
        hdr = self.names[name]
        self.seek(hdr["offset"])
        data = self.ints(hdr["length"], force=True)
        return data

    def readArrayRange(self, name, iStart, N, force=False, ints=False):
        """
        Reads N entries from given array starting at iStart
        Returns it as an array of doubles.

        Returns None if array is not found.
        """
        if name not in self.names:
            return None
        hdr = self.names[name]
        self.seek(hdr["offset"] + iStart)
        if ints:
            data = self.ints(N, force=True)
        else:
            data = self.doubles(N, force=True)
        return data

    def readArrayHeader(self):
        """
        reads array at current header
        """
        x = self.lName
        start = self.offset
        data = self.str(8 * self.lIndex)
        name = data[:x]
        index = int(struct.unpack("d", data[x : x + 8])[0])
        length = int(struct.unpack("d", data[x + 8 : x + 16])[0])
        offset = int(struct.unpack("d", data[x + 16 : x + 24])[0])
        return {
            "name": name,
            "index": index,
            "length": length,
            "offset": offset,
            "bytes": data,
        }

    def copyArrayHeader(self, src):
        """ Returns a copy of the array header """
        ret = {}
        for x in src:
            ret[x] = src[x]
        return ret

    def seek(self, offset):
        """ Seeks to offset in the input file """
        self.offset = offset
        offset = int(8 * offset)
        self.fp.seek(offset)

    def str(self, count=8, offset=None):
        """ Reads count bytes from the input file """
        count = int(count)
        if offset is not None:
            self.seek(offset)
        s = self.fp.read(count)

        # offset is counted in doubles
        self.offset += int(count / 8)
        return s

    def doubles(self, count=1, offset=None, force=False):
        """
        Reads count doubles from the input file.
        If count == 1 and force is False, then it will
        return scalars.
        """
        count = int(count)
        if offset is not None:
            self.seek(offset)
        if self.gzip:
            value = np.frombuffer(self.fp.read(count * 8), dtype="double", count=count)
        else:
            value = np.fromfile(self.fp, dtype="double", count=count)
        self.offset += count
        if count == 1 and (not force):
            return value[0]
        return value

    def ints(self, count=1, offset=None, force=False):
        """
        Reads count doubles from the input file and returns as ints.
        If count == 1 and force is False, then it will
        return scalars.
        """
        if count == 1 and (not force):
            value = int(self.doubles(count, offset, force))
        else:
            value = [int(x) for x in self.doubles(count, offset, force)]
        return value


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print(
            f"""
        {sys.argv[0]} takes at least one argument.
        """
        )
    else:
        filename = sys.argv[1]
        p = pio(filename)
        # p.updateCsrIndices("chunk_nummat", "chunk_mat", "vcell", 1)
        # print("CSR variables are:")
        for n in p.names:
            print("    ", n, p.names[n]['offset'], p.names[n]['length'])
        # c = p.readArray("hist_cycle_0")
        # print(c)
        # c = p.readArray("hist_time_0")
        # print(c)
        nbrs = [None] * (2 * p.ndim)
        for i in range(2 * p.ndim):
            nbrs[i] = p.readArray(f"cell_index_{i+1}")
            print(nbrs[i])

        print(f'numdim={p.ndim}, nummat={p.nummat}, numcell={p.numcell}')
        # Will write chunk_vol and chunk_eng for all CSR Indices
        # outName = "bigfile-dmp000000"
        # myVars = ["chunk_vol", "chunk_eng"]
        # outIndices = None
        # # p.writeWithExpandedCsrArray(outName, myVars, outIndices)

        # # # To write all CSR variables to file for all materials:
        # # p.writeWithExpandedCsrArray(outName)
