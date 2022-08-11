#!/usr/bin/env python3
"""
========================================================================================
 (C) (or copyright) 2021. Triad National Security, LLC. All rights reserved.

 This program was produced under U.S. Government contract 89233218CNA000001 for Los
 Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC
 for the U.S. Department of Energy/National Nuclear Security Administration. All rights
 in the program are reserved by Triad National Security, LLC, and the U.S. Department
 of Energy/National Nuclear Security Administration. The Government is granted for
 itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
 license in this material to reproduce, prepare derivative works, distribute copies to
 the public, perform publicly and display publicly, and to permit others to do so.
========================================================================================

Takes one argument: PIO file to convert

Reads in the PIO file, changes neighbors of coarse cell at
AMR boundaries to point at real neighbor and writes the
binary format for EAP patterns:
bytes   kind   what
16      char   signature: eap-patterns-bin
8       int64  endian marker: 1
8       int64  number of dimensions, 
8       int64  number of cells
8       int64  nLen= length of space padded name strings
8       int64  number of variables

nLen    char   name of variable 1
8       int64  offset

nLen    char   name of variable 2
8       int64  offset

....

nLen    char   name of variable N
8       int64  offset

<data for variables at offsets above>

Characters: trailer: Contents: eap-patterns-bin, 1(i64), ndim(i64), ncell(i64), name_len(i64), nVars(i64), list of names + offsets(i64), data(doubles)\n"

** REQUIRES Python3 **

Author: Sriram Swaminarayan (sriram@lanl.gov)
Date: August 10, 2022
Version: 1.0

"""
import sys, os
import numpy as np
basepath = os.path.dirname(os.path.abspath(__file__))
sys.path.append(basepath + "/..")
from pio import pio
import re
import struct


def convertPIO(fname, outfile, verbose=False):

    myNames = re.compile("cell_center_" +
                         "|cell_daughter_" +
                         "|cell_index_" +
                         "|cell_level_" +
                         "|vcell_")
    p = pio(fname)

    daughter = p.readArray("cell_daughter_0").astype(np.int64)
    offsets = [1,0,2,0,4,0]

    if(verbose):
        print(len(daughter),daughter)
    nbrs = [None] * (2 * p.ndim)
    for idx in range(2 * p.ndim):
        if(verbose):
            print(f"cell_index_{idx+1}")
        a = p.readArray(f"cell_index_{idx+1}").astype(np.int64)
        for iCell in range(p.numcell):
            mydtr = daughter[iCell]
            if mydtr == 0:
                iNbr = a[iCell] - 1
                dtr = daughter[iNbr]
                if dtr > 0:
                    a[iCell] = (daughter[iNbr] + offsets[idx])

        nbrs[idx] = a.astype(np.double)
    daughter = None
    newNames = {}
    
    for n in p.names:
        m = myNames.match(n)
        if m is not None:
            newNames[n] = 0

    # Now to write the file
    with open(outfile,'wb') as ofp:
        if(verbose):
            print('   writing header')
        ofp.write(b'eap-patterns-bin')
        ofp.write(struct.pack("q",1))
        ofp.write(struct.pack("q",p.ndim))
        ofp.write(struct.pack("q",p.numcell))
        ofp.write(struct.pack("q",p.lName))
        ofp.write(struct.pack("q",len(newNames)))

        # Offset is current position + name list
        offset = ofp.tell() + ( p.lName + 8 ) * len(newNames)

        # Write variable offsets
        if(verbose):
            print('   writing offsets')
        for n in newNames:
            newNames[n] = offset
            ofp.write(f"{n:<{p.lName}}".encode())
            ofp.write(struct.pack("q",offset))
            offset += 8 * p.numcell

        # Write data
        for n in newNames:
            if(verbose):
                print('   writing: ',n.strip())
            if n.startswith('cell_index_'):
                idx = int(n.strip()[-1]) - 1
                nbrs[idx].tofile(ofp)
            else:
                c = p.readArray(n)
                c.tofile(ofp)

        # Write variable data
        trailer = "\nContents: eap-patterns-bin, 1(i64), ndim(i64), ncell(i64), name_len(i64), nVars(i64), list of names + offsets(i64), data(doubles)\n"
        ofp.write(trailer.encode())

        if(verbose):
            print('done writing')
        
if __name__ == "__main__":
    import sys

    myFile = sys.argv[1]
    if len(sys.argv) > 2:
        outfile = sys.argv[2]
    else:
        outfile = "out.bin"

    convertPIO(myFile, outfile, True)
