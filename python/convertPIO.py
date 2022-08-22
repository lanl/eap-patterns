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
8       int64  element size 1
8       int64  offset 1

nLen    char   name of variable 2
8       int64  element size 2
8       int64  offset 2

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


def getFaceType(index, iCell, iNbr, level):
    if iCell == iNbr:
        iType = 1 + index%2
    elif level[iNbr] == level[iCell]:
        iType = 3
    elif level[iNbr] > level[iCell]:
        iType = 4 + index%2
    elif level[iNbr] < level[iCell]:
        iType = 5 - index%2
    else:
        raise ValueError("Invalid type calculation")
    return iType
        
def convertPIO(fname, outfile, verbose=False):

    p = pio(fname)

    cell_level = p.readArray("cell_level_0").astype(np.int32)
    daughter = p.readArray("cell_daughter_0").astype(np.int64)
    offsets = [1,0,2,0,4,0]
    nFaces = np.zeros((2 * p.ndim,6),'i')

    if(verbose):
        print(len(daughter),daughter)
    nbrs = np.zeros((2 * p.ndim, p.numcell), np.int64)
    face_type = np.zeros((2 * p.ndim, p.numcell), np.int8)
    numtop = 0
    for idx in range(2 * p.ndim):
        idim = int(idx/2)
        if(verbose):
            print(f"cell_index_{idx+1}")
        a = p.readArray(f"cell_index_{idx+1}").astype(np.int64)
        for iCell in range(p.numcell):
            mydtr = daughter[iCell]
            iType = 0
            if mydtr == 0:
                if idx ==0:
                    numtop += 1
                iNbr = a[iCell] - 1
                dtr = daughter[iNbr]
                if dtr > 0:
                    a[iCell] = (daughter[iNbr] + offsets[idx])
                    iNbr = a[iCell] - 1
                iType = getFaceType(idx, iCell, iNbr, cell_level)
                face_type[idx,iCell] = iType
                if iType < 3:
                    nFaces[idx,iType] += 1
                elif iType == 3 and idx%2 == 0:
                    nFaces[idx,iType] += 1
                elif iType == 4 and idx%2 > 0:
                    nFaces[idx,iType] += 1
                elif iType == 5 and idx%2 == 0:
                    nFaces[idx,iType] += 1
        nbrs[idx, :] = a
    print('numtop=', numtop)
    # variables we will write
    sizeOf = {'i08': 1, 'i32':4, 'i64':8, 'f64':8}
    newNames = {}
    newNames['cell_daughter'] = {'type':'i64', 'n':p.numcell}
    newNames['cell_level'] = {'type':'i32', 'n':p.numcell}
    newNames['cell_index'] = {'type':'i64', 'n':2 * p.ndim * p.numcell}
    newNames['face_type'] = {'type':'i08', 'n':2 * p.ndim * p.numcell}
    newNames['cell_center'] = {'type':'f64', 'n':p.ndim * p.numcell}
    newNames['vcell'] = {'type':'f64', 'n':p.numcell}

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

        # Write variable offsets        
        if(verbose):
            print('   writing offsets')
            
        # Offset is current position + size of  name list
        offset = ofp.tell() + (p.lName + 8 + 8) * len(newNames)
        for n in newNames:
            isz = sizeOf[newNames[n]['type']]
            ofp.write(f"{n:<{p.lName}}".encode())
            ofp.write(struct.pack("q",isz))
            ofp.write(struct.pack("q",offset))
            offset += isz * newNames[n]['n']

        # Write data
        if(verbose):
            print('   writing variables')
        for n in newNames:
            if(verbose):
                print('   ',n)
            if n == 'cell_center':
                for idx in range(p.ndim):
                    name = f'{n}_{idx+1:1}'
                    if(verbose):
                        print('    reading:',name)
                    c = p.readArray(name)
                    c.tofile(ofp)
                    c = None
            elif n == 'cell_daughter':
                daughter.tofile(ofp)
            elif n == 'cell_index':
                nbrs.tofile(ofp)
                nbrs = None
            elif n == 'cell_level':
                cell_level.tofile(ofp)
            elif n == 'face_type':
                face_type.tofile(ofp)
            else:
                # read from file
                name = f'{n}_0'
                if(verbose):
                    print('    reading:',name)
                c = p.readArray(name)
                c.tofile(ofp)
                c = None

        # Write variable data
        trailer = "\nContents: eap-patterns-bin, 1(i64), ndim(i64), ncell(i64), name_len(i64), nVars(i64), " + \
            "list of names + element_size(i64) + offsets(i64), data(doubles)\n"
        ofp.write(trailer.encode())

        if(verbose):
            print('done writing')
    counts = [0]*6
    for i2,t in enumerate(face_type):
        for iCell,s in enumerate(t):
            if s == 1 and i2%2 == 1:
                print('HIgh side 1:', iCell, t)
            counts[s] += 1
    print(' typecount=', sum(counts), counts)
    print('nFacecount=', sum(sum(nFaces)), [sum(nFaces[:,iType]) for iType in range(1,6)])
    for iDim in range(p.ndim):
        print('idim=', iDim, nFaces[idim,:])

    icell = 53
    # print(f'cell={icell-1} faces={face_type[:][icell-1]}')
if __name__ == "__main__":
    import sys

    myFile = sys.argv[1]
    if len(sys.argv) > 2:
        outfile = sys.argv[2]
    else:
        outfile = "out.bin"

    convertPIO(myFile, outfile, True)
