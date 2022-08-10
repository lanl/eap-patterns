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

Takes two arguments: PIO file and numprocs

Use the PIO class to estimate number of
clone cells for a simulation

Uses a slow but reliable algorithm.

** REQUIRES Python3 **

Author: Sriram Swaminarayan (sriram@lanl.gov)
Date: March 24, 2021
Version: 1.0

"""
import sys, os
import numpy as np
basepath = os.path.dirname(os.path.abspath(__file__))
sys.path.append(basepath + "/..")
from pio import pio
import re
import struct


def fixNbrs(fname, outfile="out-dmp000000", verbose=False):

    myNames = re.compile("cell_center_|cell_daughter_|cell_index_|cell_level_|vcell_")
    p = pio(fname)

    daughter = p.readArray("cell_daughter_0").astype(np.int64)

    print(len(daughter),daughter)
    nbrs = [None] * (2 * p.ndim)
    for idx in range(2 * p.ndim):
        print(f"cell_index_{idx+1}")
        a = p.readArray(f"cell_index_{idx+1}").astype(np.int64)
        for iCell in range(p.numcell):
            if daughter[a[iCell]-1] > 0:
                a[iCell] = daughter[a[iCell]-1]

        nbrs[idx] = a.astype(np.double)
    daughter = None
    newNames = []
    
    for n in p.names:
        m = myNames.match(n)
        if m is not None:
            newNames.append(n)

    # Now to write the file
    with open(outfile,'wb') as ofp:
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
        print('   writing offsets')
        for n in newNames:
            ofp.write(f"{n:<{p.lName}}".encode())
            ofp.write(struct.pack("q",offset))
            offset += 8 * p.numcell

        # Write data
        for n in newNames:
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

        print('done writing')
if __name__ == "__main__":
    import sys

    myFile = sys.argv[1]

    fixNbrs(myFile)
