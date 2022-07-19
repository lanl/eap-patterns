# EAP patterns
## Authors
- [Sriram Swaminarayan](mailto:sriram@lanl.gov)

Memory access and iteration patterns from the EAP code base with the physics removed.

This is intended to be a serial app representing memory access
patterns.  We will  will populate the data structures from EAP output files
to provide representative patterns of face and cell loops within the
EAP code base.

The arguments to the program are the EAP output files and the number
of MPI processors to emulate.  Currently there is no MPI in the
application and the number of processors represents a means of
emulating the halo (clone) cells around the domain of the given
processor.  An optional argument `processor_ID` can be provided to 
specify the processor to emulate.

In the future we plan to include the ability to run multiple MPI
ranks, not for parallel performance, but to more accurately represent
the on-node demands on memory bandwidth. 

## Quick Start: Gradients emulation

- `mkdir build`
- `cd build`
- `cmake ..`
- `make -j`
- `./mygrad <pio-dump-file> <nProcs> [processor_ID]` 

Currently this code does nothing except compile without errors.  

## Next Steps:
- [x] Populate mesh from PIO file
- [ ] Calculate derivatives for a sample variable
- [ ] Compare with results from EAP application 
- [ ] Create a HDF5 specification
- [ ] Populate mesh from HDF5
- [ ] Test HDF5
- [ ] Create a binary file specification
- [ ] Populate mesh from binary file 
- [ ] Test binary
- [ ] Get EAP blessings
- [ ] Release open source
 
***



