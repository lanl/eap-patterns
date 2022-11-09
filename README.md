# EAP patterns
## Authors
- [Sriram Swaminarayan](mailto:sriram@lanl.gov)

Memory access and iteration patterns from the EAP code base with the physics removed.

This is intended to be a serial app representing memory access
patterns.  We will  will populate the data structures from EAP output files
to provide representative patterns of face and cell loops within the
EAP code base.

The argument to the program is an EAP output files converted
using the python/convertPIO.py script in this repository.  

In the future we plan to include the ability to run openMP
ranks, not for parallel performance, but to more accurately represent
the on-node demands on memory bandwidth. 

## Quick Start: Gradients emulation
### Serial
- `mkdir build`
- `cd build`
- `cmake .. -DEP_MPI=OFF`
- `make -j`
- `./mygrad <converted_PIO_file>`

### Parallel
- `mkdir build_parallel`
- `cd build_parallel`
- `cmake .. -DEP_MPI=ON`
- `make -j`
- `mpirun -n 2048 ./mygrad <converted_PIO_file>`

 
***



