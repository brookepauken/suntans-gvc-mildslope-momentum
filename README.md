# GVC SUNTANS


Built upon SUNTANS (Fringer et al., 2006) with added GVC framework, hybrid vertical coordinate system, and a conservative momentum scheme that is stable for vanishingly small layer heights. The original SUNTANS user guide (see https://github.com/ofringer/suntans) is useful for getting started. Requires a grid generator, MPI, and ParMETIS (if run in parallel). 

**Quick start guide**
1. Download source code
2. Edit paths to mpi and ParMETIS in main/Makefile.in 
3. Navigate to main and type 'make' to compile
4. Navigate to example directory of choice
5. 'make clean all' to compile
6. 'make test' to run example script


This examples in this repository contain the code needed to run the examples in [this preprint](https://doi.org/10.48550/arXiv.2607.28356) submitted to Ocean Modelling. This README will be updated with the journal citation if published.


## Citations
O. B. Fringer, M. Gerritsen, and R. L. Street (2006), An unstructured-grid, finite-volume, nonhydrostatic, parallel coastal ocean simulator, Ocean Model., 14 (3-4), 139-173, doi:10.1016/j.ocemod.2006.03.006
