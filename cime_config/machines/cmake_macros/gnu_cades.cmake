string(APPEND FFLAGS " -O -fno-range-check")
set(HDF5_PATH "/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/hdf5-parallel/1.8.17/centos7.2_gnu5.3.0")
set(NETCDF_PATH "/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/netcdf-hdf5parallel/4.3.3.1/centos7.2_gnu5.3.0")
set(PNETCDF_PATH "/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/pnetcdf/1.9.0/centos7.2_gnu5.3.0")
set(LAPACK_LIBDIR "/software/tools/compilers/intel_2017/mkl/lib/intel64")
string(APPEND SLIBS " -L${NETCDF_PATH}/lib -Wl,-rpath=${NETCDF_PATH}/lib -lnetcdff -lnetcdf")
string(APPEND CXX_LIBS " -lstdc++")
set(MPICXX "mpic++")
set(SCXX "g++")