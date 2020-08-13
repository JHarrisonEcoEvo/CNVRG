## Test environments
* local OS X install, R 3.1.2
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran
* Debian Linux, R-devel, GCC ASAN/UBSAN

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

 GNU make is a SystemRequirements.

I use GNU make because it is used by some of the dependencies of this package. Regular make throws an error. 

