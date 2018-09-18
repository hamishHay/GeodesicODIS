# GeodesicODIS
GeodesicODIS is a finite volume fluid dynamics code for simulating global oceanic 
flow in icy satellites and other bodies. GeodesicODIS currently solves the Laplace 
Tidal Equations (LTE) assuming a 1-layer free-surface or subsurface ocean of uniform thickness. 

## GeodesicODIS vs ODIS
GeodesicODIS is a new version of the finite difference code ODIS. ODIS solved the LTE 
over a rectangular structured lat-lon mesh ([Hay and Matsuyama, 2017](https://www.sciencedirect.com/science/article/pii/S0019103516300239)), while GeodesicODIS does the same over an 
unstructured mesh based on the icosahedral-geodesic grid ([Hay and Matsuyama, 2018](https://www.sciencedirect.com/science/article/pii/S0019103518304470?via%3Dihub#!)).

## Getting Started
To get started with ODIS, copy the repository to your local machine and run make.

### Prerequisites
Three prerequisites are required for ODIS: the Fortran-95 SHTOOLS library, Intel's 
c++ compiler with MKL, and the HDF5 library.

[SHTools](https://shtools.oca.eu/shtools/) is a high-performance spherical harmonics library used in ODIS to decompose the tidal displacement at the ocean surface into spherical harmonics, which is required for subsurface ocean calculations and free-surface calculations that include ocean self-gravity.

The [Intel c++ compiler](https://software.intel.com/en-us/c-compilers) is a highly optimized compiler for numerical applications written in c++ on Intel processors. We make use of the compiler and the Intel Math Kernel Library (MKL) to perform calculations involving sparse matrices. In ODIS, interpolation from the geodesic to lat-lon grid is performed with a single sparse matrix multiplication using MKL.

[HDF5](https://www.hdfgroup.org/) is an optimized file format for storing large amounts of numerical data. ODIS makes use of this format to write *all* of its numerical output to memory.

## Running ODIS
ODIS can be run after running 'make'. Just run "./ODIS" with the default input file 
to begin subsurface ocean simulations on Enceladus.

## Authors
* **Hamish Hay** - *Lead Developer* [Lunar and Planetary 
Lab](https://www.lpl.arizona.edu/graduate/students/hamish-hay)

## Version History
GeodesicODIS 1.0: Recorded in branch HayMatsuyama2018.
