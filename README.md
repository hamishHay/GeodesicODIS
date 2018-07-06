# GeodesicODIS
GeodesicODIS is a finite volume fluid dynamics code for simulating global oceanic 
flow in icy satellites and other bodies. GeodesicODIS currently solves the Laplace 
Tidal Equations (LTE) assuming a 1-layer free-surface or subsurface ocean in a 
spatially uniform ocean. CAUTION: GeodesicODIS is still under development and may be 
unstable/unreliable.
## GeodesicODIS vs ODIS
GeodesicODIS is a new version of the finite difference code ODIS. ODIS solved the LTE 
over a rectangular structured lat-lon mesh, while GeodesicODIS does the same over an 
unstructured mesh based on the icosahedral-geodesic grid.
## Getting Started
To get started with ODIS, copy the repository to your local machine and run make.
### Prerequisites
Two prerequisites are required for ODIS: the Fortran-95 SHTOOLS library and Intel's 
c++ compiler with MKL.
## Running ODIS
ODIS can be run after running 'make'. Just run "./ODIS" with the default input file 
to begin subsurface ocean simulations on Enceladus.
## Authors
* **Hamish Hay** - *Lead Developer* [Lunar and Planetary 
Lab](https://www.lpl.arizona.edu/graduate/students/hamish-hay)
