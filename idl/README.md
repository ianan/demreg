# demreg

This directory contains the latest version of the regularized inversion DEM code. 

It is based on the mapping version of the DEM reg code (http://www.astro.gla.ac.uk/~iain/demreg/map/) - it works but is still occassionally being worked on/tweaked.

The improvements over the older code are:
* Works with any instrument - just supply the data and corresponding temperature response functions
* Can work data of various dimensions - a single pixel/data set or an array/map of pixels
* Some bug fixes to remove outlier cases where NaNs are returned

Code:
* demmap_pos.pro - Main DEM Reg code (input is spatial 1D of data)
* dem_inv_gsvdcsq.pro - Routine used by demmap_pos.pro
* dem_inv_reg_parameter_map.pro - Routine used by demmap_pos.pro
* dn2dem_pos_nb.pro - Wrapper to convert data (single set or map) for use with demmap_pos.pro
* example_*.pro - Examples using a variety of synthetic or real data
* dem2dn.pro - Used by example codes to synthesize DN from model DEMs
* aia_resp* - IDL save files of AIA temperature responses, used by example codes
