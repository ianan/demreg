Note that this is the original version of the DEM reg code that does one data set 
(whatever instruments, as long as you provide the response/contribution functions) 
and has many options. Best for a single/few data sets - not maps.

This code performs a regularized inversion (Tikhonov 1963) on multi-filter
and/or spectroscopic line intensities with the corresponding temperature
response/contribution functions to calculate the Differential Emission Measure
(DEM), as detailed in Hannah & Kontar A&A 2011. The is an modification of the
regularized inversion code already included in the xray ssw package (Kontar et
al 2004 Sol Phys) under $SSW/packages/xray/idl/inversion/.

18-Nov-2011 IGH

Update: 08-Nov-2015 IGH
Note that the DEM code now does not require SSW to run and neither does most of the example
scripts (i.e. line_example.pro, aia_example.pro and aia_example_new.pro)
Update: 19-Geb-2020 IGH 
EMD version now available - calculation done in EMD [cm^5] instead of DEM [cm^5 K^-1]

##############################################

The main driver programs is:
data2dem_reg.pro
or
data2em_reg.pro if want EMD solution

##############################################

The regularization is performed in various steps using dem_inv_*
dem_inv_make_constraint.pro	
	 ** Calculates the chosen constraint matrix
dem_inv_gsvdcsq.pro			
	  **Performs GSVD
dem_inv_reg_parameter.pro		
	  ** Calculate the regularization parameter
dem_inv_reg_parameter_pos.pro		
	  ** Additional version that calculates reg param and solution for positive case
dem_inv_reg_solution.pro		
	  ** Calculates the regularized solution
dem_inv_reg_resolution.pro		
	  ** Calculates the DEM horizontal error (temperature resolution)
dem_inv_confidence_interval.pro	
	** Calculates the DEM vertical error	

##############################################

Example scripts are given for SDO/AIA and Hinode/EIS for a variety of DEMs.
These are batch scripts so to execute just via
IDL>@aia_example

aia_example.pro		
	** SDO/AIA with a Gaussian DEM Model
aia_example_ar.pro		
	** SDO/AIA with the CHINATI Active Region DEM Model
aia_example_new.pro		
	** SDO/AIA with generic Gaussian DEM Model
line_example.pro		
	 ** Hinode/EIS with the CHINATI Quiet Sun DEM Model

For the CHIANTI model DEM examples both the SDO/AIA and chianti ssw packages need to
be installed, i.e. setenv SSW_INSTR "aia chianti"

##############################################

The code is distributed under a Creative Commons through the
Attribution-Noncommercial-Share Alike 3.0 license (can copy, distribute and
adapt the work but full attribution must be given and cannot be used for
commercial purposes)
