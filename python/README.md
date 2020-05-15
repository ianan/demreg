# demreg python
This is a test version of demreg in python - not currently recommended for public use. 

Development is by Alasdair Wilson and the most up-to-date version is available at [https://github.com/alasdairwilson/demreg-py](https://github.com/alasdairwilson/demreg-py).

This is a python version of the "mapping" demreg code which you can find in this repo in ../idl/ 

##Longer term things

If you use this code (once its recommended for public use) please reference back to the original [Hannah & Kontar A&A 539 A146 (2012)](https://doi.org/10.1051/0004-6361/201117576).

##### Potted history

* The original SSWIDL code to obtain the Differential Emission Measure (DEM) from solar data using regularised inversion in SSWIDL, DEM Reg, is the one featured in [Hannah & Kontar (2012)](https://doi.org/10.1051/0004-6361/201117576) and in the repo in ../idl_org/
* This was based on the regularization code of [Kontar et al (2004)](https://doi.org/10.1007/s11207-004-4140-x) which is for X-ray spectrum inversion. It implements [Hansen (1992)](https://doi.org/10.1088/0266-5611/8/6/005) GSVD version of [Tikhonov (1963)](https://scholar.google.com/scholar_lookup?author=Tikhonov%2C+A.+N.&journal=Soviet+Math.+Dokl.&volume=4&pages=1035&publication_year=1963) regularization.
* A faster/parallized version of [demreg](http://www.astro.gla.ac.uk/~iain/demreg/map/) was developed for working with SDO/AIA images, not just single pixel/data sets. This "mapping" version of the code was featured in [Hannah & Kontar (2013)](https://doi.org/10.1051/0004-6361/201219727). An updated version of this code, with bug fixes, and to work with a variety of data (not just SDO/AIA) is the version now available in this repo in ../idl/ (however it is not currently parallized - no IDL Bridges) This python version is based on this idl version.
* A previous non-public python version of the original DEM Reg was developed and used by Erwin Verwichte and Petra Kohutova in [Kohutova & Verwichte (2016)](https://doi.org/10.3847/0004-637X/827/1/39) and [Verwichte & Kohutova (2017)](https://doi.org/10.1051/0004-6361/201730675).
* There might be other versions out there and not sure of their current state of features: [dstansby/demregpy/](https://github.com/dstansby/demregpy), ... ? 

