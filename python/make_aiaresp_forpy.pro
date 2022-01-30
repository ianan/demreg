pro make_aiaresp_forpy

  ; IDL script to produce the AIA temeprature response functions,
  ; then output it is a simple format that python will be happy with.
  ;
  ; This is the v10 responses with no degradation correction in them, since ran in Jan 2022
  ; https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/response/V10_release_notes.txt
  ;
  ; So need to degradation correct your AIA data (using ssw or aiapy)
  ; before combining with these responses
  ;
  ; 17-Aug-2020 IGH
  ; 25-Feb-2021 IGH - Renamed output file as not necessarily v9
  ; 30-Jan-2022 IGH - Added in noblend output as well
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ;  tresp=aia_get_response(/temperature,/dn,/eve,timedepend_date='01-Jul-2010')
  tresp=aia_get_response(/temperature,/dn,/eve)
  ids=[0,1,2,3,4,6]
  channels=tresp.channels[ids]
  logt=tresp.logte
  tr=tresp.all[*,ids]
  units=tresp.units
  save,file='aia_tresp_en.dat',channels,logt,tr,units
   
  ; Do the noblend version so don't include crosstalk between 131+335 and 94+304
  tresp=aia_get_response(/temperature,/dn,/eve,/noblend)
    ids=[0,1,2,3,4,6]
  channels=tresp.channels[ids]
  logt=tresp.logte
  tr=tresp.all[*,ids]
  units=tresp.units
  save,file='aia_tresp_en_nb.dat',channels,logt,tr,units

end