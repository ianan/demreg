pro make_aiaresp_forpy

; IDL script to produce the AIA temeprature response functions for the chosen date,
; then output it is a simple format that python will be happy with.
;
; 15-May-2020 IGH
;
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  tresp=aia_get_response(/temperature,/dn,/eve,timedepend_date='01-Jul-2010')
  
  ids=[0,1,2,3,4,6]
  channels=tresp.channels[ids]
  logt=tresp.logte
  tr=tresp.all[*,ids]
  units=tresp.units
  
  save,file='aia_resp.dat',channels,logt,tr,units


stop
end