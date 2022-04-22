pro example_dem_pxl

  ; Example script to recover the DEM from AIA single pixel data
  ; The AIA data is synthetic for specified Gaussian DEM model
  ;
  ;
  ; 13-Apr-2016 IGH
  ; 27-Apr-2016 IGH   - Changed the naming of the temperatures to make things clearer:
  ;                     tr_logt is the binning of the response function
  ;                     temps is the bin edges you want for the DEM
  ;                     logtemps is the log of the above
  ;                     mlogt is the mid_point of the above bins
  ; 20-May-2019 IGH   - Minor update: must use the timedepend_date option when getting aia response  
  ; 15-Aug-2019 IGH   - Renamed from _1D to _pxl to better reflect what this example code does 
  ;                     Changed file name *_demmap_* to *_map_*
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ; Gaussian DEM model parameters
  d1=4d22
  m1=6.5
  s1=0.15

  ; What temperature binning do you want of the DEM?
  ; These are the bin edges
  tt=5.7+findgen(30)/20.
  temps=10d^tt;[0.5,1,1.5,2,3,4,6,8,11,14,19,25,32]*1e6
  logtemps=alog10(temps)
  ; This is is the temperature bin mid-points
  mlogt=get_edges(logtemps,/mean)

  ; Need to make the response functions?
  if (file_test('aia_resp.dat') eq 0) then begin
    tresp=aia_get_response(/temperature,/dn,/chianti,/noblend,/evenorm,timedepend_date='01-Jul-2010')
    save,file='aia_resp.dat',tresp
  endif else begin
    restore,file='aia_resp.dat'
  endelse
  
  ; Only want the coronal ones without 304A
  idc=[0,1,2,3,4,6]

  tr_logt=tresp.logte
  ; Don't need the response outside of the T range we want for the DEM
  gdt=where(tr_logt ge min(logtemps) and tr_logt le max(logtemps),ngd)
  tr_logt=tr_logt[gdt]
  TRmatrix=tresp.all[*,idc]
  TRmatrix=TRmatrix[gdt,*]
  root2pi=sqrt(2.*!PI)
  dem_mod=(d1/(root2pi*s1))*exp(-(tr_logT-m1)^2/(2*s1^2))
  dn_mod=dem2dn(tr_logT, dem_mod, TRmatrix)

  dn=dn_mod.dn
  nf=n_elements(dn)
  ; workout the error on the data
  edn_in=fltarr(nf)
  gains=[18.3,17.6,17.7,18.3,18.3,17.6]
  dn2ph=gains*[94,131,171,193,211,335]/3397.
  rdnse=[1.14,1.18,1.15,1.20,1.20,1.18]
  ; assume all obs were 2.9s long
  dn0=dn*2.9
  shotnoise=sqrt(dn2ph*dn0)/dn2ph/2.9
  ; error in DN/s/px
  edn=sqrt(rdnse^2+shotnoise^2)
  
  dn2dem_pos_nb, dn, edn,TRmatrix,tr_logt,temps,dem,edem,elogt,chisq,dn_reg,/timed;,/gloci,glcindx=[0,1,1,1,1,1]

  yr=d1*[5e-5,1e1]
  !p.thick=2
  plot,tr_logt,dem_mod,/ylog,chars=2,xtit='Log!D10!N T', ytit='DEM [cm!U-5!N K!U-1!N]',yrange=yr,ystyle=17,xstyle=17
  loadct,39,/silent
  for i=0,n_elements(mlogt)-1 do oplot,logtemps[i:i+1],dem[i]*[1,1],color=250
  demax=(dem+edem) < yr[1]
  demin=(dem-edem) > yr[0]
  for i=0,n_elements(mlogt)-1 do oplot, mlogt[i]*[1,1],[demin[i],demax[i]],color=250

  stop
end
