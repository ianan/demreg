pro example_demmap_1d_fe18

  ; Example script to recover the DEM from AIA single pixel data
  ; The AIA data is synthetic for specified Gaussian DEM model
  ;
  ; Also uses the Fe18 pseudo-channel via Del Zanna 2013 A&A [94Å, 171Å, 211Å]
  ;     - it might help as an addition (more with actual data than synthetic?
  ; 
  ; 13-Apr-2016 IGH
  ; 27-Apr-2016 IGH   - Changed the naming of the temperatures to make things clearer:
  ;                     tr_logt is the binning of the response function
  ;                     temps is the bin edges you want for the DEM
  ;                     logtemps is the log of the above
  ;                     mlogt is the mid_point of the above bins
  ; 20-May-2019 IGH   - Minor update: must use the timedepend_date option when getting aia response   
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ; Gaussian DEM model parameters
  d1=4d22
  m1=6.4
  s1=0.1

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
  dem_mod=(d1/(root2pi*s1))*exp(-(tr_logt-m1)^2/(2*s1^2))
  dn_mod=dem2dn(tr_logt, dem_mod, TRmatrix)

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

  ; Calculate the Fe XVIII
  ; Using empirical params of del Zanna and filter combo of del Zanna i.e
  ; Del Zanna 2013 A&A [94Å, 171Å, 211Å]
  ; assuming all have the same duration (i.e. can use DN/s instead of DN)
  dn_fe18z=(dn[0]-dn[4]/120.-dn[2]/450.) >0.
  dn_all18=[dn,dn_fe18z]
  edn_all18=sqrt(edn[0]^2+(edn[4]/120.)^2+(edn[2]/450.)^2)
  TRfe18z= (TRmatrix[*,0]-TRmatrix[*,4]/120.-TRmatrix[*,2]/450.) >0.
  ; Make sure nothing in the low T peak
  ; A fairly crude way of doing this.....
  TRfe18z[where(tr_logt le 6.4)]=0.
  TRall18=dblarr(n_elements(TRmatrix[*,0]),7)
  TRall18[*,0:5]=TRmatrix
  TRall18[*,6]=TRfe18z

  dn2dem_pos_nb, dn_all18, edn_all18,TRall18,tr_logt,temps,dem18,edem18,elogt18,chisq18,dn_reg18,/timed

  yr=d1*[1e-3,1e1]
  plot,tr_logt,dem_mod,/ylog,chars=2,xtit='Log T', ytit='DEM [cm!U-5!NK!U-1!N]',yrange=yr,ystyle=17,xstyle=17
  loadct,39,/silent
  for i=0,n_elements(mlogt)-1 do oplot,logtemps[i:i+1],dem[i]*[1,1],color=250
  demax=(dem+edem) < yr[1]
  demin=(dem-edem) > yr[0]
  for i=0,n_elements(mlogt)-1 do oplot, mlogt[i]*[1,1],[demin[i],demax[i]],color=250

  for i=0,n_elements(mlogt)-1 do oplot,logtemps[i:i+1],dem18[i]*[1,1],color=120
  demax18=(dem18+edem18) < yr[1]
  demin18=(dem18-edem18) > yr[0]
  for i=0,n_elements(mlogt)-1 do oplot, mlogt[i]*[1,1],[demin18[i],demax18[i]],color=120

  xyouts, tr_logt[0]+0.1,0.5*yr[1],'AIA 6',color=250,/data,chars=2
  xyouts, tr_logt[0]+0.1,0.2*yr[1],'AIA 6 + Fe18',color=120,/data,chars=2

  stop
end
