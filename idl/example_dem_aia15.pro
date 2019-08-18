pro example_dem_aia15

  ; Example script to recover the DEM from AIA Lvl1.5 fits files
  ; The specific AIA fits used here are not include with the code
  ;
  ; 14-Apr-2016 IGH
  ; 27-Apr-2016 IGH   - Changed the naming of the temperatures to make things clearer:
  ;                     tr_logt is the binning of the response function
  ;                     temps is the bin edges you want for the DEM
  ;                     logtemps is the log of the above
  ;                     mlogt is the mid_point of the above bins
  ; 28-Apr-2016       - Still testing: not optimised T bins, initial weighting or errors
  ; 20-May-2019 IGH   - Minor update: must use the timedepend_date option when getting aia response
  ; 15-Aug-2019 IGH   - Changed file name *_demmap_* to *_map_*
  ; 18-Aug-2019 IGH   - example_get_aia15.pro now used to create input data
  ;                     Now accounts for number of original 0.6" pixels when calculating noise
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ; Need to have a fits file of the 6 AIA maps properly prep'd and rebinned (if needed)
  ; This one is produved by example_get_aia15.pro
  out_dir='/Users/iain/Desktop/temp/'
  fits2map,out_dir+'test_maps_rebin_20151225.fits',mm

  ; How many pixels combined to make this map?
  npix=(mm[0].dx/0.6)*(mm[0].dy/0.6)

  ; For this example just do part of the map
  ;  sub_map,mm,mms,xrange=[-650,-250],yrange=[-500,-200]
  ;  sub_map,mm,mms,xrange=[550,1000],yrange=[-450,-150]
  sub_map,mm,mms,xrange=[550,700],yrange=[-350,-200]
  mm=mms

  ; Setup the data for input to DEMREG code
  dn0=mm.data
  durs=mm.dur
  ; Get into DN/s/px
  nf=n_elements(durs)
  for i=0, nf-1 do dn0[*,*,i]=dn0[*,*,i]/durs[i]
  na=n_elements(dn0[*,0,0])
  nb=n_elements(dn0[0,*,0])

  ; Work out the errors on the data
  ; Ignoring systematic at the moment
  ; This can also be done (and better it do it?) via aia_bp_estimate_error.pro

  ; workout the error on the data
  edn0=fltarr(na,nb,nf)
  gains=[18.3,17.6,17.7,18.3,18.3,17.6]
  dn2ph=gains*[94,131,171,193,211,335]/3397.
  rdnse=[1.14,1.18,1.15,1.20,1.20,1.18]
  ; error in DN/s/px
  for i=0, nf-1 do begin
    shotnoise=sqrt(dn2ph[i]*abs(dn0[*,*,i])*durs[i]*npix)/dn2ph[i]
    edn0[*,*,i]=sqrt(rdnse[i]^2.+shotnoise^2.)/durs[i]/npix
  endfor

  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; What temperature binning do you want of the DEM?
  ; These are the bin edges
  ;  temps=[0.5,1,1.5,2,3,4,6,8,11,14,19,25,32]*1e6
  temps=[0.5,1,1.5,2,3,4,6,8,11,14,19]*1d6
  ;  temps=[0.5,1,2,4,6,9,14]*1d6
  logtemps=alog10(temps)
  ; This is is the temperature bin mid-points
  mlogt=get_edges(logtemps,/mean)
  nt=n_elements(mlogt)

  ; Need to make the response functions?
  if (file_test('aia_resp15.dat') eq 0) then begin
    tresp=aia_get_response(/temperature,/dn,/chianti,/noblend,/evenorm,timedepend_date=mm[0].time)
    save,file='aia_resp15.dat',tresp
  endif else begin
    restore,file='aia_resp15.dat'
  endelse

  ; Only want the coronal ones without 304A
  idc=[0,1,2,3,4,6]

  tr_logt=tresp.logte
  ; Don't need the response outside of the T range we want for the DEM
  gdt=where(tr_logt ge min(logtemps) and tr_logt le max(logtemps),ngd)
  tr_logt=tr_logt[gdt]
  TRmatrix=tresp.all[*,idc]
  TRmatrix=TRmatrix[gdt,*]
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  dn2dem_pos_nb, dn0, edn0,TRmatrix,tr_logt,temps,dem,edem,elogt,chisq,dn_reg,/timed

  loadct,39,/silent

  ;  ; Plot them all with the same scaling
  ;  ; needs ssw plot_image for this to work
  ;  !p.multi=[0,4,nt/3]
  ;  for t=0,nt-1 do plot_image,alog10(dem[*,*,t]),chars=2,max=23,min=19,$
  ;    title=string(temps[t]*1d-6,format='(f4.1)')+' to '+string(temps[t+1]*1d-6,format='(f4.1)')+' MK'


  ; Or plot as ssw maps of DEM
  mdem0=make_map(dem[*,*,0]*0,dx=mm[0].dx,dy=mm[0].dy,xc=mm[0].xc,yc=mm[0].yc,time=anytim(mm[0].time,/yoh,/trunc))
  mdem=replicate(mdem0,nt)
  for i=0,nt-1 do begin
    mdem[i].data=dem[*,*,0]
    ;    mdem[i].id=string(logtemps[i],format='(f4.2)')+' to '+string(logtemps[i+1],format='(f4.2)')+' Log!D10!N MK'
    mdem[i].id=string(temps[i]*1d-6,format='(f4.1)')+' to '+string(temps[i+1]*1d-6,format='(f4.1)')+' MK'
  endfor

  !p.multi=0
  plot_map,mdem[5],tit=mdem[5].id



  stop
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;
  ;  ; Rough initial normalisation really just for testing the code
  ;   ; norm0 needs to be same number of bins as output dem (and temps)
  ;   norm0=[1e-2,3e3,4.2e3,5e3,5e3,1e3,1e2,1e2,1e-2,1e-2]
  ;   ;;norm0=[1e-2,4.2e3,5e3,1e3,1e2,1e-2]
  ;
  ;   nxn=n_elements(dn0[*,0,0])
  ;   nyn=n_elements(dn0[0,*,0])
  ;   dem_norm0=dblarr(nxn,nyn,nt)
  ;   for xx=0,nxn-1 do begin
  ;     for yy=0,nyn-1 do begin
  ;       dem_norm0[xx,yy,*]=norm0
  ;     endfor
  ;   endfor
  ;
  ;  dn2dem_pos_nb, dn0, edn0,TRmatrix,tr_logt,temps,dem,edem,elogt,chisq,dn_reg,/timed,dem_norm0=dem_norm0
  ;
  ;  loadct,39,/silent
  ;  !p.multi=[0,4,nt/3]
  ;  ; Plot them all with the same scaling
  ;  ; needs ssw plot_image for this to work
  ;  for t=0,nt-1 do plot_image,alog10(dem[*,*,t]),chars=2,max=23,min=19,$
  ;    title=string(temps[t]*1d-6,format='(f4.1)')+' to '+string(temps[t+1]*1d-6,format='(f4.1)')+' MK'


  stop
end
