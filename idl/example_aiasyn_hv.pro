pro example_aiasyn_hv,day=day,hour=hour,min=mins,$
  min_snr=min_snr,sat_lvl=sat_lvl,serr_per=serr_per,$
  tempdir=tempdir,respdir=respdir, use_diw=use_diw,em=em

  ; Example to produce full disk DEM map from synoptic AIA data for a given day, hour, mins
  ;
  ; Note that this code will go and get the synoptic date from jsoc and
  ; save it in tempdir (default is a temp/ in current working directory)
  ;
  ; Note that this code will calculate AIA temperature response functions for the month/year of the obs
  ; and save it in respdir (default is a resp/ in current working directory)
  ;
  ; Options to:
  ; min_snr         use pixels with this min SNR (defulat 3)
  ; sat_lvl         use pixels with DN/px less than this value (default 15,000)
  ;                     Note that streaks from a saturated region won't be removed unless >sat_lvl
  ; use_diw         use an initial DEM weighting at start of calc (defaul 0, no)
  ; serr_per        add a systematic uncertainty (%) to the data (default is 0%)
  ; em              final maps in cm^-5 instead of cm^-5 K^-1
  ;
  ; Todos:
  ;   ***   If badly saturated data or noisy (very small ind.exptime) then automatically use neighbouring time instead?
  ;   ***   Re-bridge the DEM code to speed things up again ?
  ;   ***   Best initial DEM weighting to use (if any ??)
  ;   ***   Remove DEMs with large uncertainties?
  ;   ***   Average/interpolate the temperature bins?
  ;   ***   Output as jpeg2000s
  ;
  ; May 2016      IGH   Started
  ; 17-Nov-2016   IGH   Tidied up code and added comments
  ; 18-Nov-2016   IGH   Added option for EM output via /em
  ;
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;  ;   About X3, Impulsive start of X-flare
;    if (n_elements(day) ne 1) then day='24-Oct-2014'
;    if (n_elements(hour) ne 1) then hour='21'
;    if (n_elements(mins) ne 1) then mins='26'

  ; About B7, few  ARs, day of C-flares
    if (n_elements(day) ne 1) then day='25-Dec-2015'
    if (n_elements(hour) ne 1) then hour='12'
    if (n_elements(mins) ne 1) then mins='00'

  ;  ; About A5, quiet day
  ;  if (n_elements(day) ne 1) then day='01-Jul-2010'
  ;  if (n_elements(hour) ne 1) then hour='12'
  ;  if (n_elements(mins) ne 1) then mins='00'

  date=day+' '+hour+':'+mins+':00'

  ; Only calculate the DEM for pixels which has this minimum SNR
  if (n_elements(min_snr) lt 1) then min_snr=3.0
  ; Ignore data with values above this level
  ; This is the default AIA saturation of >15,000 DN/px
  if (n_elements(sat_lvl) lt 1) then sat_lvl=1.5e4

  ;Systematic uncertainty (in % terms) to add to the data
  if (n_elements(serr_per) lt 1) then serr_per=0.0

  ; Where to store the synoptic images downloaded from jsoc?
  if (n_elements(temp_dir) ne 1) then temp_dir=curdir()+'/temp/'

  ; Where to store the calculated temperature response functions?
  if (n_elements(resp_dir) ne 1) then resp_dir=curdir()+'/resp/'

  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; For this can we use load in the 1024x1024 synoptic from jsoc (this ok??)
  ystr=strmid(anytim(date,/ccsds),0,4)
  mstr=strmid(anytim(date,/ccsds),5,2)
  dstr=strmid(anytim(date,/ccsds),8,2)
  durl='http://jsoc.stanford.edu/data/aia/synoptic/'+ystr+'/'+mstr+'/'+dstr+'/H'+hour+'00/'
  wavenum=['94','131','171','193','211','335']
  ffs='AIA'+ystr+mstr+dstr+'_'+hour+mins+'_0'+['094','131','171','193','211','335']+'.fits'
  furls=durl+ffs

  ; Already have the files?
  ftest=file_test(temp_dir+ffs)
  ; If not go and download them
  if (product(ftest) ne 1) then wgetc=ssw_wget_mirror(furls,temp_dir,/spawn)

  ;Read in the data
  ; Don't need to aia_prep as alrady 1.5 (?)
  read_sdo,temp_dir+ffs,ind,data,/uncomp_delete

  ; Just double check in the correct order of 94, 131, 171, 193, 211, 335
  wv_srt=sort(ind.wavelnth)
  ind=ind[wv_srt]
  data=data[*,*,wv_srt]

  nx=n_elements(data[*,0,0])
  ny=n_elements(data[0,*,0])
  nf=n_elements(data[0,0,*])

  xc=ind[0].xcen
  yc=ind[0].ycen
  dx=ind[0].cdelt1
  dy=ind[0].cdelt2

  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Filter out the "bad" pixels
  ; Only need to set the actual bad ones to 0 (not all filters at that location)
  ; as the DEM code will only perform the calc if all filters are non-zero at that pixel location

  ; If the value is above sat_lvl get rid of it
  ; assuming bad from saturation
  id=where(data ge sat_lvl,nid)
  if (nid gt 1) then data[id]=0.0

  ; Anything less than 0 just set to 0
  id=where(data le 0,nid)
  if (nid gt 1) then data[id]=0.0

  ; Set the edge few pixels to 0 to minimise problems - from prep'd data that has been shifted (?)
  edg0=10
  data[0:edg0-1,*,*]=0.0
  data[*,0:edg0-1,*]=0.0
  data[nx-edg0:nx-1,*,*]=0.0
  data[*,nx-edg0:nx-1,*]=0.0

  ; Work out the errors
  ; The synoptic data is in DN/px
  edata=fltarr(nx,ny,nf)
  ; Proper way but a bit slow?
  ;  for i=0,nf-1 do edata[*,*,i]=aia_bp_estimate_error(reform(data[*,*,i]),replicate(wavenum[i],nx,ny),n_sample=16)
  ; So instead quick approx based on aia_bp_estimate_error.pro
  npix=4096.^2/(nx*ny)
  edata=fltarr(nx,ny,nf)
  gains=[18.3,17.6,17.7,18.3,18.3,17.6]
  dn2ph=gains*[94,131,171,193,211,335]/3397.
  rdnse=1.15*sqrt(npix)/npix
  drknse=0.17
  qntnse=0.288819*sqrt(npix)/npix
  ; error in DN/px
  for i=0, nf-1 do begin
    etemp=sqrt(rdnse^2.+drknse^2.+qntnse^2.+(dn2ph[i]*abs(data[*,*,i]))/(npix*dn2ph[i]^2))
    esys=serr_per*data[*,*,i]/100.
    edata[*,*,i]=sqrt(etemp^2. + esys^2.)
  endfor

  ; Get rid of data with too large an uncertaintity
  id=where(data/edata le min_snr,nid)
  if (nid gt 1) then data[id]=0.0

  ; For the DEM code need the data in DN/px/s
  durs=ind.exptime
  for i=0, nf-1 do data[*,*,i]=data[*,*,i]/durs[i]
  for i=0, nf-1 do edata[*,*,i]=edata[*,*,i]/durs[i]

  ; What temperature binning do you want for the DEM?
  ; temps variable are the bin edges
  ; These are the bin edges
  ;  temps=[0.5,1,1.5,2,3,4,6,8,11,14,19,25,32]*1e6
  ;  temps=[0.5,1,2,4,6,8,11,14,19]*1d6
  ;  logtemps=alog10(temps)

  ;  ;or more bins and then rebin at the end?
  logtemps=5.7+findgen(17)*0.1
  temps=10d^logtemps

  ; This is is the temperature bin mid-points
  mlogt=get_edges(logtemps,/mean)
  nt=n_elements(mlogt)

  ; Given that the responses are time dependent, use one on a monthly basis
  ; Calculate and save if not already done so, then load in
  respfile='resp/aia_resp_'+mstr+ystr+'.dat'
  ; Need to make the response functions?
  if (file_test(respfile) eq 0) then begin
    tresp=aia_get_response(/temperature,/dn,/chianti,/noblend,/evenorm,timedepend_date=ystr+'/'+mstr+'/01')
    save,file=respfile,tresp
  endif
  restore,file=respfile

  ; Only want the coronal ones without 304A
  idc=[0,1,2,3,4,6]

  tr_logt=tresp.logte
  ; Don't need the response outside of the T range we want for the DEM
  gdt=where(tr_logt ge min(logtemps) and tr_logt le max(logtemps),ngd)
  tr_logt=tr_logt[gdt]
  TRmatrix=tresp.all[*,idc]
  TRmatrix=TRmatrix[gdt,*]

  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Use the initial weighting to the DEM calc or not (default not)
  if keyword_set(use_diw) then begin

    dem_norm0=dblarr(nx,ny,nt)
    ; Not sure best form for this
    ; this example is just a gaussian then smoothed
    ; absoulte value doesn't matter (within reason) just the shape
    dem_norm_temp=exp(-(findgen(nt)+1-nt*0.5)^2/12)
    dem_norm_temp=smooth(dem_norm_temp,3)

    ; Very rough initial normalisation but just testing if the code works
    for xx=0,nx-1 do begin
      for yy=0,ny-1 do begin
        dem_norm0[xx,yy,*]=dem_norm_temp
      endfor
      ; Do DEM calculation
      dn2dem_pos_nb, data,edata,TRmatrix,tr_logt,temps,dem,edem,elogt,chisq,dn_reg,/timed,dem_norm0=dem_norm0
    endfor
  endif else begin
    ; Do DEM calculation without the specified intial weighting
    dn2dem_pos_nb, data,edata,TRmatrix,tr_logt,temps,dem,edem,elogt,chisq,dn_reg,/timed
  endelse

  ; Convert the DEM output into maps
  dem_maps=replicate(make_map(fltarr(nx,ny),xc=xc,yc=yc,dx=dx,dy=dy,date=date),nt)

  ; Want the final map to be in cm^-5 (EM) or cm^-5 K^-1 (DEM)
  if keyword_set(em) then begin
    dt=fltarr(nt)
    for i=0,nt-1 do dt[i]=10d^(logtemps[i+1])-10d^(logtemps[i])
    em=dblarr(nx,ny,nt)
    for xx=0,nx-1 do begin
      for yy=0,ny-1 do begin
        em[xx,yy,*]=dem[xx,yy,*]*dt
      endfor
    endfor
    dem_maps.data=em
    idd='EM '
  endif else begin
    dem_maps.data=dem
    idd='DEM '
  endelse

  temp_ids=strarr(nt)
  for i=0,nt-1 do temp_ids[i]=string(logtemps[i],format='(f3.1)')+' - '+string(logtemps[i+1],format='(f3.1)')
  dem_maps.id=idd+'logT: '+temp_ids

  ; plot them - assuming using default 16 temp bins and no interpolating
  !p.multi=[0,4,4]
  loadct,5,/silent
  
  if keyword_set(em) then drang=[1d25,1d30] else drang=[1d19,1d23]
  
  for i=0,nt-1 do plot_map,dem_maps[i],title=dem_maps[i].id,/log,dmin=drang[0],dmax=drang[1],chars=2.0
  xyouts,10,10,date,chars=1.2,/device

  stop

  ; best to avoid file deletion at the moment as still work in progress
  ;  ssw_file_delete,temp_files

  stop
end