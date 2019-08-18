pro example_get_aia15, vso_get=vso_get, prep_aia=prep_aia

  ; Example script to get and prep some AIA lvl1.5 fits files
  ; Output of this code can be used with example_dem_aia15
  ; 
  ; Note - occasionally need to .r vso_get before compiling this code (dont know why)
  ;
  ; Option:
  ;   vso_get   Need to get some files (default no)
  ;   prep_aia  Need to prep the files you have got? (default no)
  ;
  ; 15-Aug-2019 IGH   - Started
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;

  ; A short time range to get the data from
  tt='2015-12-25 '+['12:00:00','12:00:11']
  out_dir='/Users/iain/Desktop/temp/'

  ; Deafult setup of what to do
  vso_get=0
  prep_aia=0
  
  ; Get some AIA data from vso
  if keyword_set(vso_get) then begin
    ; Do the vso search for sdo/aia over the time
    ; The wave search filter can only do min to max, so will recover 304 which do not want (optically thick)
    ss = vso_search(tt[0],tt[1], source='sdo',instr='aia',wave='94-335')
    ss = ss[where(ss.wave.min ne 304.)]

    ; Just double check we have one of each filter type
    ss=ss[uniq(ss.wave.min)]

    ; Get the data and put it somewhere
    status=vso_get(ss, out_dir=out_dir,/rice,sire='NSO')
  endif

  ; Prep the aia data and save out as a single fits
  if keyword_set(prep_aia) then begin
    ; Load in the fits and then prep them
    files=file_search(out_dir, '*.fits')
    read_sdo,files,ind,data

    aia_prep,ind,data,indp,datap,/use_hdr_pnt
    sid=sort(indp.wavelnth)

    index2map,indp,datap > 0.,maps
    maps=maps[sid]

    map2fits,maps,out_dir+'test_maps_20151225.fits'

  endif
  
  ; Load in the processed fits
  fits2map,out_dir+'test_maps_20151225.fits',maps
  
  wvid=strmid(maps.id,9,4)
;
;  ; Plot them
;  !p.multi=[0,3,2]
;  for i=0, 5 do begin
;    aia_lct,wave=wvid[i],/load
;    plot_map,maps[i],/log,charsize=1.5
;  endfor

  ; rebin the maps
  mapsr=rebin_map(maps,1024,1024)
  ; Plot them
  !p.multi=[0,3,2]
  for i=0, 5 do begin
    aia_lct,wave=wvid[i],/load
    plot_map,mapsr[i],/log,charsize=2,$
      dmin=0.5*median(mapsr[i].data),dmax=0.8*max(mapsr[i].data),$
      title=wvid[i]+STRING(197B)
  endfor
  
  map2fits,mapsr,out_dir+'test_maps_rebin_20151225.fits'


  stop

end