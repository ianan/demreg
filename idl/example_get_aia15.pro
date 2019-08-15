pro example_get_aia15

  ; Example script to get and prep some AIA Lvl1.5 fits files
  ; Output of this code can be used with example_dem_aia15
  ; 
  ; 15-Aug-2019 IGH   - Started
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; 
  ; A short time range to get the data from
  tt='2015-12-25 '+['12:00:00','12:00:11']

  ; Do the vso search for sdo/aia over the time
  ; The wave search filter can only do min to max, so will recover 304 which don't want (optically thick)
  ss = vso_search(tt[0],tt[1], source='sdo',instr='aia',wave='94-335')
  ss = ss[where(ss.wave.min ne 304.)]
  
  ; Just double check we have one of each filter type
  ss=ss[uniq(ss.wave.min)]
  
  ; Get the data and put it somewhere
  out_dir='/Users/iain/Desktop/temp/'
  status=vso_get(ss,out_dir=out_dir,/rice,site='NSO')
  
  ; Load in the fits and then prep them
  files=file_search(out_dir, '*.fits')
  read_sdo,files,ind,data
  
  aia_prep,ind,data,indp,datap
  sid=sort(indp.wavelnth)
  
  index2map,indp,datap,maps
  maps=maps[sid]
  
  !p.multi=[0,3,2]
  for i=0, 5 do plot_map,maps[0]
  
  
 
 stop
 
 end