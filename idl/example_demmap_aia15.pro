pro example_demmap_aia15

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
  ;
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ; Initial get the data, remove neagtives and rebin to smaller resolution for testing
  ; Also note the int to float for the AIA data is done during the rebinning
  ; If not rebinning still need to do this.
  fdir='~/Downloads/aia15_test_full/'
  waves=['094','131','171','193','211','335']
  nf=n_elements(waves)

  ;  f094=file_search(fdir, '*94_.fts')
  ;  f131=file_search(fdir, '*131_.fts')
  ;  f171=file_search(fdir, '*171_.fts')
  ;  f193=file_search(fdir, '*193_.fts')
  ;  f211=file_search(fdir, '*211_.fts')
  ;  f335=file_search(fdir, '*335_.fts')
  ;  ff=[f094[0],f131[0],f171[0],f193[0],f211[0],f335[0]]
  ;
  ;  for i=0, nf-1 do begin
  ;     fits2map,ff[i],map
  ;     ; Get rid of negative values before we begin
  ;     idn=where(map.data lt 0,nid)
  ;     if (nid gt 1) then map.data[idn]=0
  ;
  ;     ; make the map smaller - easier to handle for testing
  ;     ; output is still DN/px? Now with bigger pixels
  ;     rmap=rebin_map(map,1024,1024)
  ;
  ;     ; save it out
  ;     map2fits,rmap,fdir+'test_aia15_1024_'+waves[i]+'A.fts'
  ;  endfor

  ; Load back in the rebinned data
  f2=file_search(fdir,'*A.fts')
  fits2map, f2,mm

  ; Setup the data for input to DEMREG code
  dn0=mm.data
  durs=mm.dur
  ; Get into DN/s/px
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
    shotnoise=sqrt(dn2ph[i]*abs(dn0[*,*,i])*durs[i])/dn2ph[i]
    edn0[*,*,i]=sqrt(rdnse[i]^2.+shotnoise^2.)/durs[i]
  endfor

  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; What temperature binning do you want of the DEM?
  ; These are the bin edges
  ;  temps=[0.5,1,1.5,2,3,4,6,8,11,14,19,25,32]*1e6
  ;  temps=[0.5,1,1.5,2,3,4,6,8,11,14,19]*1d6
  temps=[0.5,1,2,4,6,9,14]*1d6
  logtemps=alog10(temps)
  ; This is is the temperature bin mid-points
  mlogt=get_edges(logtemps,/mean)
  nt=n_elements(mlogt)

  ; Need to make the response functions?
  if (file_test('aia_resp.dat') eq 0) then begin
    tresp=aia_get_response(/temperature,/dn,/chianti,/noblend,/evenorm)
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
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ; Just do a sub-part of the image for testing purposes
  dn=dn0[261:360,311:410,*]
  edn=edn0[261:360,311:410,*]
  
  nxn=n_elements(dn[*,0,0])
  nyn=n_elements(dn[0,*,0])
  dem_norm0=dblarr(nxn,nyn,nt)

  ; Bad initial normalisation but just testing if the code works
  for xx=0,nxn-1 do begin
    for yy=0,nyn-1 do begin
      dem_norm0[xx,yy,*]=[1e-2,4.2e3,5e3,1e3,1e2,1e-2]
    endfor
  endfor

  dn2dem_pos_nb, dn, edn,TRmatrix,tr_logt,temps,dem,edem,elogt,chisq,dn_reg,/timed;,dem_norm0=dem_norm0

  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ;  ; Plot one of the temperature bins
  ;  nt=n_elements(dem[0,0,*])
  ;  mdem0=make_map(dblarr(n_elements(dem[*,0,0])),dblarr(n_elements(dem[0,*,0])))
  ;
  ;  mdem=replicate(mm[0],nt)
  ;  for i=0,nt-1 do begin
  ;    mdem[i].data=dem[*,*,0]
  ;    mdem[i].id=string(logtemps[i],format='(f4.2)')+' to '+string(logtemps[i+1],format='(f4.2)')+' Log!D10!N MK'
  ;  endfor


  loadct,39

  !p.multi=[0,3,nt/3]
  ; Plot them all with the same scaling
  for t=0,nt-1 do plot_image,alog10(dem[*,*,t]),chars=2,max=23,min=19,$
    title=string(temps[t]*1d-6,format='(f4.1)')+' to '+string(temps[t+1]*1d-6,format='(f4.1)')+' MK'

  stop
end
