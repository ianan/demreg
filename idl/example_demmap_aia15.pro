pro example_demmap_aia15

  ; Example script to recover the DEM from AIA Lvl1.5 fits files
  ; The specific AIA fits used here are not include with the code
  ;
  ; 14-Apr-2015 IGH
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
  temps=[0.5,1,1.5,2,3,4,6,8,11,14,19,25,32]*1e6
  logt0=alog10(temps)

  ; Need to make the response functions?
  if (file_test('aia_resp.dat') eq 0) then begin
    tresp=aia_get_response(/temperature,/dn,/chianti,/noblend,/evenorm)
    save,file='aia_resp.dat',tresp
  endif else begin
    restore,file='aia_resp.dat'
  endelse

  ; Only want the coronal ones without 304A
  idc=[0,1,2,3,4,6]

  logT=tresp.logte
  gdt=where(logt ge min(logt0) and logt le max(logt0),ngd)
  logt=logt[gdt]
  TRmatrix=tresp.all[*,idc]
  TRmatrix=TRmatrix[gdt,*]
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ; Just do a sub-part of the image for testing purposes
  dn=dn0[261:360,311:410,*]
  edn=edn0[261:360,311:410,*]
  dn2dem_pos_nb, dn, edn,TRmatrix,logt,temps,dem,edem,elogt,chisq,dn_reg,/timed

  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ;  ; Plot one of the temperature bins
  ;  nt=n_elements(dem[0,0,*])
  ;  mdem0=make_map(dblarr(n_elements(dem[*,0,0])),dblarr(n_elements(dem[0,*,0])))
  ;
  ;  mdem=replicate(mm[0],nt)
  ;  for i=0,nt-1 do begin
  ;    mdem[i].data=dem[*,*,0]
  ;    mdem[i].id=string(logt0[i],format='(f4.2)')+' to '+string(logt0[i+1],format='(f4.2)')+' Log!D10!N MK'
  ;  endfor


  loadct,39
  t=8
  plot_image,dem[*,*,t],$
    title=string(logt0[t],format='(f4.2)')+' to '+string(logt0[t+1],format='(f4.2)')+' Log!D10!N MK'

  stop
end