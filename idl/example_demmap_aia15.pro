pro example_demmap_aia15

  ; Example script to recover the DEM from AIA Lvl1.5 fits files
  ; The specific AIA fits used here are not include with the code
  ; 
  ; 14-Apr-2015 IGH
  ;
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  fdir='~/Downloads/aia15_test_full/'
  f094=file_search(fdir, '*94_.fts')
  f131=file_search(fdir, '*131_.fts')
  f171=file_search(fdir, '*171_.fts')
  f193=file_search(fdir, '*193_.fts')
  f211=file_search(fdir, '*211_.fts')
  f335=file_search(fdir, '*335_.fts')
  
  fits2map,f094[0],m094
  ; Get rid of negative values before we begin
  idn=where(m094.data lt 0,nid)
  if (nid gt 1) then m094.data[idn]=0
  
  ; make the map smaller - easier to handle for testing
  ; output is still DN/s/px? Now with bigger pixels
  rm094=rebin_map(m094,1024,1024)

stop

  d10=1d22
  m10=6.5
  s10=0.05
  
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
  root2pi=sqrt(2.*!PI)

  nx=30
  ny=30
  dn=dblarr(nx,ny,6)
  edn=dblarr(nx,ny,6)

  m1s=fltarr(nx)
  s1s=fltarr(ny)

  ; Make the Gausian model be different in each pixel
  for ii=0,nx-1 do begin
    for jj=0, ny-1 do begin
      m1=m10+ii*(1/(nx-1.0))
      s1=s10+jj*(0.15/(ny-1.0))
      m1s[ii]=m1
      s1s[jj]=s1
      d1=d10
      dem_mod=(d1/(root2pi*s1))*exp(-(logT-m1)^2/(2*s1^2))
      dn_mod=dem2dn(logT, dem_mod, TRmatrix)
      dn[ii,jj,*]=dn_mod.dn
      edn[ii,jj,*]=0.1*dn_mod.dn
    endfor
  endfor

  dn2dem_pos_nb, dn, edn,TRmatrix,logt,temps,dem,edem,elogt,chisq,dn_reg,/timed
 

  ; Plot one of the temperature bins
  loadct,39
  t=8
  plot_image,dem[*,*,t],xtitle='m1',ytitle='s1',origin=[min(m1s),min(s1s)],scale=[m1s[1]-m1s[0],s1s[1]-s1s[0]],$
   title=string(logt0[t],format='(f4.2)')+' to '+string(logt0[t+1],format='(f4.2)')+' Log!D10!N MK'

  stop
end