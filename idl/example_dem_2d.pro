pro example_dem_2d

  ; Example script to recover the DEM from AIA array of data
  ; The AIA data is synthetic for specified Gaussian DEM model
  ;
  ; 13-Apr-2016 IGH
  ; 27-Apr-2016 IGH   - Changed the naming of the temperatures to make things clearer:
  ;                     tr_logt is the binning of the response function
  ;                     temps is the bin edges you want for the DEM
  ;                     logtemps is the log of the above
  ;                     mlogt is the mid_point of the above bins
  ; 20-May-2019 IGH   - Minor update: must use the timedepend_date option when getting aia response   
  ; 15-Aug-2019 IGH   - Changed file name *_demmap_* to *_map_*
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  d10=4d22
  m10=6.2
  s10=0.05

  ; What temperature binning do you want of the DEM?
  ; These are the bin edges
  ;  temps=[0.5,1,1.5,2,3,4,6,8,11,14,19,25,32]*1e6
  temps=[0.5,1,1.5,2,3,4,6,8,11,14,19]*1d6
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

  nx=75
  ny=75
  dn=dblarr(nx,ny,6)
  edn=dblarr(nx,ny,6)
  demmods=dblarr(nx,ny,n_elements(tr_logt))

  m1s=fltarr(nx)
  s1s=fltarr(ny)

  ; Make the Gausian model be different in each pixel

  ; For the error calc
  gains=[18.3,17.6,17.7,18.3,18.3,17.6]
  dn2ph=gains*[94,131,171,193,211,335]/3397.
  rdnse=[1.14,1.18,1.15,1.20,1.20,1.18]

  for ii=0,nx-1 do begin
    for jj=0, ny-1 do begin
      m1=m10+ii*(1/(nx-1.0))
      s1=s10+jj*(0.15/(ny-1.0))
      m1s[ii]=m1
      s1s[jj]=s1
      d1=d10
      dem_mod=(d1/(root2pi*s1))*exp(-(tr_logt-m1)^2/(2*s1^2))
      demmods[ii,jj,*]=dem_mod
      dn_mod=dem2dn(tr_logt, dem_mod, TRmatrix)
      dn[ii,jj,*]=dn_mod.dn
      ; assuming 2.9s duration
      shotnoise=sqrt(dn2ph*dn_mod.dn*2.9)/dn2ph/2.9
      ; add in a small systematic as well ?
      ; small errors closer to model but ends up with more negative solutions ;(
      edn[ii,jj,*]=sqrt(rdnse^2+shotnoise^2)+0.01*dn_mod.dn
    endfor
  endfor

  dn2dem_pos_nb, dn, edn,TRmatrix,tr_logt,temps,dem,edem,elogt,chisq,dn_reg,/timed;,rgt_fact=2

  ; Plot one of the temperature bins
  loadct,39
  t=5
  window,0,xsize=650,ysize=500
  plot_image,dem[*,*,t],xtitle='log T',ytitle='sig',origin=[min(m1s),min(s1s)],scale=[m1s[1]-m1s[0],s1s[1]-s1s[0]],$
    title=string(logtemps[t],format='(f4.2)')+' to '+string(logtemps[t+1],format='(f4.2)')+' Log!D10!N MK',/nosquare


  window,1,xsize=650,ysize=500,title='Regularized DEM'
  x=5
  y=7
  plot,mlogt,dem[x,y,*],/ylog,yrange=max(dem[x,y,*])*[1e-2,5], tit='logT: '+string(m1s[x],format='(f4.1)')+', sig: '+string(s1s[y],format='(f5.3)')
  for i=0, n_elements(mlogt)-1 do oplot,alog10(temps[i:i+1]),dem[x,y,i]*[1,1],color=200
  for i=0, n_elements(mlogt)-1 do oplot,mlogt[i]*[1,1],dem[x,y,i]+edem[x,y,i]*[-1.,1.],color=200
  oplot,tr_logt,demmods[x,y,*],color=120

;  print,reform(dem[x,y,*])


  stop
end
