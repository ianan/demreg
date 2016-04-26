pro example_demmap_2d

  ; Example script to recover the DEM from AIA array of data
  ; The AIA data is synthetic for specified Gaussian DEM model
  ;
  ; 13-Apr-2015 IGH
  ;
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  d10=4d22
  m10=6.2
  s10=0.05

  ; What temperature binning do you want of the DEM?
  temps=[0.5,1,1.5,2,3,4,6,8,11,14,19,25,32]*1e6
  temps=[0.5,1,1.5,2,3,4,6,8,11,14,19]*1d6
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

  nx=50
  ny=50
  dn=dblarr(nx,ny,6)
  edn=dblarr(nx,ny,6)
  demmods=dblarr(nx,ny,n_elements(logt))

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
      dem_mod=(d1/(root2pi*s1))*exp(-(logT-m1)^2/(2*s1^2))
      demmods[ii,jj,*]=dem_mod
      dn_mod=dem2dn(logT, dem_mod, TRmatrix)
      dn[ii,jj,*]=dn_mod.dn
      ; assuming 2.9s duration
      shotnoise=sqrt(dn2ph*dn_mod.dn*2.9)/dn2ph/2.9
      ; add in a small systematic as well ?
      ; small errors closer to model but ends up with more negative solutions ;(
      edn[ii,jj,*]=sqrt(rdnse^2+shotnoise^2)+0.01*dn_mod.dn
    endfor
  endfor

  dn2dem_pos_nb, dn, edn,TRmatrix,logt,temps,dem,edem,elogt,chisq,dn_reg,/timed;,rgt_fact=2

  ; Plot one of the temperature bins
  loadct,39
  t=5
  window,0,xsize=650,ysize=500
  plot_image,dem[*,*,t],xtitle='log T',ytitle='sig',origin=[min(m1s),min(s1s)],scale=[m1s[1]-m1s[0],s1s[1]-s1s[0]],$
    title=string(logt0[t],format='(f4.2)')+' to '+string(logt0[t+1],format='(f4.2)')+' Log!D10!N MK',/nosquare


  window,1,xsize=650,ysize=500,title='Regularized DEM'
  lgt_out=get_edges(alog10(temps),/mean)
  x=5
  y=7
  plot,lgt_out,dem[x,y,*],/ylog,yrange=max(dem[x,y,*])*[1e-2,5], tit='logT: '+string(m1s[x],format='(f4.1)')+', sig: '+string(s1s[y],format='(f5.3)')
  for i=0, n_elements(lgt_out)-1 do oplot,alog10(temps[i:i+1]),dem[x,y,i]*[1,1],color=200
  for i=0, n_elements(lgt_out)-1 do oplot,lgt_out[i]*[1,1],dem[x,y,i]+edem[x,y,i]*[-1.,1.],color=200
  oplot,logt,demmods[x,y,*],color=120

  print,reform(dem[x,y,*])


  stop
end