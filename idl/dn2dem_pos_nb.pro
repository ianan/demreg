pro dn2dem_pos_nb, dn_in, edn_in,tresp,tresp_logt,temps,$
  dem,edem,elogt,chisq,dn_reg,$
  timed=timed,gloci=gloci,glcindx=glcindx,$
  rgt_fact=rgt_fact, max_iter=max_iter,reg_tweak=reg_tweak,dem_norm0=dem_norm0

  ; Performs a Regularization on solar data, returning the Differential Emission Measure (DEM)
  ; using the method of Hannah & Kontar A&A 553 2013
  ;
  ; Basically getting DEM(T) out of g(f)=K(f,T)#DEM(T)
  ;
  ; This code is a wrapper for your input data to the regularization code demmap_pos.pro
  ;
  ; Previously the wrapper was AIA map specific and parallised via IDL-Bridges, i.e.
  ; http://www.astro.gla.ac.uk/~iain/demreg/map/
  ; But this version takes any data in terms of multiple filters or lines and single pixel or maps/arrays
  ; also no parallisation at the moment
  ;
  ; Inputs:
  ;                   Number of x-pixels          nx
  ;                   Number of y-pixels          ny
  ;                   Number of lines/filters     nf
  ;                   Number of Resp temp bins    nt0
  ;                   Number of DEM temp bins     nt
  ;
  ;
  ;   dn_in      -    Data in appropriate units for response function (i.e. DN/s/px)
  ;                         Either a single pixel of multiple filters/lines (nf) or 2d array/map (nx,ny,nf)
  ;   edn_in     -    Error in dn_in
  ;                         Either a single pixel of multiple filters/lines (nf) or 2d array/map (nx,ny,nf)
  ;   tresp      -    Temperature response function in structure (nt0,nf) in appropiate units (i.e. DN/cm^5/s/px)
  ;   tresp_logt -    Temperature (log) binning of the response function (nt0)
  ;   temps      -    Temperature bin edges of the output DEM (nt+1) - NOT LOG10
  ;
  ;
  ; Outputs:
  ;
  ;   dem      -    Regularized DEM in units of 1/cm^5/K and array of (nx,ny,nt)
  ;   edem     -    Error (vertical) in the DEM
  ;   elogt    -    Error (horiztonal) in the DEM (temperature resolution/deviation from identiy matrix)
  ;   chisq    -    Chisq in data space of DEM (folded through responses) and observations
  ;   dn_reg   -    The data space value of the DEM (folded through responses)
  ;
  ;   Optional Inputs:
  ;
  ;   timed    -    /timed if you want to print out DEM/s
  ;   gloci    -    /gloci if you want to use the min of the EM loci curve for the calculation ???
  ;                       Could if line data or narrow filter responses
  ;   glcindx  -    Can choose just to use gloci on certain filters, array of (nf)
  ;                       i.e. use 4th of 6 filters glcindx=[0,0,0,1,0,0]
  ;   reg_tweak - What chisq trying to get for final solution in DN space (default 1)
  ;   rgt_fact -  Multiplying factor of how much to increase reg tweak each iteration when trying to make positive (default 1.5, shouldn't be more than a few)
  ;   max_iter -  Max number of iterations to try and get a positive solution (default 10)
  ;
  ;   dem_norm0 - Initial guess/constraing/weighting DEM, normalised and array of (nx,ny,nt)
  ;
  ;
  ;    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;   13-Apr-2016 IGH   Updated and tidied version to start further development
  ;                         (Need to optimise L calculation and weighting, provide option to input weighting ??)
  ;   26-Apr-2016 IGH   Added options to change, reg_tweak, rgt_fact and max_iter
  ;   27-Apr-2016 IGH   Added in option to supply initial guess/constraint normalized DEM to weight L
  ;                     Fixed bug where only the first pixel of dem_norm0 was sent to demmap_pos
  ;
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  if (n_elements(temps) lt 6) then temps=[0.5,1,1.5,2,3,4,6,8,11,14,19,25,32]*1e6

  logT=get_edges(alog10(temps),/mean)
  dlogT=get_edges(alog10(temps),/width)
  nt=n_elements(logT)

  ; Hopefully this means you can input a variety of data (single pixel or x,y arrays)
  sze=size(dn_in)
  if (sze[0] eq 1) then begin
    nx=1
    ny=1
    nf=sze[1]
    dn=dblarr(1,1,nf)
    dn[0,0,*]=dn_in
    edn=dblarr(1,1,nf)
    edn[0,0,*]=edn_in
    if ((size(dem_norm0))[0] gt 0) then begin
      dem0=dblarr(1,1,nt)
      dem0[0,0,*]=dem_norm0
    endif
  endif
  if (sze[0] eq 2) then begin
    nx=sze[1]
    ny=1
    nf=sze[2]
    dn=dblarr(nx,1,nf)
    dn[*,0,*]=dn_in
    edn=dblarr(nx,1,nf)
    edn[*,0,*]=edn_in
    if ((size(dem_norm0))[0] gt 0) then begin
      dem0=dblarr(nx,1,nt)
      dem0[*,0,*]=dem_norm0
    endif
  endif
  if (sze[0] eq 3) then begin
    nx=sze[1]
    ny=sze[2]
    nf=sze[3]
    dn=dn_in
    edn=edn_in
    if ((size(dem_norm0))[0] gt 0) then dem0=dem_norm0
  endif

  ;**************
  ; Do the gloci for the initial L calculation ?
  ; Either do it for all the filters/lines or just those specified

  if keyword_set(gloci) then begin
    if (n_elements(glcindx) eq nf) then glc=glcindx else glc=intarr(nf)+1
  endif else begin
    glc=intarr(nf)
  endelse

  if (n_elements(tresp[0,*]) ne nf) then $
    print, 'Tresp needs to be the same number of wavelengths/filters as the data.'

  truse=dblarr(n_elements(tresp[*,0]),nf)
  for i=0,nf-1 do begin
    goodtr=where(tresp[*,i] gt 0.)
    badtr=where(tresp[*,i] le 0.,nbadtr)
    truse[goodtr,i]=tresp[goodtr,i]
    if (nbadtr gt 0) then truse[badtr,i]=min(tresp[goodtr,i])
  endfor
  TR=dblarr(nt,nf)
  for i=0, nf-1 do TR[*,i]=interpol(truse[*,i], tresp_logt, logT)

  RMatrix=dblarr(nt,nf)
  ; Put in the 1/K factor (remember doing it in logT not T hence the extra terms)
  for i=0, nf-1  do RMatrix[*,i]=TR[*,i]*10d^logT*alog(10d^dlogT)
  ; Just scale so not dealing with tiny numbers
  sclf=1d15
  RMatrix=RMatrix*sclf

  ; Just set some other parms up before sending to the actual demmap_pos.pro
  if keyword_set(timed) then tnow=systime(1)
  if (n_elements(reg_tweak) lt 1) then reg_tweak=1.0
  if (n_elements(max_iter) lt 1) then max_iter=10
  if (n_elements(rgt_fact) lt 1) then rgt_fact=1.5

  ; need to convert dn and edn from [x,y,f] to [x*y,f]
  dn1d=dblarr(nx*ny,nf)
  edn1d=dblarr(nx*ny,nf)

  ; Convert 2d maps per filter into 1d array
  ; Was like this for the parallisation/bridges even though current version not using bridges
  for yy=0, ny-1 do begin
    for xx=0, nx-1 do begin
      dn1d[yy*nx+xx,*]=reform(dn[xx,yy,*])
      edn1d[yy*nx+xx,*]=reform(edn[xx,yy,*])
    endfor
  endfor

  dem1d=dblarr(nx*ny,nt)
  chisq1d=dblarr(nx*ny)
  edem1d=dblarr(nx*ny,nt)
  elogt1d=dblarr(nx*ny,nt)
  dn_reg1d=dblarr(nx*ny,nf)

 ;*****************************************************
; Actually doing the DEM calculations
;*****************************************************
  ; Do we have an initial DEM guess/constraint to send to demmap_pos as well?
  if ( (size(dem0))[0] eq (size(dn))[0]) then begin
    dem01d=dblarr(nx*ny,nt)
    for yy=0, ny-1 do begin
      for xx=0, nx-1 do begin
        dem01d[yy*nx+xx,*]=reform(dem0[xx,yy,*])
      endfor
    endfor
    demmap_pos,dn1d,edn1d,RMatrix,logt,dlogt,glc,$
      dem1d,chisq1d,edem1d,elogt1d,dn_reg1d,$
      reg_tweak=reg_tweak,max_iter=max_iter,rgt_fact=rgt_fact,$
      dem_norm0=dem01d
  endif else begin
    demmap_pos,dn1d,edn1d,RMatrix,logt,dlogt,glc,$
      dem1d,chisq1d,edem1d,elogt1d,dn_reg1d,$
      reg_tweak=reg_tweak,max_iter=max_iter,rgt_fact=rgt_fact
  endelse
    ;*****************************************************
    ;*****************************************************

  dem=dblarr(nx,ny,nt)
  chisq=dblarr(nx,ny)
  edem=dblarr(nx,ny,nt)
  elogt=dblarr(nx,ny,nt)
  dn_reg=dblarr(nx,ny,nf)

  ; Change back into 2d from 1d
  for yy=0, ny-1 do begin
    for xx=0, nx-1 do begin
      dem[xx,yy,*]=reform(dem1d[yy*nx+xx,*])*sclf
      edem[xx,yy,*]=reform(edem1d[yy*nx+xx,*])*sclf
      elogt[xx,yy,*]=reform(elogt1d[yy*nx+xx,*])/(2.0*sqrt(2.*alog(2.)))
      chisq[xx,yy]=reform(chisq1d[yy*nx+xx])
      dn_reg[xx,yy,*]=reform(dn_reg1d[yy*nx+xx,*])
    endfor
  endfor

  dem=reform(dem)
  edem=reform(edem)
  elogt=reform(elogt)
  chisq=reform(chisq)
  dn_reg=reform(dn_reg)


  if keyword_set(timed) then tend=systime(1)
  if keyword_set(timed) then print,string((nx*ny*1.)/(tend-tnow),format='(f9.1)'), ' DEM/s'

end



