FUNCTION data2dem_reg, logT ,TRmatrix ,data ,edata ,$
  mint=mint,maxt=maxt,nt=nt,$
  order=order , guess=guess ,reg_tweak=reg_tweak,$
  channels=channels,debug=debug,gloci=gloci, pos=pos

  ; Recovers an underlying Differential Emission Measure DEM(T) given
  ; multi-filter observations and the corresponding filter's temperature responses
  ; or multi-line observations and the corresponding contribution functions
  ;
  ; The inversion is done using regulariation & GSVD method to obtain a
  ; model-independent solution
  ;
  ; For details see:
  ; 	Hannah and Kontar, A&A 2012 539 146
  ;   http://adsabs.harvard.edu/abs/2012A%26A...539A.146H
  ; 	http://www.astro.gla.ac.uk/~iain/demreg/
  ;
  ; Inversion method used is an updated and optimised version of method to invert a RHESSI DN spectrum
  ; given the RHESSI DRM to the source photon spectrum and uses some of the same software in ssw
  ; $SSW/packages/xray/idl/inversion/function_DN2photons.pro
  ;
  ; INPUTS:
  ;
  ; 	logT  		- Temperatures where the filter responses/contribution functions are known [log10(K)]
  ; 	TRmatrix 	- Filter Responses as a function of temperature for each filter (nt by nf) [DN/s/px]
  ; 				- or contribution functions if working with spectroscopic data [erg cm^3 s^-1 sr^-1]
  ; 	data 		- Data Number Rate per filter (nf) [DN/s/px]
  ; 				- or line intensity [erg cm^-2 s^-1 sr^-1]
  ; 	edata 		- Error in Data Number Rate per filter (nf) [DN/s/px]
  ; 			  	- or error in line intensity [erg cm^-2 s^-1 sr^-1]
  ; 	channels	- String array labelling each filter/line i.e ['94',131' etc]
  ;
  ;
  ; 	OPTIONAL INPUTS:
  ; 	minT     		- float, minimum temperature in logT space: Default=5.7
  ; 	maxT     		- float, maximum temperature in logT space: Default=7.3
  ; 	nt     			- integer, number of temperature bins: Default=33
  ; 	order     		- integer, regularization order. Can be 0 (default), 1, or 2.
  ; 	reg_tweak		- float,  parameter to adjust the regularization parameter, chisq in DEM space (default is 1).
  ; 	guess     		- integer, 1 or 0 (default) - to use or NOT to use guess solution in second run
  ;   gloci 			- integer, 1 or 0 (default) - Use the min of EM Loci curves for the guess DEM_model
  ; 	pos 			- integer, 1 or 0 (default) - Produce regularized and regularized positive solution?
  ;  					- note that if the default approach gives a postive soultion then get dem and dem_pos outputted
  ;
  ; OUTPUTS:

  ; The results of the inversion are returned as a structure
  ; main ones of interest:
  ;
  ; 	DEM      		- regularised DEM, [cm^-5 K^-1]
  ; 	eDEM         	- 1 sigma (vertical) error, [cm^-5 K^-1]
  ; 	eLOGT       	- energy resolution of the regularised solution (horizontal error) [log T(K)]
  ;
  ;	If the pos option has been used the additional positive solution is denoted by *_pos in the result structure
  ;
  ; OPTIONAL OUTPUTS:
  ;   none
  ;
  ; KEYWORDS:
  ;   none
  ;
  ; COMMON BLOCKS:
  ;   none
  ;
  ; SIDE EFFECTS:
  ;
  ;
  ; RESTRICTIONS:
  ;
  ;
  ; MODIFICATION HISTORY:
  ; 06-Sep-2010 Rewrote and tidied up orginal IGH
  ; 21-Jul-2011 Major rewrite and tidy IGH
  ; 29-Jul-2011 interpol on log(TRmatrix) not (TRmatrix) as smoother IGH
  ; 03-Aug-2011 added gloci as alternative guess method IGH
  ; 17-Aug-2011 added additional temp res calculation (fwhm2, elogt2) IGH
  ; 18-Nov-2011 merged broadband and line codes into single routine (data2dem_reg.pro) IGH
  ; 18-Nov-2011 added positivity constraint option IGH
  ; 25-Jan-2012 code no longer returns some unneccesary diagnostics IGH
  ; 12-Nov-2012 output data_cont_t in res structure (contribution to data per each T) IGH
  ; 15-May-2013 correct typo in header file IGH
  ; 02-Nov-2013 added the conversion factors DEM to EM (EM=DEM*DEMtoEM) as output in reg structure
  ; 13-Nov-2013 added warning if any of errors=0 (GSVD will fail as well)
  ; 13-Nov-2013 added test for data=0 when using emloci (replaces 0 with tiny fraction of error)
  ; 13-Sep-2013 removed nt in demtoem factor and emloci scaling
  ; 11-Oct-2016 removed get_edges() call and calculated logT bin-mids directly
  ; 20-Dec-2017 minor bug fix for missing data_cont_t_pos definition
  ; 28-Apr-2020  Changed any fltarr() to dblarr()
  ;
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (n_elements(order) lt 1) then order=0
  if (n_elements(guess) lt 1) then guess=0
  if (n_elements(reg_tweak) lt 1) then reg_tweak=1.
  if (n_elements(nt) lt 1) then nt=33
  if (n_elements(mint) lt 1) then mint=5.7
  if (n_elements(maxt) lt 1) then maxt=7.3
  nf=n_elements(data)
  if (n_elements(channels) lt nf) then channels=strarr(nf)
  if (n_elements(pos) lt 1) then pos=0

  if (nf ge nt) then print, '!!!!! Warning: n_T must be bigger than n_F'
  if (nf ne n_elements(edata)) then print, '!!!!! Warning: Data and eData must be the same size'
  if ((where(edata le 0.))[0] ne -1) then print, '!!!!! Warning: None of eData can be 0'

  ;interpolate temperature responses to binning desired of DEM
  logTint=mint+(maxt-mint)*findgen(nt)/(nt-1.0)
  TRmatint=dblarr(nt,nf)
  ;better interpolation if done in log space though need check not -NaN
  for i=0, nf-1 do TRmatint[*,i]=10d^interpol(alog10(TRmatrix[*,i]), logT, logTint) > 0.
  ; also get response in correct units and scale to make numerics simplier
  dlogTint  =logTint[1:nT-1]-logTint[0:nT-2]
  dlogTint=[dlogTint,dlogTint[nt-2]]

  ;Before we begin can work out the conversion factors for DEM to EM
  lgt_edg=dblarr(nt-1)
  for ll=0,nt-2 do lgt_edg[ll]=(logTint[ll+1]+logTint[ll])*0.5
  lgt_edg=[lgt_edg[0]-(lgt_edg[1]-lgt_edg[0]),lgt_edg,lgt_edg[nt-2]+(lgt_edg[nt-2]-lgt_edg[nt-3])]
  dlgT_edg=lgt_edg[1:nt]-lgt_edg[0:nt-1]
  DEMtoEM=10d^lgt_edg*alog(10d^dlgt_edg)

  ; Now intrepolate the response functions to the temperature binning of DEM output
  RMatrix=dblarr(nT,nF)
  for i=0, nF-1  do RMatrix[*,i]=TRmatint[*,i]*10d^logTint*alog(10d^dlogTint)
  RMatrix=RMatrix*1d20
  RMatrix_org=Rmatrix
  DEM_model =dblarr(nt)

  ; normalize everything by the errors
  data_in=data/edata
  edata_in=edata/edata
  for i=0, nf-1 do RMatrix[*,i]=RMatrix[*,i]/edata[i]

  if keyword_set(gloci) then begin
    ; can choose to only run regularization once and that the guess solution to be min of EM loci
    ; this is then used for weighting the constraint matrix
    emloci=dblarr(nt,nf)
    ; Test to see if any of the data is 0, if so use a tiny fraction of error for emloci
    dd4eml=data
    test_0_data=where(data le 0.,nt0d)
    if (nt0d gt 0) then dd4eml[test_0_data]=1e-6*edata[test_0_data]
    for ii=0, nf-1 do emloci[*,ii]=dd4eml[ii]/(TRmatint[*,ii]*10d^logtint*alog(10d^dlogTint))
    for ii=0, nf-1 do emloci[*,ii]=data[ii]/(TRmatint[*,ii]*10d^logtint*alog(10d^dlogTint))
    for jj=0, nt-1 do dem_model[jj]=min(emloci[jj,*])
    dem_model=smooth(dem_model,3)*1d-20
  endif else begin
    ;********************************************************************************************************
    ; Regularization is run twice
    ; the first time the constraint matrix is taken as the identity matrix normalized by dlogT
    ; must use this 0th order as no dem_model guess
    L=dblarr(nT,nT)
    for i=0, nT-1 do L[i,i]=1.0/sqrt(dlogTint[i])

    ; GSVD on temperature responses (Rmatrix) and constraint matrix (L)
    dem_inv_gsvdcsq,RMatrix,L,Alpha,Betta,U,V,W

    ; Determine the regularization parameter
    ; for first run using weekly regularized with reg_tweak=sqt(nt)
    dem_inv_reg_parameter,Alpha,Betta,U,W,data_in,edata_in,transpose(DEM_model)*Guess,sqrt(nt*1.0),opt

    ; Now work out the regularized DEM solution
    dem_inv_reg_solution,Alpha,Betta,U,W,data_in,opt,DEM_model*Guess,DEM_reg
    ;********************************************************************************************************
    ; For second run use found regularized solution as weighting for constraint matrix and possible guess solution
    DEM_reg=DEM_reg *(DEM_reg GT 0)+1e-4*max(Dem_reg)*(DEM_reg LT 0)
    DEM_reg=smooth(DEM_reg,3)
    DEM_model=DEM_reg
  endelse

  ;This time make constraint to specified order
  dem_inv_make_constraint,L,logTint,dlogTint,DEM_model,order

  ; GSVD on temperature responses (Rmatrix) and constraint matrix (L)
  dem_inv_gsvdcsq,RMatrix,L,Alpha,Betta,U,V,W

  if keyword_set(pos) then begin

    ; To find the positive solution all solutions are found for a range of regularization parameters
    ; Then one with chi^2 closest to reg_tweak and gives positive DEM is given as _pos solution

    ; Determine the regularization parameter
    ; for second run regularize to level specified with reg_tweak
    dem_inv_reg_parameter_pos,Alpha,Betta,U,W,data_in,edata_in,transpose(DEM_model)*Guess,$
      reg_tweak,opt,dem_reg,opt_pos,dem_reg_pos

    ;********************************************************************************************************
    ; now work out the temperature resolution/horizontal spread
    dem_inv_reg_resolution,Alpha,Betta,opt,W,logTint,dlogTint,FWHM,cent,RK,fwhm2
    dem_inv_reg_resolution,Alpha,Betta,opt_pos,W,logTint,dlogTint,FWHM_pos,cent_pos,RK_pos,fwhm2_pos

    ; now work out the DEM error (vertical error)
    npass=300.
    dem_inv_confidence_interval,DEM_reg,data_in,edata_in,Alpha,Betta,U,W,opt,DEM_model,Guess,Npass,reg_sol_err
    dem_inv_confidence_interval,DEM_reg_pos,data_in,edata_in,Alpha,Betta,U,W,opt_pos,DEM_model,Guess,Npass,reg_sol_err_pos
    ;********************************************************************************************************
    ;Calculate the data signal that found regularized DEM gives you and work out data residuals and chisq
    data_reg=transpose(Rmatrix_org##dem_reg)
    residuals=(data-data_reg)/edata
    chisq=total(residuals^2)/(nf*1.0)

    data_reg_pos=transpose(Rmatrix_org##dem_reg_pos)
    residuals_pos=(data-data_reg_pos)/edata
    chisq_pos=total(residuals_pos^2)/(nf*1.0)

    data_cont_t=dblarr(nt,nf)
    for i=0, nf-1 do data_cont_t[*,i]=Rmatrix_org[*,i]*dem_reg
    data_cont_t_pos=dblarr(nt,nf)
    for i=0, nf-1 do data_cont_t_pos[*,i]=Rmatrix_org[*,i]*dem_reg_pos

    reg_solution={data:data,edata:edata, Tresp:Trmatint,channels:channels,$
      DEM:DEM_reg[0:nT-1]*1d20, eDEM:reg_sol_err[0:nT-1]*1d20, $
      logT:logTint,elogT:fwhm/(2.0*sqrt(2*alog(2))),RK:RK,$
      data_reg:data_reg,residuals:residuals,chisq:chisq,data_cont_t:data_cont_t,$
      DEM_pos:DEM_reg_pos[0:nT-1]*1d20, eDEM_pos:reg_sol_err_pos[0:nT-1]*1d20, $
      elogT_pos:fwhm_pos/(2.0*sqrt(2*alog(2))),RK_pos:RK_pos,$
      data_reg_pos:data_reg_pos,residuals_pos:residuals_pos,chisq_pos:chisq_pos,$
      data_cont_t_pos:data_cont_t_pos,$
      reg_tweak:reg_tweak,guess:guess,order:order,DEMtoEM:DEMtoEM}

  endif else begin
    ; Here do not require positive solution so first find regularization parameter than gives chosen chi^2
    ; Then use it to work out the regularized solution

    ; Determine the regularization parameter
    ; for second run regularize to level specified with reg_tweak
    dem_inv_reg_parameter,Alpha,Betta,U,W,data_in,edata_in,transpose(DEM_model)*Guess,reg_tweak,opt

    ; Now work out the regularized DEM solution
    dem_inv_reg_solution,Alpha,Betta,U,W,data_in,opt,DEM_model*Guess,DEM_reg

    ;********************************************************************************************************
    ; now work out the temperature resolution/horizontal spread
    dem_inv_reg_resolution,Alpha,Betta,opt,W,logTint,dlogTint,FWHM,cent,RK,fwhm2

    ; now work out the DEM error (vertical error)
    npass=300.
    dem_inv_confidence_interval,DEM_reg,data_in,edata_in,Alpha,Betta,U,W,opt,DEM_model,Guess,Npass,reg_sol_err
    ;********************************************************************************************************
    ;Calculate the data signal that found regularized DEM gives you and work out data residuals and chisq
    data_reg=transpose(Rmatrix_org##dem_reg)
    residuals=(data-data_reg)/edata
    chisq=total(residuals^2)/(nf*1.0)

    data_cont_t=dblarr(nt,nf)
    for i=0, nf-1 do data_cont_t[*,i]=Rmatrix_org[*,i]*dem_reg

    reg_solution={data:data,edata:edata, Tresp:Trmatint,channels:channels,$
      DEM:DEM_reg[0:nT-1]*1d20, eDEM:reg_sol_err[0:nT-1]*1d20, $
      logT:logTint,elogT:fwhm/(2.0*sqrt(2*alog(2))),RK:RK,$
      data_reg:data_reg,residuals:residuals,chisq:chisq,data_cont_t:data_cont_t,$
      reg_tweak:reg_tweak,guess:guess,order:order,DEMtoEM:DEMtoEM}

  endelse

  return,reg_solution

  if keyword_set(debug) then stop

end

