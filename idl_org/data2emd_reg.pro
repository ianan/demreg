FUNCTION data2emd_reg, logT ,TRmatrix ,data ,edata ,$
  mint=mint,maxt=maxt,nt=nt,$
  order=order , guess=guess ,reg_tweak=reg_tweak,$
  channels=channels,debug=debug,gloci=gloci, pos=pos

  ; Recovers an underlying Emission Measure Distribution EMD(T) given
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
  ; 	logT  		- Temperatures where the filter responses/contribution functions are known (nf) [log10(K)]
  ; 	TRmatrix 	- Temperature response as a dunction of each channel/filter, nt by nf
  ; 	            units need to be such that data/Tresp has units of cm^5
  ; 	            i.e. AIA is DN cm^5 s^-1 px^-1, so data input needs to be DN s^-1 px^-1
  ; 	data 		- Data value per channel (nf): any units as long as data/Tresp has units of cm^5
  ; 	edata 		- Error in data value (nf), same units as data
  ; 	channels	- String array labelling each filter/line (nf) i.e ['94',131' etc]
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
  ; 	EMD      		- regularised EMD, [cm^-5]
  ; 	eEMD         	- 1 sigma (vertical) error, [cm^-5]
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
  ; 18-Feb-2020   EMD version of data2dem_reg.pro - so returns cm^-5 instead of cm^-5 K^-1
  ; 28-Apr-2020   Changed any fltarr() to dblarr() 
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

  ;interpolate temperature responses to binning desired of EMD
  logTint=mint+(maxt-mint)*findgen(nt)/(nt-1.0)
  dlogTint  =logTint[1:nT-1]-logTint[0:nT-2]
  dlogTint=[dlogTint,dlogTint[nt-2]]
  TRmatint=dblarr(nt,nf)
  ;better interpolation if done in log space though need check not -NaN
  for i=0, nf-1 do TRmatint[*,i]=10d^interpol(alog10(TRmatrix[*,i]), logT, logTint) > 0.
  ; just apply scaling so not dealing which such small numbers, will be removed at the end
  sclfc=1d20
  RMatrix=TRmatint*sclfc
  RMatrix_org=Rmatrix
  EMD_model =dblarr(nt)

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
    for ii=0, nf-1 do emloci[*,ii]=dd4eml[ii]/RMatrix[*,ii]
    for ii=0, nf-1 do emloci[*,ii]=data[ii]/RMatrix[*,ii]
    for jj=0, nt-1 do EMD_model[jj]=min(emloci[jj,*])
    EMD_model=smooth(EMD_model,3)
  endif else begin
    ;********************************************************************************************************
    ; Regularization is run twice
    ; the first time the constraint matrix is taken as the identity matrix normalized by dlogT
    ; must use this 0th order as no EMD_model guess
    L=dblarr(nT,nT)
    for i=0, nT-1 do L[i,i]=1.0/sqrt(dlogTint[i])

    ; GSVD on temperature responses (Rmatrix) and constraint matrix (L)
    dem_inv_gsvdcsq,RMatrix,L,Alpha,Betta,U,V,W

    ; Determine the regularization parameter
    ; for first run using weekly regularized with reg_tweak=sqt(nt)
    dem_inv_reg_parameter,Alpha,Betta,U,W,data_in,edata_in,transpose(EMD_model)*Guess,sqrt(nt*1.0),opt

    ; Now work out the regularized EMD solution
    dem_inv_reg_solution,Alpha,Betta,U,W,data_in,opt,EMD_model*Guess,EMD_reg
    ;********************************************************************************************************
    ; For second run use found regularized solution as weighting for constraint matrix and possible guess solution
    EMD_reg=EMD_reg *(EMD_reg GT 0)+1e-4*max(EMD_reg)*(EMD_reg LT 0)
    EMD_reg=smooth(EMD_reg,3)
    EMD_model=EMD_reg
  endelse

  ;This time make constraint to specified order
  dem_inv_make_constraint,L,logTint,dlogTint,EMD_model,order

  ; GSVD on temperature responses (Rmatrix) and constraint matrix (L)
  dem_inv_gsvdcsq,RMatrix,L,Alpha,Betta,U,V,W

  if keyword_set(pos) then begin

    ; To find the positive solution all solutions are found for a range of regularization parameters
    ; Then one with chi^2 closest to reg_tweak and gives positive EMD is given as _pos solution

    ; Determine the regularization parameter
    ; for second run regularize to level specified with reg_tweak
    dem_inv_reg_parameter_pos,Alpha,Betta,U,W,data_in,edata_in,transpose(EMD_model)*Guess,$
      reg_tweak,opt,EMD_reg,opt_pos,EMD_reg_pos

    ;********************************************************************************************************
    ; now work out the temperature resolution/horizontal spread
    dem_inv_reg_resolution,Alpha,Betta,opt,W,logTint,dlogTint,FWHM,cent,RK,fwhm2
    dem_inv_reg_resolution,Alpha,Betta,opt_pos,W,logTint,dlogTint,FWHM_pos,cent_pos,RK_pos,fwhm2_pos

    ; now work out the EMD error (vertical error)
    npass=300.
    dem_inv_confidence_interval,EMD_reg,data_in,edata_in,Alpha,Betta,U,W,opt,EMD_model,Guess,Npass,reg_sol_err
    dem_inv_confidence_interval,EMD_reg_pos,data_in,edata_in,Alpha,Betta,U,W,opt_pos,EMD_model,Guess,Npass,reg_sol_err_pos
    ;********************************************************************************************************
    ;Calculate the data signal that found regularized EMD gives you and work out data residuals and chisq
    data_reg=transpose(Rmatrix_org##EMD_reg)
    residuals=(data-data_reg)/edata
    chisq=total(residuals^2)/(nf*1.0)

    data_reg_pos=transpose(Rmatrix_org##EMD_reg_pos)
    residuals_pos=(data-data_reg_pos)/edata
    chisq_pos=total(residuals_pos^2)/(nf*1.0)

    data_cont_t=dblarr(nt,nf)
    for i=0, nf-1 do data_cont_t[*,i]=Rmatrix_org[*,i]*EMD_reg
    data_cont_t_pos=dblarr(nt,nf)
    for i=0, nf-1 do data_cont_t_pos[*,i]=Rmatrix_org[*,i]*EMD_reg_pos


    reg_solution={data:data,edata:edata, Tresp:Trmatint,channels:channels,$
      EMD:EMD_reg[0:nT-1]*sclfc, eEMD:reg_sol_err[0:nT-1]*sclfc, $
      logT:logTint,elogT:fwhm/(2.0*sqrt(2*alog(2))),RK:RK,$
      data_reg:data_reg,residuals:residuals,chisq:chisq,data_cont_t:data_cont_t,$
      EMD_pos:EMD_reg_pos[0:nT-1]*sclfc, eEMD_pos:reg_sol_err_pos[0:nT-1]*sclfc, $
      elogT_pos:fwhm_pos/(2.0*sqrt(2*alog(2))),RK_pos:RK_pos,$
      data_reg_pos:data_reg_pos,residuals_pos:residuals_pos,chisq_pos:chisq_pos,$
      data_cont_t_pos:data_cont_t_pos,$
      reg_tweak:reg_tweak,guess:guess,order:order}

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

    ; now work out the EMD error (vertical error)
    npass=300.
    dem_inv_confidence_interval,EMD_reg,data_in,edata_in,Alpha,Betta,U,W,opt,EMD_model,Guess,Npass,reg_sol_err
    ;********************************************************************************************************
    ;Calculate the data signal that found regularized EMD gives you and work out data residuals and chisq
    data_reg=transpose(Rmatrix_org##EMD_reg)
    residuals=(data-data_reg)/edata
    chisq=total(residuals^2)/(nf*1.0)

    data_cont_t=dblarr(nt,nf)
    for i=0, nf-1 do data_cont_t[*,i]=Rmatrix_org[*,i]*EMD_reg

    reg_solution={data:data,edata:edata, Tresp:Trmatint,channels:channels,$
      EMD:EMD_reg[0:nT-1]*sclfc, eEMD:reg_sol_err[0:nT-1]*sclfc, $
      logT:logTint,elogT:fwhm/(2.0*sqrt(2*alog(2))),RK:RK,$
      data_reg:data_reg,residuals:residuals,chisq:chisq,data_cont_t:data_cont_t,$
      reg_tweak:reg_tweak,guess:guess,order:order}

  endelse

  return,reg_solution

  if keyword_set(debug) then stop

end

