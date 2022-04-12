
; NAME:
;   dem_inv_reg_parameter
;
; PURPOSE:
;   to compute regularization parameter
;
; CALLING SEQUENCE:
;
;  dem_inv_reg_parameter,sigmaA,SigmaB,U,W,Data,Err,dem_guess,reg_tweak,opt
;
; CALLS:
;   none
;
; INPUTS:
;  SigmaA 	- vector, generalised singular values
;  SigmaB 	- vector, generalised singular values
;  U      		- matrix, GSVD matrix
;  W      		- matrix, GSVD matrix
;  Data   		- vector, containing data (eg dn)
;  Err    		- vector, uncertanty on data (same units and dimension)
;  dem_guess		-vector, guess dem solution
; reg_tweak	-scalar, parameter to adjusting regularization (chisq)
;
; OPTIONAL INPUTS:
;   none
;
; OUTPUTS:
;    opt -  regularization parameter
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
;   eduard(at)astro.gla.ac.uk, 23 May, 2005
;   16- Sept-2005: now plots to display and eps files
;  21-Jul-2011	Program and Variable names changed    IGH
;  21-Jul-2011	Commented out plotting of picard condition IGH
;  12-Apr-2022  Commented out un-needed arg2, SV, data_u and c
;-


pro dem_inv_reg_parameter,sigmaA,SigmaB,U,W,Data,Err,dem_guess,reg_tweak,opt
  ;calculates regularisation parameter

  Nmu=9900

  Ndata=n_elements(Data)
  Nreg =n_elements(SigmaA)

  arg=dblarr(Nreg,Nmu)
  arg2=dblarr(Nreg,Nmu)
  discr=dblarr(Nmu)

  maxx=max(SigmaA/SigmaB)*1d3
  minx=max(SigmaA/SigmaB)*1d-15

  step=(alog(maxx)-alog(minx))/(Nmu-1.)
  mu=exp(findgen(Nmu)*step)*minx

  omega=invert(w)##dem_guess

  for k=0,Ndata-1 do begin
    coef=data##u[k,*]-SigmaA[k]*omega[k]
    for i=0,Nmu-1 do begin
      arg[k,i]=(mu[i]*SigmaB[k]*SigmaB[k]*coef/(SigmaA[k]*SigmaA[k]+mu[i]*SigmaB[k]*SigmaB[k]))^2
    end
;    arg2[k,*]=sqrt(arg[k,*])/Err[k]
  end
  discr=total(arg,1)-total(err*err)*reg_tweak

  minimum=min(abs(discr),Min_index)
  opt=mu[Min_index]

;  SV=abs(sigmaA)/abs(SigmaA^2+opt*SigmaB^2)
;  Data_U=abs(data##u)
;  C=data_u*SV

  ;****************************************************************************************
  ;window,4,xsize=600,ysize=600,title='Regularization parameter and Picard condition'
  ;!P.MULTI=[0,1,2]
  ;plot_oo, mu, abs(discr), xrange=[minx,maxx],psym=1,Title='Tikhonov regularization',$
  ;ytitle='Descripancy',xtitle='Regularization parameter, !4k!3'
  ;plot,findgen(Nreg)+1.,C, xrange=[1,Ndata],xstyle=1,/xlog,Psym=10,$
  ;Title='Picard Condition',xtitle='Singular Value number',/ylog
  ;oplot,findgen(Nreg)+1.,data_u,line=1
  ;oplot,findgen(Nreg)+1.,SV,line=2
  ;!P.MULTI=0
  ;****************************************************************************************

  print, 'Regularization parameter (discrepancy): ', opt
  ;****************************************************************************************


end

