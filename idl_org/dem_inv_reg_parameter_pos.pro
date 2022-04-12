
; NAME:
;   dem_inv_reg_parameter_pos
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
;  26-Sep-2011	This modified version produces an extra positive solution IGH
;  28-Apr-2020  Changed any fltarr() to dblarr() IGH
;  12-Apr-2022  Commented out un-needed end calc (from testing)
;-

pro dem_inv_reg_parameter_pos,sigmaA,SigmaB,U,W,Data,Err,dem_guess,reg_tweak,opt,reg,opt_pos,reg_pos
  ;calculates regularisation parameter

  Nmu=1000.

  Ndata=n_elements(Data)
  Nreg =n_elements(SigmaA)

  arg=dblarr(Nreg,Nmu)
  ;arg2=dblarr(Nreg,Nmu)
  discr=dblarr(Nmu)

  xi=dblarr(nreg,nmu)
  xi_np=dblarr(nmu)

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
  end

  for i=0, nmu-1 do begin
    dem_inv_reg_solution,sigmaA,sigmaB,U,W,data,mu[i],transpose(dem_guess),regout
    xi[*,i]=regout
    ii=where(regout lt 0.,nii)
    xi_np[i]=nii/(nreg*1.0)
  endfor

  discr=total(arg,1)-total(err*err)*reg_tweak

  pos=where(xi_np eq 0, npos)
  reg_pos=dblarr(nreg)
  if (npos gt 0) then begin
    minimum=min(abs(discr[pos]),Min_index)
    opt_pos=mu[pos[Min_index]]
    reg_pos=xi[*,pos[min_index]]
  endif
  minimum=min(abs(discr),Min_index)
  opt=mu[Min_index]
  reg=xi[*,min_index]

  ; SV=abs(sigmaA)/abs(SigmaA^2+opt*SigmaB^2)
  ; Data_U=abs(data##u)
  ; C=data_u*SV

;  ar=dblarr(n_elements(data),n_elements(SigmaA))
;  arf=dblarr(n_elements(data),n_elements(SigmaA))
;  fact=(sigmaA/(sigmaA*sigmaA+opt*sigmaB*sigmaB))[0:5]
;  fact_pos=(sigmaA/(sigmaA*sigmaA+opt_pos*sigmaB*sigmaB))[0:5]
;
;  for k=0,n_elements(data)-1 do begin
;    scal =data##u[k,*]
;    for j=0,n_elements(SigmaA)-1 do begin
;      ar[k,j]=scal*w[k,j]
;      arf[k,j]=(sigmaA[k]*scal*w[k,j])/(sigmaA[k]*sigmaA[k]+opt*sigmaB[k]*sigmaB[k])
;    end
;  end


  print, 'Regularization parameter (discrepancy and pos): ', opt, opt_pos
  ;****************************************************************************************


end

