;+
;
; NAME:
;
;   dem_inv_reg_solution
;
; PURPOSE:
;
;   calculates regulalarised solution
;

; CALLING SEQUENCE:
;
;  dem_inv_reg_solution,sigmaA,SigmaB,U,W,data,opt,dem_guess,reg_sol
;
; CALLS:
;   none
;
; INPUTS:
;  SigmaA 	- vector, generalised singular values
;  SigmaB 	- vector, generalised singular values
;  U     		 - matrix, GSVD matrix
;  W      		- matrix, GSVD matrix
;  data   		- vector, containing data
;  opt    		- scalar, regularization parameter
;  dem_guess		-vector, optional guess solution
;
; OPTIONAL INPUTS:
;   none
;
; OUTPUTS:
;   reg_sol -  regularization parameter
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
;   eduard(at)astro.gla.ac.uk, 16 Sept, 2005
;  21-Jul-2011	Program and Variable names changed    IGH
;  28-Apr-2020  Changed any fltarr() to dblarr()      IGH


pro dem_inv_reg_solution,sigmaA,SigmaB,U,W,data,opt,dem_guess,reg_sol
  ;regularised solution

  ar=dblarr(n_elements(data),n_elements(sigmaA))

  omega=invert(w)##dem_guess

  for k=0,n_elements(data)-1 do begin
    scal =data##u[k,*]
    for j=0,n_elements(SigmaA)-1 do begin
      ar[k,j]=(sigmaA[k]*scal+opt*omega[k]*SigmaB[k]^2)*w[k,j]/$
        (sigmaA[k]*sigmaA[k]+opt*sigmaB[k]*sigmaB[k])
    end
  end

  reg_sol=total(ar,1)

end
