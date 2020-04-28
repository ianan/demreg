;+
;
; NAME:
;  dem_inv_confidence_interval
;
;
; PURPOSE:
;   Calculates a confidence interval - 1 sigma uncertanty on
;   regularised solution (vertical error)
;
;
; CALLING SEQUENCE:
;
;  dem_inv_confidence_interval,reg_sol,data,edata,Alpha,Betta,U,W,opt,dem_guess,Guess,Npass,reg_sol_err
;
; CALLS:
;
;   dem_inv_reg_solution.pro
;
; INPUTS:
;
;   reg_sol 		- vector, regularised solution
;   data 			 - vector, filter or line data
;   edata			 - vector, 1 sigma uncertanty on data
;
;  Results of Generalised Singular Value Decomposition:
;   alpha  			- diagonal elements of SA
;   betta  			- diagonal elements of SB
;   U,W    			- matrixes consisting of decomposition products
;
;   opt    			- scalar, regularisation parameter
;   dem_gues 		- vector, guess DEM
;   Guess  		- scalar (0 or 1.), use or not to use guess DEM in the constraint
;   Npass  		- scalar, number of Monte Carlo simulations to produce confidence interval
;
; OPTIONAL INPUTS:
;   none
;
; OUTPUTS:
;
;  reg_sol_err		- vector, 1 sigma uncertanty on regularised solution
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
;   eduard(at)astro.gla.ac.uk, 30 September, 2005
;  21-Jul-2011	Program and Variable names changed    IGH
;  28-Apr-2020  Changed any fltarr() to dblarr()    IGH
;-


pro dem_inv_confidence_interval,reg_sol,data,edata,Alpha,Betta,U,W,opt,dem_guess,Guess,Npass,reg_sol_err

  seed = 1001L
  ; Initial seed for a repeatable sequence

  Strip_array=dblarr(Npass,n_elements(reg_sol))

  for i=0, Npass-1 do begin
    rand_arr=(randomu(seed,n_elements(reg_sol))-0.5D0)*2.*edata
    dataX=data+rand_arr
    dem_inv_reg_solution,Alpha,Betta,U,W,dataX,opt,dem_guess*Guess,Reg_solX
    Strip_array[i,*]=reg_solX
  end

  reg_sol_err_plus=dblarr(n_elements(reg_sol))
  reg_sol_err_minus=dblarr(n_elements(reg_sol))
  for j=0, n_elements(reg_sol)-1 do begin
    reg_sol_err_plus[j]=max((strip_array[*,j]-reg_sol[j]))
    reg_sol_err_minus[j]=max(-(strip_array[*,j]-reg_sol[j]))
  end

  reg_sol_err=dblarr(n_elements(reg_sol))
  for j=0, n_elements(reg_sol)-1 do reg_sol_err[j]=max(abs(strip_array[*,j]-reg_sol[j]))

  ;*****************************************************************************

end
