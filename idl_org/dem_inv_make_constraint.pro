;+
;
; NAME:
;  dem_ inv_make_constraint
;
; PURPOSE:
;   to make contraint matrix
;
;
; CALLING SeqUENCE:
;  dem_inv_make_constraint,L,logT,dlogT,dem_guess,order
;
; CALLS:
;   none
;
; INPUTS:
;   logT    			- array, middle temperature bin
;   dlogT   		 	- array, width of temperature bins
;   dem_guess 		- array, DEM guess
;   order 			- scalar, oder of regularization
;
; OPTIONAL INPUTS:
;   none
;
; OUTPUTS:
;    L -  constraint matrix
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
;  second order regularization works somewhat worse than zero/first
;
;
; RESTRICTIONS:
;
;
; MODifICATION HISTORY:
;   eduard(at)astro.gla.ac.uk, 23 May, 2005
;  21-Jul-2011	Program and Variable names changed    IGH
;  28-Apr-2020  Changed any fltarr() to dblarr()  IGH
;

pro dem_inv_make_constraint,L,logT,dlogT,dem_guess,order

  M=n_elements(dlogT)

  ; second order, constraint is the second derivative
  D2=dblarr(M,M)
  for i=0, M-1 do D2(i,i)  =-2.*sqrt(dlogT[i])/sqrt(abs(dem_guess[i])/logT[i]^2)
  for i=0, M-2 do D2(i,i+1)= 1.*sqrt(dlogT[i])/sqrt(abs(dem_guess[i])/logT[i]^2)
  for i=1, M-1 do D2(i,i-1)= 1.*sqrt(dlogT[i])/sqrt(abs(dem_guess[i])/logT[i]^2)

  ; first order, constraint is the first derivative
  D1=dblarr(M,M)
  for i=0, M-2 do D1(i,i+1)=sqrt(dlogT[i])/sqrt(abs(dem_guess[i])/logT[i])
  for i=0, M-1 do D1(i,i)  =-sqrt(dlogT[i])/sqrt(abs(dem_guess[i])/logT[i])

  ; zeroth order, constraint is the identity matrix
  D0=dblarr(M,M)
  For i=0, M-1 do D0(i,i)=sqrt(dlogT[i])/sqrt(abs(dem_guess[i]))


  if (order ne 0) and (order ne 1) and (order ne 2) then  begin
    order=0
    message,'WARNING: order should be 0,1,or 2 only.regularization order set to zero !'
  end

  if (order eq 0) then L=D0
  if (order eq 1) then L=D1
  if (order eq 2) then L=D2

end


