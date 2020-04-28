
; NAME:
;   dem_inv_reg_resolution
;
; PURPOSE:
;   Calculates resolution of regularised solution (horizontal error)
;
; CALLING SEQUENCE:
;;  deminv_reg_resolution,Alpha,Betta,opt,w,logT,dT,FWHM,cent, RK,fwhm2
;
; CALLS:
;   none
;
; INPUTS:
;
;  Alpha  		- vector, generalised singular values
;  Betta  		- vector, generalised singular values
;  opt    		- regularisation parameter
;  w      		- matrix, generalised singular value decomposition
;  logT     		- temperture bin average energy, [K]
;  dT     		- log temperature bin width, [K]
;
; OPTIONAL INPUTS:
;   none
;
; OUTPUTS:
;   FWHM -  temeprature resolution as full width at half maximum
;			of the resolution matrix for each logT bin
;   FWHM1 -  temeprature resolution as 2* maximum logT distance
;			from digional to half max for each logT bin
;   cent -  the centroid temperature of the resolution matrix
; 			for each logT bins
;   RK - 	the resolution matrix (R_\lambda K) itself
; 			regularization worked where diagional
;			resolution of ith temperature bin is from RK[i,*]
;
;; MODIFICATION HISTORY:
;   eduard(at)astro.gla.ac.uk, 16 Sept, 2005
;  21-Jul-2011	Program and Variable names changed    IGH
;  21-Jul-2011	Rewritten to work in more transparent way IGH
;  17-Aug-2011  Added FWHM2 as double max distance from diagional to half max
;  28-Apr-2020  Changed any fltarr() to dblarr()

pro dem_inv_reg_resolution,Alpha,Betta,opt,w,logT,dT,fwhm, cent,RK,fwhm2

  M=n_elements(Alpha)

  Filter=dblarr(M,M)
  For i=0, M-1 do Filter[i,i]=alpha[i]^2/(alpha[i]^2+opt*betta[i]^2)

  RK =w##(Filter##invert(w))

  FWHM=dblarr(M)
  FWHM2=dblarr(M)
  cent=dblarr(M)

  ltt=min(logt)+(max(logt)-min(logt))*findgen(1001)/(1001-1.0)

  for i=0, M-1 do BEGIN

    rr=interpol(transpose(RK[i,*]),logt,ltt)
    ; where is the row bigger than the half maximum
    hm_index=where(rr ge max(RK[i,*])/2.,nhmi)

    ;if fewer than 2 temperature bins just set FWHM as 1.5dt
    fwhm[i]=1.5*dt[i]
    fwhm2[i]=1.5*dt[i]
    if (nhmi gt 1) then begin
      fwhm[i]=(ltt[max(hm_index)]-ltt[min(hm_index)]) >dt[i]
      fwhm2[i]=2*max([abs(logt[i]-ltt[min(hm_index)]),abs(logt[i]-ltt[max(hm_index)])])
    endif

    ; also stored the assumed centorid
    ; if the regularisation worked at this temperature this shold be the diagional values
    ; i.e. roughly the same logT and this temperature bin
    if (nhmi ge 1 ) then cent[i]=ltt[mean(hm_index)]

  endfor

end



