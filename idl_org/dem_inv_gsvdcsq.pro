;+
; PROJECT:
;   RHESSI/INVERSION
;
; NAME:
;  dem_inv_gsvdcsq
;
; PURPOSE:
;   To perform Generalised Singular Value Decomposition of two matrixes with
;   double precision math
;
;
; CALLING SEQUENCE:
;
;  dem_inv_gsvdcsq, A,B, alpha,betta,U,V,W
;
; CALLS:
;   none
;
; INPUTS:
;
;   A  - cross-section matrix
;   B  - regularization matrix (should be square)
;
; OPTIONAL INPUTS:
;   none
;
; OUTPUTS:

;  Results of Generalised Singular Value Decomposition:
;
; A=U*SA*W^-1
; B=V*SB*W^-1
;
; where SA and SB are diagonal and SA^2+SB^2=1
; matrixes U and V are ortogonal e.g. U*U^T=U^T*U=1, V*V^T=V^T*V=1
;
;   alpha  - diagonal elements of SA
;   betta  - diagonal elements of SB
;   U,V,W  - matrixes consisting of decomposition products
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
;   20- Sept-2005: corrected SVDC routine output - re-ordering SValues and SVectors
;   Note: SVDC produces non-decreasing singular values
;  21-Jul-2011	Program and Variable names changed    IGH
;  22-Apr-2022  Commented out unnecessary OneA calc
;-

pro dem_inv_gsvdcsq, A,B, alpha,betta,U,V,W
  ;; Generalised Singular Value Decomposition with double precision ariphmetics


  AB1=A##invert(B,/DOUBLE)
  ; when B is a square matrix

  SVDC, AB1,S,U,V,/double
  ; using SVDC routine
  Gamma=S
  Gamma=S[REVERSE(SORT(S))]
  U    =U(REVERSE(SORT(S)),*)
  V    =V(REVERSE(SORT(S)),*)
  Sigma=Gamma

  Betta=1.D0/sqrt(1.D0+Sigma*Sigma)
  Alpha=Sigma*betta
  ;singular values

  OneB=dblarr(N_elements(Sigma),N_elements(Sigma))
;  OneA=dblarr(N_elements(Sigma),N_elements(Sigma))

  for i=0, N_elements(Sigma)-1 DO OneB(i,i)=betta(i)
;  for i=0, N_elements(Sigma)-1 DO OneA(i,i)=alpha(i)
  ; scaling

  W=invert(invert(oneB)##transpose(V)##B)

end