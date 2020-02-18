function dem2dn, logT,dem, Tresp

  ; for a given DEM, Temperature response function and temperature binning
  ; return the dn/s per filter and the contribution from each temperature bin

  ; Input:
  ; DEM [cm-5 K-1]
  ; logT [K]
  ; Tresp [DN cm^5 s^-1 pix^-1]

  ;Output
  ; logT [k]
  ; dn [DN/s/px]
  ; TC [DN/s/px]

  ; All should use the same temperature binning
  ; Should work with any Tresp as long and units correct and T binning the same
  ; between all input parameters

  ; 15/09/2010 IGH started

  nt=n_elements(logT)
  nf=n_elements(Tresp[0,*])
  TC=dblarr(nt,nf)
  dn=dblarr(nf)

  dlogT  =logT[1:nT-1]-logT[0:nT-2]
  dlogT=[dlogT,dlogT[nt-2]]

  good=indgen(nt)
  ngd=nt

  if (ngd[0] gt 0) then for f=0, nf -1 do TC[good,f]=dem[good]*Tresp[good,f]*10d^logT[good]*alog(10d^dlogt[good])

  dn=total(TC,1)

  result={logt:logt,dn:dn, TC:TC}

  return, result

end