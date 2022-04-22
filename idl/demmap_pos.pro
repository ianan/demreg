pro demmap_pos,dd,ed,rmatrix,logt,dlogt,glc,dem,chisq,$
  edem,elogt,dn_reg,reg_tweak=reg_tweak,max_iter=max_iter,rgt_fact=rgt_fact,$
  dem_norm0=dem_norm0

  ; This is an updated/optimised/bug fixed version of demmap_pos.pro that was included in
  ; the AIA map specific version of the Regularized DEM maps code
  ; http://www.astro.gla.ac.uk/~iain/demreg/map/ and Hannah & Kontar 2013 A&A 553
  ;
  ; The original version is the Regularized DEM code from
  ; Hannah & Kontar 2012 A&A 539
  ; http://www.astro.gla.ac.uk/~iain/demreg or
  ; https://github.com/ianan/demreg/tree/master/idl_org
  ; Which provides more options (order of constraint matrix, any data type, guess solution) but was slower than the
  ; AIA map version.
  ;
  ; This version is trying to maintain the speed of the AIA map version but give some of the features
  ; of the original non-map version e.g. regarding different data types and the guess solution
  ;
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;
  ;  We are trying to get the DEM from the response functions (K), and the data (g) as related via
  ;
  ;     g=K.DEM
  ;
  ;  (Note that the vector and matrix multiplcation is a bit messy in this explanation, but correct in the code and paper)
  ;
  ;  Regularized approach solves this via
  ;
  ;     ||K.DEM-g||^2 + lamb ||L.DEM||^2=min
  ;
  ;  where the extra bits are the constraint matrix (L) and regularization parameter (lamb).
  ;  L is taken as a "zeroth order" constraint, something like diag(L)=sqrt(dLogT)/sqrt(dem_guess)
  ;  As we might not have an initial dem_guess solution we can find one by running the regularization
  ;  using diag(L)=1/sqrt(dLogT) and use the output as dem_reg to make a new L and then run the
  ;  regularization a second time
  ;
  ;  The actual regularization is solved via GSVD of K and L.
  ;  This outputs singular values sva and svb (with sva^2+svb^2=1) and vectors u, v, w
  ;  with properties U^T K W=diag(sva) and V^T L W=diag(svb)
  ;
  ;  The DEM solution is then given by
  ;
  ;     DEM_lamb = Sum_i (sva_i/(sva_i^2+svb_i^1*lamb)) * (g.u) w
  ;
  ;     or
  ;
  ;     K^-1=K^dag= Sum_i (sva_i/(sva_i^2+svb_i^1*lamb)) * u.w
  ;
  ;  We know all the bits of it apart from lamb. We get this from the Discrepancy principle (Morozon, 1967)
  ;  such that the lamb chosen gives a DEM_lamb that produces a specified reduced chisq in data space which
  ;  we call the "regularization parameter" (or reg_tweak) and we normally take this to be 1. As we also want a
  ;  physically real solution (e.g. a DEM_lamb that is positive) we iteratively increase reg_tweak until a
  ;  positive solution is found (or a max number of iterations is reached).
  ;
  ;  Once the solution is found we work out the uncertainies: vertical (DEM uncertainty) is a linear propagation
  ;  of errors on DN through solution; horizontal (T resolution) is how much K^dag#K deviates from I, so measuring
  ;  spread from diagonal but also if regularization failed at that T.
  ;
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Each child process called by dn2dem_map which actually does the DEM calculation
  ; 25-May-2012 IGH - Created
  ; 03-Jul-2012 IGH - Changed DEM weighting as to not used smoothed & sqrt version
  ; 03-Jul-2012 IGH - Now returns dn_reg to be used to produce DN_reg maps
  ; 18-Nov-2012 IGH - Updated version of demmap_broen.pro but works with 2D (px vs ft)
  ;                                     instead of 3D (x,y,FT) input from fn2dem_map2.pro
  ; 20-Jun-2013 IGH - Added reg_tweak (desired chisq of solution) as optional input
  ;
  ; 20-Jun-2013 IGH - New version based on demmap_broen2.pro
  ;                  Iterates until gets posive dem (or gives up) by increasing reg_tweak (chisq)
  ;                  Optional inputs:
  ;                         max_iter (number of iteration to try before giving up)
  ;                         rgt_fact (factor to increase reg_tweak by each iteration)
  ;                         [so max reg_tweak used would be reg_tweak*rgt_fact^max_iter]
  ; 10-Mar-2015 IGH - Changed from demmap_broen_pos so can be tweaked with dn2dem_po_nb.pro
  ;
  ;                  Need to input logt and dlogt
  ;                  The initial dem_reg guess calculation in previous versions can create NaNs
  ;                     - Work in progress to fix it
  ;
  ;                  Removed the first calc of dem_reg from the pos loop
  ;                      - No effect on code output other than a speed up
  ;
  ;                  Calculate Lorg then initial dem_reg or gloci and initial dem_reg in here now
  ;                     - don't pass in Lorg anymore
  ;                     - if doing gloci do it using all filters or just the selected via glc ne 0
  ;
  ; 14-Apr-2016 IGH - Corrected bug with wrong dem_reg (should be dem_reg_out) being used to calculate dn_reg and chisq
  ; 25-Apr-2016 IGH - Updated some of the internal variable names and increased comments (though more to do!)
  ; 27-Apr-2016 IGH - Added in option to supply initial guess/constraint normalized DEM to weight L
  ; 19-May-2016 IGH - Added in check for dem_norm0, if supplied but any <=0 then ignore
  ;                       Also tweaked testing that data in all filters >0 via product()
  ; 02-Aug-2016 IGH - Renamed variable using dem_norm to dem_reg_wght to avoid bug/conflict with solution dem_reg
  ; 20-Sep-2019 IGH - Added check for Nan/Inf as well as 0s before doing the calculation
  ; 22-Apr-2022 IGH - Moved final L calc and GSVD out of the pos loop, faster if just calc before
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  na=n_elements(dd[*,0])
  nf=n_elements(RMatrix[0,*])
  nt=n_elements(logt)

  dem=dblarr(na,nt)
  edem=dblarr(na,nt)
  elogt=dblarr(na,nt)
  RMatrixin=dblarr(nt,nf)
  kdag=dblarr(nf,nt)
  filter=dblarr(nf,nt)
  chisq=dblarr(na)
  kdagk=dblarr(nt,nt)
  dn_reg=dblarr(na,nf)

  ltt=min(logt)+(max(logt)-min(logt))*findgen(51)/(51-1.0)

  if (n_elements(reg_tweak) ne 1) then reg_tweak=1
  if (n_elements(max_iter) ne 1) then max_iter=10
  if (n_elements(rgt_fact) ne 1) then rgt_fact=1.5

  nmu=42 ; but of course
  for i=0, na-1 do begin
    dnin=reform(dd[i,*])
    ednin=reform(ed[i,*])
    if ((size(dem_norm0))[0] gt 0) then dem_reg_wght=reform(dem_norm0[i,*])
    for kk=0,nf-1 do RMatrixin[*,kk]=RMatrix[*,kk]/eDNin[kk]

    dn=dnin/ednin
    edn=ednin/ednin

    ; Test if any of the data is 0, if so can just ignore this one
    ;    if (product(dn) gt 0.) then begin
    ; Added extra text in case any of the data is NaNs or Infs
    if (product(dn,/nan) gt 0. and product(finite(dn)) gt 0.) then begin
      ; reset the positive check
      ndem=1
      piter=0
      rgt=reg_tweak

      L=dblarr(nT,nT)
      ; If you have supplied an initial guess/constraint normalized DEM then don't
      ; need to calculate one (either from L=1/sqrt(dLogT) or min of EM loci)

      ; Though need to check what you have supplied is correct dimension
      ; and no element 0 or less.
      test_dem_reg=0

      if (n_elements(dem_reg_wght) eq nt) then begin
        if (product(dem_reg_wght) gt 0) then test_dem_reg=1
      endif

      if (test_dem_reg eq 0) then begin
        if (total(glc) gt 0.) then begin
          ; use the min of the emloci as the initial dem_reg
          gdglc=where(glc gt 0,ngdglc)
          emloci=dblarr(nt,ngdglc)
          for ee=0, ngdglc-1 do emloci[*,ee]=dnin[gdglc[ee]]/(RMatrix[*,gdglc[ee]])
          dem_model=dblarr(nt)
          for ttt=0, nt-1 do dem_model[ttt]=min(emloci[ttt,*] > 0.)
          dem_model=smooth(dem_model,3)
          dem_reg=dem_model/max(dem_model)
          dem_reg=dem_reg*(dem_reg gt 0.) +1d-10
        endif else begin
          ; Calculate the initial constraint matrix
          ; Just a diagional matrix scaled by dlogT
          for gg=0, nt-1 do L[gg,gg]=1.0/sqrt(dlogT[gg])
          ; Better to use dT not dlogT - probably not from synthetic tests.
          ;  for gg=0, nt-1 do L[gg,gg]=1.0/sqrt((10^(logt[gg]+0.5*dlogt[gg])-10^(logt[gg]-0.5*dlogt[gg])))

          ;################ Work out the 1st DEM_reg ###########################
          dem_inv_gsvdcsq,RMatrixin,L,sva,svb,U,V,W
          dem_inv_reg_parameter_map,sva,svb,U,W,DN,eDN,rgt,lamb,nmu
          for kk=0, nf-1 do filter[kk,kk]=sva[kk]/(sva[kk]*sva[kk]+$
            svb[kk]*svb[kk]*lamb)
          kdag=W##matrix_multiply(U[0:nf-1,0:nf-1],filter,/atrans)
          dr0=reform(kdag##dn)

          ; only take the positive with ceratin amount (fcofmx) of max, then make rest small positive
          fcofmx=1d-4
          dem_reg=dr0*(dr0 gt 0 and dr0 gt fcofmx*max(dr0))+1*(dr0 lt 0 or dr0 lt fcofmx*max(dr0))
          dem_reg=dem_reg/(fcofmx*max(dr0))
          ;  Don't need the smoothed version anymore? - seems to help with synthetic tests
          dem_reg=smooth(dem_reg,3)
        endelse
      endif else begin
        dem_reg=dem_reg_wght
      endelse
      ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      ; Faster and similar results in just do the GSVD on R and L here
      ; And not updated L each time
      for kk=0, nt-1 do L[kk,kk]=sqrt(dlogT[kk])/sqrt(abs(dem_reg[kk]))
      dem_inv_gsvdcsq,RMatrixin,L,sva,svb,U,V,W

      ; ######## Still don't have a positive DEM_reg (or reached max_iter?) ########
      while(ndem gt 0 and piter lt max_iter) do begin

        ;################ Use first DEM_reg to weight L ###########################
        ;################ or from last loop   ######################################
;        for kk=0, nt-1 do L[kk,kk]=sqrt(dlogT[kk])/sqrt(abs(dem_reg[kk]))
;        dem_inv_gsvdcsq,RMatrixin,L,sva,svb,U,V,W
        ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dem_inv_reg_parameter_map,sva,svb,U,W,DN,eDN,rgt,lamb,nmu

        ;################ Work out the inverse of K (Rmatrixin) ####################
        for kk=0, nf-1 do filter[kk,kk]=sva[kk]/(sva[kk]*sva[kk]+svb[kk]*svb[kk]*lamb)
        kdag=W##matrix_multiply(U[0:nf-1,0:nf-1],filter,/atrans)

        ;################ Work out the final DEM_reg ####################

        DEM_reg_out=reform(kdag##dn)

        ; Is any of this solution negative?????
        nn=where(DEM_reg_out lt 0, ndem)
        rgt=rgt_fact*rgt
        piter=piter+1

        ;       IS RECALCULATING DEM_REG FOR THE L AS THE LOOP STARTS AGAIN IMPORTANT?
        ;       IT SEEMS THAT IT IS ALL JUST DOWN TO INCREASING RGT TO GET A +VE DEM....
        ;       So having the below in doesn't do much other than slow the code down
        ;        ; just in case we need dem_reg for the next loop and a new L
        ;        ; only take the positive with ceratin amount (fcofmx) of max, then make rest small positive
        ;        fcofmx=1d-4
        ;        dr0=dem_reg_out
        ;        dem_reg=dr0*(dr0 gt 0 and dr0 gt fcofmx*max(dr0))+1*(dr0 lt 0 or dr0 lt fcofmx*max(dr0))
        ;        dem_reg=dem_reg/(fcofmx*max(dr0))
        ;        ;  Don't need the smoothed version anymore - seems to help with synthetic tests
        ;        dem_reg=smooth(dem_reg,3)


      endwhile
      ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ;############ if positive or reached max_iter work rest out ###########

      dem[i,*]=DEM_reg_out

      ; Work out DN and residuals the DEM solution gives
      dn_reg0=reform(rmatrix##DEM_reg_out)
      dn_reg[i,*]=dn_reg0
      residuals=(dnin-dn_reg0)/ednin
      chisq[i]=total(residuals^2)/(nf)


      ;################ Do the error calcualtion ######################

      ; Error matrix for eDEM
      delxi2=matrix_multiply(kdag,kdag,/atrans)
      edem[i,*]=sqrt(diag_matrix(delxi2))

      ; Resolution matrix for elogt
      ; If everything worked then kdag so be the perfect inverse of k
      ; so kdag#k=I would just be a diagional matrix
      kdagk=kdag##rmatrixin

      for kk=0, nt-1 do begin
        rr=interpol(transpose(kdagk[kk,*]),logt,ltt)
        ; Where is the row bigger than the half maximum
        hm_index=where(rr ge max(kdagk[kk,*])/2.,nhmi)
        elogt[i,kk]=dlogt[kk]
        if (nhmi gt 1) then $
          elogt[i,kk]=(ltt[max(hm_index)]-ltt[min(hm_index)]) >dlogt[kk]
      endfor
    endif
    if ((i mod 5000) eq 0)  then print,string(i,format='(i7)')+' of '+string(na,format='(i7)') +$
      ', '+string((i*100./na*1.),format='(f5.1)')+'%'
  endfor
  print,string(100,format='(i3)')+'% Done'

end
