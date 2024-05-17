import numpy as np
from numpy import diag
from dem_inv_gsvd import dem_inv_gsvd
from dem_reg_map import dem_reg_map
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from tqdm import tqdm
from threadpoolctl import threadpool_limits

def demmap_pos(dd,ed,rmatrix,logt,dlogt,glc,reg_tweak=1.0,max_iter=10,rgt_fact=1.5,dem_norm0=None,nmu=42,warn=False,l_emd=False,rscl=False):
    """
    demmap_pos
    computes the dems for a 1 d array of length na with nf filters using the dn (g) counts and the temperature
    response matrix (K) for each filter.
    where 

        g=K.DEM

    Regularized approach solves this via
  
        ||K.DEM-g||^2 + lamb ||L.DEM||^2=min

    L is a zeroth order constraint matrix and lamb is the rrgularisation parameter

    The regularisation is solved via the GSVD of K and L (using dem_inv_gsvd)
    which provides the singular values (sva,svb) and the vectors u,v and w
    witht he properties U.T*K*W=sva*I and V.T L W = svb*I

    The dem is then obtained by:

        DEM_lamb = Sum_i (sva_i/(sva_i^2+svb_i^1*lamb)) * (g.u) w

    or

        K^-1=K^dag= Sum_i (sva_i/(sva_i^2+svb_i^1*lamb)) * u.w    

    We know all the bits of it apart from lamb. We get this from the Discrepancy principle (Morozon, 1967)
    such that the lamb chosen gives a DEM_lamb that produces a specified reduced chisq in data space which
    we call the "regularization parameter" (or reg_tweak) and we normally take this to be 1. As we also want a
    physically real solution (e.g. a DEM_lamb that is positive) we iteratively increase reg_tweak until a
    positive solution is found (or a max number of iterations is reached).

    Once a solution that satisfies this requirement is found the uncertainties are worked out:
    the vertical errors (on the DEM) are obtained by propagation of errors on dn through the
    solution; horizontal (T resolution) is how much K^dag#K deviates from I, so measuring
    spread from diagonal but also if regularization failed at that T.

    Inputs

    dd
        the dn counts for each channel
    ed
        the error on the dn counts
    rmatrix
        the trmatrix for each channel 
    logt
        log of the temperature bin averages
    dlogt
        size of the temperature bins
    glc
        indices of the filters for which gloci curves should be used to set the initial L constraint
        (if called from dn2dem_pos, then all 1s or 0s)

    Optional inputs

    reg_tweak
        initial chisq target
    rgt_fact
        scale factor for the increase in chi-sqaured target for each iteration
    max_iter
        maximum number of times to attempt the gsvd before giving up, returns the last attempt if max_iter reached
    dem_norm0
        provides a "guess" dem as a starting point, if none is supplied one is created.
    nmu
        number of reg param samples to use
    warn
        print out warnings
    l_emd
        remove sqrt from constraint matrix (best with EMD)
    
    Outputs

    
    dem
        The DEM(T)
    edem
        the error on the DEM(T)
    elogt
        the error on logt    
    chisq
        the chisq for the dem compared to the dn
    dn_reg
        the simulated dn for each filter for the recovered DEM    
    """
 
    na=dd.shape[0]
    nf=rmatrix.shape[1]
    nt=logt.shape[0]
    #set up some arrays
    dem=np.zeros([na,nt])
    edem=np.zeros([na,nt])
    elogt=np.zeros([na,nt])
    rmatrixin=np.zeros([nt,nf])
    kdag=np.zeros([nf,nt])
    filt=np.zeros([nf,nt])
    chisq=np.zeros([na])
    kdagk=np.zeros([nt,nt])
    dn_reg=np.zeros([na,nf])
    ednin=np.zeros([nf])
 
    # do we have enough DEM's to make parallel make sense?
    if (na>=200):
        n_par = 100
        niter=(int(np.floor((na)/n_par)))
#       Put this here to make sure running dem calc in parallel, not the underlying np/gsvd stuff (this correct/needed?)  
        with threadpool_limits(limits=1):
            with ProcessPoolExecutor() as exe:
                futures=[exe.submit(dem_unwrap, dd[i*n_par:(i+1)*n_par,:],ed[i*n_par:(i+1)*n_par,:],rmatrix,logt,dlogt,glc, \
                    reg_tweak=reg_tweak,max_iter=max_iter,rgt_fact=rgt_fact,dem_norm0=dem_norm0[i*n_par:(i+1)*n_par,:],\
                        nmu=nmu,warn=warn,l_emd=l_emd,rscl=rscl) \
                        for i in np.arange(niter)]
                kwargs = {
                    'total': len(futures),
                    'unit': ' x10^2 DEM',
                    'unit_scale': True,
                    'leave': True
                    }
                for f in tqdm(as_completed(futures), **kwargs):
                    pass
            for i,f in enumerate(futures):
            #store the outputs in arrays
                dem[i*n_par:(i+1)*n_par,:]=f.result()[0]
                edem[i*n_par:(i+1)*n_par,:]=f.result()[1]
                elogt[i*n_par:(i+1)*n_par,:]=f.result()[2]
                chisq[i*n_par:(i+1)*n_par]=f.result()[3]
                dn_reg[i*n_par:(i+1)*n_par,:]=f.result()[4]
            #if there are any remaining dems then execute remainder in serial
            if (np.mod(na,niter*n_par) != 0):
                i_start=niter*n_par
                for i in range(na-i_start):
                    result=dem_pix(dd[i_start+i,:],ed[i_start+i,:],rmatrix,logt,dlogt,glc, \
                        reg_tweak=reg_tweak,max_iter=max_iter,rgt_fact=rgt_fact,dem_norm0=dem_norm0[i_start+i,:],\
                            nmu=nmu,warn=warn,l_emd=l_emd,rscl=rscl)
                    dem[i_start+i,:]=result[0]
                    edem[i_start+i,:]=result[1]
                    elogt[i_start+i,:]=result[2]
                    chisq[i_start+i]=result[3]
                    dn_reg[i_start+i,:]=result[4]
        
    #else we execute in serial
    else:   
        for i in range(na):
            result=dem_pix(dd[i,:],ed[i,:],rmatrix,logt,dlogt,glc, \
                reg_tweak=reg_tweak,max_iter=max_iter,rgt_fact=rgt_fact,dem_norm0=dem_norm0[i,:],nmu=nmu,warn=warn,l_emd=l_emd,rscl=rscl)
            dem[i,:]=result[0]
            edem[i,:]=result[1]
            elogt[i,:]=result[2]
            chisq[i]=result[3]
            dn_reg[i,:]=result[4]

    return dem,edem,elogt,chisq,dn_reg

def dem_unwrap(dn,ed,rmatrix,logt,dlogt,glc,reg_tweak=1.0,max_iter=10,rgt_fact=1.5,dem_norm0=0,nmu=42,warn=False,l_emd=False,rscl=False):
    #this nasty function serialises the parallel blocks
    ndem=dn.shape[0]
    nt=logt.shape[0]
    nf=dn.shape[1]
    dem=np.zeros([ndem,nt])
    edem=np.zeros([ndem,nt])
    elogt=np.zeros([ndem,nt])
    chisq=np.zeros([ndem])
    dn_reg=np.zeros([ndem,nf])
    for i in range(ndem):
        result=dem_pix(dn[i,:],ed[i,:],rmatrix,logt,dlogt,glc, \
                reg_tweak=reg_tweak,max_iter=max_iter,rgt_fact=rgt_fact,dem_norm0=dem_norm0[i,:],nmu=nmu,warn=warn,l_emd=l_emd,rscl=rscl)
        dem[i,:]=result[0]
        edem[i,:]=result[1]
        elogt[i,:]=result[2]
        chisq[i]=result[3]
        dn_reg[i,:]=result[4]
    return dem,edem,elogt,chisq,dn_reg

def dem_pix(dnin,ednin,rmatrix,logt,dlogt,glc,reg_tweak=1.0,max_iter=10,rgt_fact=1.5,dem_norm0=0,nmu=42,warn=True,l_emd=False,rscl=False):
    
    nf=rmatrix.shape[1]
    nt=logt.shape[0]
    # nmu=42
    ltt=min(logt)+1e-8+(max(logt)-min(logt))*np.arange(51)/(52-1.0)
    dem=np.zeros(nt)
    edem=np.zeros(nt)
    elogt=np.zeros(nt)
    chisq=0
    dn_reg=np.zeros(nf)

    rmatrixin=np.zeros([nt,nf])
    filt=np.zeros([nf,nt])  
    
    for kk in np.arange(nf):
        #response matrix
        rmatrixin[:,kk]=rmatrix[:,kk]/ednin[kk]
    dn=dnin/ednin
    edn=ednin/ednin

    # checking for Inf and NaN
    if ( sum(np.isnan(dn)) == 0 and sum(np.isinf(dn)) == 0 and np.prod(dn) > 0):
        ndem=1
        piter=0
        rgt=reg_tweak

        L=np.zeros([nt,nt])

        test_dem_reg=(np.zeros(1)).astype(int)
        
#  If you have supplied an initial guess/constraint normalized DEM then don't
#  need to calculate one (either from L=1/sqrt(dLogT) or min of EM loci)

# As the call to this now sets dem_norm to array of 1s if nothing provided by user can also test for that

#     Before calling this dem_norm0 is set to array of 1s if nothing provided by user
#     So we need to work out some weighting for L or is one provided as dem_norm0 (not 0 or array of 1s)?
        if (np.prod(dem_norm0) == 1.0 or dem_norm0[0] == 0):
# Need to work out a weighting here then, have two appraoches:
#         1. Do it via the min of em loci - chooses this if gloci, glc=1 from user
            if (np.sum(glc) > 0.0):
                gdglc=(glc>0).nonzero()[0]
                emloci=np.zeros((nt,gdglc.shape[0]))
                #for each gloci take the minimum and work out the emission measure
                for ee in np.arange(gdglc.shape[0]):
                    emloci[:,ee]=dnin[gdglc[ee]]/(rmatrix[:,gdglc[ee]])
                #for each temp we take the min of the loci curves as the estimate of the dem
                dem_model=np.zeros(nt)
                for ttt in np.arange(nt):
                    dem_model[ttt]=np.min(emloci[ttt,np.nonzero(emloci[ttt,:])])
                dem_reg_lwght=dem_model
#                ~~~~~~~~~~~~~~~~~ 
#             2. Or if nothing selected will run reg once, and use solution as weighting (self norm appraoch)
            else:
                # Calculate the initial constraint matrix
                # Just a diagional matrix scaled by dlogT
                L=np.diag(1.0/np.sqrt(dlogt[:]))
                #run gsvd
                sva,svb,U,V,W=dem_inv_gsvd(rmatrixin.T,L)
                #run reg map
                lamb=dem_reg_map(sva,svb,U,W,dn,edn,rgt,nmu)
                #filt, diagonal matrix
                for kk in np.arange(nf):
                    filt[kk,kk]=(sva[kk]/(sva[kk]**2+svb[kk]**2*lamb))
                kdag=W@(filt.T@U[:nf,:nf])
                dr0=(kdag@dn).squeeze()
                # only take the positive with certain amount (fcofmx) of max, then make rest small positive
                fcofmax=1e-4
                mask=np.where(dr0 > 0) and (dr0 > fcofmax*np.max(dr0))
                dem_reg_lwght=np.ones(nt)
                dem_reg_lwght[mask]=dr0[mask]
#                ~~~~~~~~~~~~~~~~~ 
#            Just smooth these inital dem_reg_lwght and max sure no value is too small
#             dem_reg_lwght=(np.convolve(dem_reg_lwght,np.ones(3)/3))[1:-1]/np.max(dem_reg_lwght[:])     
            dem_reg_lwght=(np.convolve(dem_reg_lwght[1:-1],np.ones(5)/5))[1:-1]/np.max(dem_reg_lwght[:])       
            dem_reg_lwght[dem_reg_lwght<=1e-8]=1e-8
        else:
#             Otherwise just set dem_reg to inputted weight   
            dem_reg_lwght=dem_norm0

#       Now actually do the dem regularisation using the L weighting from above
#       Faster to do this and the GSVD on R and L before the pos loop 
        if l_emd:
            # this works better with EMD calc, instead of DEM
            L=np.diag(1/abs(dem_reg_lwght)) 
        else:
            L=np.diag(np.sqrt(dlogt)/np.sqrt(abs(dem_reg_lwght))) 
        sva,svb,U,V,W = dem_inv_gsvd(rmatrixin.T,L)
#  Now loop until positive solution or max_iter reached
        while((ndem > 0) and (piter < max_iter)):
            # #make L from 1/dem reg scaled by dlogt and diagonalise
            # L=np.diag(np.sqrt(dlogt)/np.sqrt(abs(dem_reg_lwght))) 
            # #call gsvd and reg map
            # sva,svb,U,V,W = dem_inv_gsvd(rmatrixin.T,L)
            lamb=dem_reg_map(sva,svb,U,W,dn,edn,rgt,nmu)
            for kk in np.arange(nf):
                filt[kk,kk]=(sva[kk]/(sva[kk]**2+svb[kk]**2*lamb))
            kdag=W@(filt.T@U[:nf,:nf])
    
            dem_reg_out=(kdag@dn).squeeze()

            ndem=len(dem_reg_out[dem_reg_out < 0])
            rgt=rgt_fact*rgt
            piter+=1
        
        if (warn and (piter == max_iter)):
            print('Warning, positivity loop hit max iterations, so increase max_iter? Or rgt_fact too small?')
      
        dem=dem_reg_out

        #work out the theoretical dn and compare to the input dn
        dn_reg=(rmatrix.T @ dem_reg_out).squeeze()
        residuals=(dnin-dn_reg)/ednin
        #work out the chisquared
        chisq=np.sum(residuals**2)/(nf)

        #do error calculations on dem
        delxi2=kdag@kdag.T
        edem=np.sqrt(np.diag(delxi2))

        kdagk=kdag@rmatrixin.T

        elogt=np.zeros(nt)
        for kk in np.arange(nt):
            rr=np.interp(ltt,logt,kdagk[:,kk])               
            hm_mask=(rr >= max(kdagk[:,kk])/2.)
            elogt[kk]=dlogt[kk]
            if (np.sum(hm_mask) > 0):
                elogt[kk]=(ltt[hm_mask][-1]-ltt[hm_mask][0])/2

        if rscl:
            mnrat=np.mean(dnin/dn_reg)
            dem=dem*mnrat
            edem=edem*mnrat
            dn_reg=(rmatrix.T @ dem).squeeze()
            chisq=np.sum(((dnin-dn_reg)/ednin)**2)/nf
    return dem,edem,elogt,chisq,dn_reg