import numpy as np
from demmap_pos import demmap_pos

def dn2dem_pos(dn_in,edn_in,tresp,tresp_logt,temps,reg_tweak=1.0,max_iter=10,gloci=0,\
    rgt_fact=1.5,dem_norm0=None,nmu=40,warn=False,emd_int=False,emd_ret=False,l_emd=False,non_pos=False,rscl=False):
    """
    Performs a Regularization on solar data, returning the Differential Emission Measure (DEM)
    using the method of Hannah & Kontar A&A 553 2013
    Basically getting DEM(T) out of g(f)=K(f,T)#DEM(T)

    --------------------
    Inputs:
    --------------------

    dn_in:
        The dn counts in dn/px/s for each filter is shape nx*ny*nf (nf=number of filters nx,ny = spatial dimensions, 
        one or both of which can be size 0 for 0d/1d problems. (or nt,nf or nx,nt,nf etc etc to get time series)
    edn_in:
        The error on the dn values in the same units and same dimensions.
    tresp:
        The temperature response matrix size n_tresp by nf
    tresp_logt:
        The temperatures in log t which the temperature response matrix corresponds to. E.G if your tresp matrix 
        runs from 5.0 to 8.0 in steps of 0.05 then this is the input to tresp_logt
    temps:
        The temperatures at which to calculate a DEM, array of length nt.

    --------------------
    Optional Inputs:
    --------------------

    dem_norm0:
        This is an array of length nt which contains an initial guess of the DEM solution providing a weighting 
        for the inversion process (L constraint matrix). The actual values of the normalisation 
        do not matter, only their relative values. 
        If no dem_norm0 given then L weighting based on value of gloci (0 is default)
    gloci:
        If no dem_norm0 given (or dem_norm0 array of 1s) then set gloci 1 or 0 (default 0) to choose weighting for the 
        inversion process (L constraint matrix).
        1: uses the min of EM loci curves to weight L.
        0: uses two reg runs - first with L=diag(1/dT) and DEM result from this used to weight L for second run. 
    reg_tweak:
        The initial normalised chisq target.
    max_iter:
        The maximum number of iterations to attempt, code iterates if negative DEM is produced. If max iter is reached before
        a suitable solution is found then the current solution is returned instead (which may contain negative values)
        (Default is only 10 - although non_pos=True will set as 1)
    rgt_fact:
        The factor by which rgt_tweak increases each iteration. As the target chisq increases there is more flexibility allowed 
        on the DEM
    nmu:
        Number of reg param samples to calculate (default (or <=40) 500 for 0D, 42 for map)
    warn:
        Print out any warnings (always warn for 1D, default no for higher dim data)
    emd_int:
        Do the regularization in EMD [cm^-5] instead of DEM [cm^-5 K^-1] space? (default False). In some circumstances this 
        does seem to help (particularly at higher T), but needs additional tweaking, so why it is not the default.
    emd_ret:
        Return EMD solution instead of EMD [cm^-5] instead of DEM [cm^-5 K^-1] (default False)
    l_emd:
        Remove sqrt factor in constraint matrix, provides better solutions with EMD (and if higher T issues?) 
        (default False, but True with emd_int=True)
    non_pos:
        Return the first solution irrespective of it being positive or not (default False). 
        Done by setting max_iter=1, so user max_iter value ignored
    

    --------------------
    Outputs:
    --------------------

    dem:
        The DEM, has shape nx*ny*nt and units out depends on the input units of tresp and setting of emd_ret
    edem:
        Vertical errors on the DEM, same units as DEM.
    elogt:
        Horizontal errors on temperature, as the name suggests in logT.
    chisq:
        The final chisq, shape nx*ny. Pixels which have undergone more iterations will in general have higher chisq.
    dn_reg:
        The simulated dn counts, shape nx*ny*nf. This is obtained by multiplying the DEM(T) by the filter 
        response K(f,T) for each channel, very important for comparing with the initial data.
    
 
    """
    #create our bin averages:
    logt=([np.mean([(np.log10(temps[i])),np.log10((temps[i+1]))]) for i in np.arange(0,len(temps)-1)])
    #and widths
    dlogt=(np.log10(temps[1:])-np.log10(temps[:-1]))
    nt=len(dlogt)
    logt=(np.array([np.log10(temps[0])+(dlogt[i]*(float(i)+0.5)) for i in np.arange(nt)]))
    #number of DEM entries
    
    #hopefully we can deal with a variety of data, nx,ny,nf
    sze=dn_in.shape

    #for a single pixel
    if (np.any(dem_norm0)==None):
#         If no dem0 wght given just set them all to 1
        dem_norm0=np.ones(np.hstack((dn_in.shape[0:-1],nt)).astype(int))
    if len(sze)==1:
        nx=1
        ny=1
        nf=sze[0]
        dn=np.zeros([1,1,nf])
        dn[0,0,:]=dn_in
        edn=np.zeros([1,1,nf])
        edn[0,0,:]=edn_in
        if (np.all(dem_norm0) != None):
            dem0=np.zeros([1,1,nt])
            dem0[0,0,:]=dem_norm0
        if (warn == False):
            warn=True
        if (nmu <= 40):
            nmu=500
    #for a row of pixels
    if len(sze)==2:
        nx=sze[0]
        ny=1
        nf=sze[1]
        dn=np.zeros([nx,1,nf])
        dn[:,0,:]=dn_in
        edn=np.zeros([nx,1,nf])
        edn[:,0,:]=edn_in
        if (np.all(dem_norm0) != None):
            dem0=np.zeros([nx,1,nt])
            dem0[:,0,:]=dem_norm0
        if (nmu <= 40):
            nmu=42
    #for 2d image
    if len(sze)==3:
        nx=sze[0]
        ny=sze[1]
        nf=sze[2]
        dn=np.zeros([nx,ny,nf])
        dn[:,:,:]=dn_in
        edn=np.zeros([nx,ny,nf])
        edn[:,:,:]=edn_in
        if (np.all(dem_norm0) != None):
            dem0=np.zeros([nx,ny,nt])
            dem0[:,:,:]=dem_norm0
        if (nmu <= 40):
            nmu=42  

    # If want to ignore positivity constraint then set max_iter=1 and no need for the warnings
    if non_pos:
        max_iter=1
        warn=False    

    # If rgt_fact <=1 then the positivity loop wont work so warn about it
    if (warn and (rgt_fact <= 1)):
        print('Warning, rgt_fact should be > 1, for postivity loop to iterate properly.')
    
    # Set glc to either none or all, based on gloci input (default none/not using)
# IDL version of code allows selective use of gloci, i.e [1,1,0,0,1,1] to chose 4 of 6 filters for EM loci
# dem_pix() in demmap_pos.py does allow this, but not sure will work through these wrapper functions
# also not sure if this functionality is actually needed, just stick with all filter or none?
    if gloci == 1:
        glc=np.ones(nf)
        glc.astype(int)
    else:
        glc=np.zeros(nf)
        glc.astype(int)      
        
    if len(tresp[0,:])!=nf:
        print('Tresp needs to be the same number of wavelengths/filters as the data.')
    
    truse=np.zeros([tresp[:,0].shape[0],nf])
    #check the tresp has no elements <0
    #replace any it finds with the mimimum tresp from the same filter
    for i in np.arange(0,nf):
        #keep good TR data
        truse[tresp[:,i] > 0]=tresp[tresp[:,i] > 0]
        #set bad data to the minimum
        truse[tresp[:,i] <= 0,i]=np.min(tresp[tresp[:,i] > 0],axis=0)[i]

    tr=np.zeros([nt,nf])
    for i in np.arange(nf):
#       Ideally should be interp in log-space, so changed
# Not as big an issue for purely AIA filters, but more of an issue for steeper X-ray ones
        tr[:,i]=10**np.interp(logt,tresp_logt,np.log10(truse[:,i]))
#     Previous version
#         tr[:,i]=np.interp(logt,tresp_logt,truse[:,i])

    rmatrix=np.zeros([nt,nf])
    #Put in the 1/K factor (remember doing it in logT not T hence the extra terms)
    dlogTfac=10.0**logt*np.log(10.0**dlogt)
    # Do regularization of EMD or DEM 
    if emd_int:
        l_emd=True
        for i in np.arange(nf):
            rmatrix[:,i]=tr[:,i]
    else:
        for i in np.arange(nf):
            rmatrix[:,i]=tr[:,i]*dlogTfac
    #Just scale so not dealing with tiny numbers
    sclf=1E15
    rmatrix=rmatrix*sclf

    dn1d=np.reshape(dn,[nx*ny,nf])
    edn1d=np.reshape(edn,[nx*ny,nf])
#create our 1d arrays for output
    dem1d=np.zeros([nx*ny,nt])
    chisq1d=np.zeros([nx*ny])
    edem1d=np.zeros([nx*ny,nt])
    elogt1d=np.zeros([nx*ny,nt])
    dn_reg1d=np.zeros([nx*ny,nf])


# *****************************************************
#  Actually doing the DEM calculations
# *****************************************************
# Should always be just running the first part of if here as setting dem01d to array of 1s if nothing given
# So now more a check dimensions of things are correct 
    if ( dem0.ndim==dn.ndim ):
        dem01d=np.reshape(dem0,[nx*ny,nt])
        dem1d,edem1d,elogt1d,chisq1d,dn_reg1d=demmap_pos(dn1d,edn1d,rmatrix,logt,dlogt,glc,\
            reg_tweak=reg_tweak,max_iter=max_iter,rgt_fact=rgt_fact,dem_norm0=dem01d,nmu=nmu,warn=warn,l_emd=l_emd,rscl=rscl)
    else:
        dem1d,edem1d,elogt1d,chisq1d,dn_reg1d=demmap_pos(dn1d,edn1d,rmatrix,logt,\
            dlogt,glc,reg_tweak=reg_tweak,max_iter=max_iter,\
                rgt_fact=rgt_fact,dem_norm0=0,nmu=nmu,warn=warn,l_emd=l_emd,rscl=rscl)
    #reshape the 1d arrays to original dimensions and squeeze extra dimensions
    dem=((np.reshape(dem1d,[nx,ny,nt]))*sclf).squeeze()
    edem=((np.reshape(edem1d,[nx,ny,nt]))*sclf).squeeze()
    elogt=(np.reshape(elogt1d,[ny,nx,nt])/(2.0*np.sqrt(2.*np.log(2.)))).squeeze()
    chisq=(np.reshape(chisq1d,[nx,ny])).squeeze()
    dn_reg=(np.reshape(dn_reg1d,[nx,ny,nf])).squeeze()

    # There's probably a neater way of doing this (and maybe provide info of what was done as well?)
    # but fine for now as it works
    if emd_int and emd_ret:
        return dem,edem,elogt,chisq,dn_reg
    if emd_int and not emd_ret:
        return dem/dlogTfac,edem/dlogTfac,elogt,chisq,dn_reg
    if not emd_int and emd_ret:
        return dem*dlogTfac,edem*dlogTfac,elogt,chisq,dn_reg
    if not emd_int and not emd_ret:
        return dem,edem,elogt,chisq,dn_reg