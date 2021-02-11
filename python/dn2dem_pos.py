import numpy as np
import scipy.interpolate
import astropy.units as u
from astropy.units import imperial
from astropy import time
from demmap_pos import demmap_pos
imperial.enable()


def dn2dem_pos(dn_in,edn_in,tresp,tresp_logt,temps,reg_tweak=1.0,max_iter=10,gloci=0,rgt_fact=1.5,dem_norm0=None):
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
        the temperature response matrix size n_tresp by nf
    tresp_logt:
        the temperatures in log t which the temperature response matrix corresponds to. E.G if your tresp matrix 
        runs from 5.0 to 8.0 in steps of 0.05 then this is the input to tresp_logt
    temps:
        the temperatures at which to calculate a DEM, array of length nt.

    --------------------
    Optional Inputs:
    --------------------

    dem_norm0:
        Highly recommended optional input!! This is an array of length nt which contains an initial guess of the DEM 
        solution providing a weighting for the inversion process. The actual values of the normalisation do not matter,
        only their relative values. Without dem_norm0 set the initial guess is taken as an arrays of 1 (i.e. all temperatures
        weighted equally)
    reg_tweak:
        the initial normalised chisq target.
    max_iter:
        the maximum number of iterations to attempt, code iterates if negative DEM os reached. If max iter is reached before
        a suitable solution is found then the current solution is returned instead (which may contain negative values)
    gloci:
        array of maximum length nf containing the indexes of the filters of which to use a loci curve as the dem_normalisation.
        can be used in place of dem_norm0
    rgt_fact:
        the factor by which rgt_tweak increases each iteration. As the target chisq increases there is more flexibility allowed 
        on the DEM

    --------------------
    Outputs:
    --------------------

    dem:
        The DEM, has shape nx*ny*nt and units cm^-5 K^-1
    edem:
        vertical errors on the DEM, same units.
    elogt:
        Horizontal errors on temperature.
    chisq:
        The final chisq, shape nx*ny. Pixels which have undergone more iterations will in general have higher chisq.
    dn_reg:
        The simulated dn counts, shape nx*ny*nf. This is obtained by multiplying the DEM(T) by the filter response K(f,T) for each channel
        useful for comparing with the initial data.
 
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
        tr[:,i]=np.interp(logt,tresp_logt,truse[:,i])

    rmatrix=np.zeros([nt,nf])
    #Put in the 1/K factor (remember doing it in logT not T hence the extra terms)
    for i in np.arange(nf):
        rmatrix[:,i]=tr[:,i]*10.0**logt*np.log(10.0**dlogt)
    #Just scale so not dealing with tiny numbers
    sclf=1E15
    rmatrix=rmatrix*sclf
    #time it
    t_start = time.Time.now()


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
# Do we have an initial DEM guess/constraint to send to demmap_pos as well?
    if ( dem0.ndim==dn.ndim ):
        dem01d=np.reshape(dem0,[nx*ny,nt])
        dem1d,edem1d,elogt1d,chisq1d,dn_reg1d=demmap_pos(dn1d,edn1d,rmatrix,logt,dlogt,glc,reg_tweak=reg_tweak,max_iter=max_iter,\
                rgt_fact=rgt_fact,dem_norm0=dem01d)
    else:
        dem1d,edem1d,elogt1d,chisq1d,dn_reg1d=demmap_pos(dn1d,edn1d,rmatrix,logt,\
            dlogt,glc,reg_tweak=reg_tweak,max_iter=max_iter,\
                rgt_fact=rgt_fact,dem_norm0=0)
    #reshape the 1d arrays to original dimensions and squeeze extra dimensions
    dem=((np.reshape(dem1d,[nx,ny,nt]))*sclf).squeeze()
    edem=((np.reshape(edem1d,[nx,ny,nt]))*sclf).squeeze()
    elogt=(np.reshape(elogt1d,[ny,nx,nt])/(2.0*np.sqrt(2.*np.log(2.)))).squeeze()
    chisq=(np.reshape(chisq1d,[nx,ny])).squeeze()
    dn_reg=(np.reshape(dn_reg1d,[nx,ny,nf])).squeeze()
    #end the timing
    t_end = time.Time.now()
    print('total elapsed time =', time.Time(t_end-t_start,format='datetime'))
    return dem,edem,elogt,chisq,dn_reg
