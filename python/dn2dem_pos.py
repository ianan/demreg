import numpy as np
import scipy.interpolate
import astropy.units as u
from astropy.units import imperial
from astropy import time
from demmap_pos import demmap_pos
imperial.enable()


def dn2dem_pos(dn_in,edn_in,tresp,tresp_logt,temps,reg_tweak=1.0,max_iter=10,gloci=0,rgt_fact=1.5,dem_norm0=None):
    # Performs a Regularization on solar data, returning the Differential Emission Measure (DEM)
    # using the method of Hannah & Kontar A&A 553 2013
    # Basically getting DEM(T) out of g(f)=K(f,T)#DEM(T)


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
        truse[tresp[:,i] <= 0]=np.min(tresp[tresp[:,i] > 0])

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
