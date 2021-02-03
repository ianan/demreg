from dn2dem_pos import dn2dem_pos
import numpy as np
def dn2dem_pos_selfnorm(data,edata,tresp,tresp_logt,temps,reg_tweak=1.0,max_iter=10,gloci=0,rgt_fact=1.5,dem_norm0=None):
    """
    Function to run dn2dem twice, using the output from the first attempt as a initial guess for a second attempt.
    You can still provide a dem_norm0 as an initial guess for the first attempt, if you do not dn2dem_pos will assume
    an array of ones.
    """
    #first call dn2dem as a first estimate of the dem
    dem,edem,elogt,chisq,dn_reg=dn2dem_pos(data,edata,tresp,tresp_logt,temps,dem_norm0=dem_norm0,max_iter=max_iter)
    dem_norm0=np.zeros(dem.shape)
    if len(data.shape)==1:
        #normalise the dem to the max value for that pixel and convolve with  a 5-wide window to smooth
        dem_norm0[:]=(np.convolve(dem[1:-1],np.ones(5)/5))[1:-1]/np.max(dem[:])
    if len(data.shape)==2:
        #for 1d arrays of data
        for ii in np.arange(dem.shape[0]):
           dem_norm0[ii,:]=(np.convolve(dem[ii,1:-1],np.ones(5)/5))[1:-1]/np.max(data[ii,:]) 
    if len(data.shape)==3:
        #for 2d maps of data
        for ii in np.arange(dem.shape[0]):
            for jj in np.arange(dem.shape[1]):
                dem_norm0[ii,jj,:]=(np.convolve(dem[ii,jj,1:-1],np.ones(5)/5))[1:-1]/np.max(dem[ii,jj,:])
    dem_norm0[dem_norm0<=1e-8]=1e-8
    #call dn2dem again with the new norm and return the dem.
    dem,edem,elogt,chisq,dn_reg=dn2dem_pos(data,edata,tresp,tresp_logt,temps,dem_norm0=dem_norm0,max_iter=max_iter)
    return dem,edem,elogt,chisq,dn_reg