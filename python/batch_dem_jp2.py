def batch_dem_jp2(t_start,cadence,nobs,fits_dir,jp2_dir,get_fits=0,serr_per=10,min_snr=2,fe_min=2,sat_lvl=1.5e4)
    from matplotlib import pyplot as plt
    plt.rcParams['figure.figsize'] = [10, 9]  # make plots larger
    from astropy.time import Time, TimeDelta
    import sunpy.map
    from sunpy.instr.aia import aiaprep
    from sunpy.net import Fido, attrs as a
    import numpy as np
    import pprint

    from astropy.coordinates import SkyCoord
    from astropy import units as u

    import warnings
    warnings.filterwarnings("ignore")
    import dateutil.parser

    #set directory for fits files
    fits_dir='C:/Users/Alasdair/Documents/reginvpy/test/'
    #start datetime
    t_start='2014-01-01 00:00:00.000'
    #we only want optically thin coronal wavelengths
    wavenum=['94','131','171','193','211','335']

    #convert our string into a datetime object
    t=(dateutil.parser.parse(t_start))

    #deconstruct the datetime object into a synoptica data filename
    file_str=[('AIA'+str(t.year).zfill(4)+str(t.month).zfill(2)+str(t.day).zfill(2)+'_'+str(t.hour).zfill(2)+str(t.minute).zfill(2)+'_'+"{}".format(wave.zfill(4))+'.fits') for j,wave in enumerate(wavenum)]
    #find the files in their directory
    fits_files=[fits_dir+file_str[j] for j in np.arange(len(file_str))]

    #load the fits with sunpy
    aia = sunpy.map.Map(fits_files)

    #read dimensions from the header
    nx=int(aia[0].dimensions.x.value)
    ny=int(aia[0].dimensions.y.value)
    nf=len(file_str)
    pprint.pprint(vars(aia[0]))#.__dict__)

    # dx=aia[0].cdelt1
    # dy=aia[0].cdelt2

    #create data array
    data=np.zeros([nx,ny,nf])
    #convert from our list to an array of data
    for j in np.arange(nf):
        data[:,:,j]=aia[j].data

    #calculate the hot component of aia 94
    a94_fe18=np.zeros([nx,ny])
    a94_fe18[:,:]=data[:,:,0]-data[:,:,4]/120.0-data[:,:,2]/450.0
    #threshold of fe_min for the hot component
    fe_min=2
    a94_fe18[np.where(a94_fe18 < fe_min)]=0
    # plt.imshow(np.log10(a94_fe18))
    # plt.show()