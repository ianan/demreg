import numpy as np
from matplotlib import pyplot as plt
from demmap_pos import demmap_pos
from dn2dem_pos import dn2dem_pos
import scipy.interpolate
from aiapy.calibrate import degradation, register, update_pointing
from aiapy.calibrate.util import get_correction_table
import aiapy.response
from astropy import time
import astropy.units as u
from astropy.units import imperial
from astropy.visualization import time_support
from pandas import read_csv
import os
from sunpy.net import Fido, attrs
from sunpy.map import Map
from sunpy.instr.aia import aiaprep
import dateutil.parser
import cProfile
import pstats
from io import StringIO
import threadpoolctl
threadpoolctl.threadpool_limits(1)

nt=14
fits_dir='/mnt/c/Users/Alasdair/Documents/reginvpy/test/'
os.chdir('/mnt/c/Users/Alasdair/Documents/reginvpy')
# os.chdir('C:/Users/Alasdair/Documents/reginvpy')
# fits_dir="C:/Users/Alasdair/Documents/reginvpy/test/"
# fits_dir='/home/awilson/code/DEM/demreg-py/demreg-py/test/'
# os.chdir('/home/awilson/code/DEM/demreg-py/demreg-py/')
temperatures=10**np.linspace(5.7,7.1,num=nt+1)
tresp = read_csv('tresp.csv').to_numpy()
# print(tresp_logt.keys())
# data=np.ones([nx,ny,nf])
# edata=np.ones([nx,ny,nf])/10
# dem_norm=np.ones([nx,ny,nt])
data=np.array([3.4,13.8,184,338,219.55,12.22])
edata=np.array([0.2,0.43,7.83,12.9,5.80,0.23])
dem_norm=np.array([ 0.082588151,0.18005607,0.30832890,0.47582966, 0.66201794,0.83059740,0.93994260,0.95951378 ,0.88358527,0.73393929, 0.54981130, 0.37136465,0.22609001 , 0.11025056])


correction_table = get_correction_table()
wavenum=['94','131','171','193','211','335']
channels = []
for i in np.arange(len(wavenum)):
    channels.append(float(wavenum[i])*u.angstrom)

time_calibration = time.Time('2014-01-01T00:00:00', scale='utc')

time_test = time.Time('2014-01-01T00:00:00', scale='utc')

# deg_calibration = {}
deg_calibration = np.zeros([len(channels)])
# deg = {}
deg = np.zeros([len(channels)])

for i,c in enumerate(channels):
    deg_calibration[i] = degradation(c,time_calibration, correction_table=correction_table)
    deg[i] = degradation(c,time_test, correction_table=correction_table)

tresp_logt=tresp[:,0]
tresp_calibration=tresp[:,1:]/deg_calibration
trmatrix=deg[:]*tresp_calibration

# time_0 = astropy.time.Time('2010-06-01T00:00:00', scale='utc')
# now = astropy.time.Time.now()
# time = time_0 + np.arange(0, (now - time_0).to(u.day).value, 7) * u.day

# deg = {}
# for c in channels:
#     deg[c] = [degradation(c*u.angstrom, t, correction_table=correction_table) for t in time]
# time_support()  # Pass astropy.time.Time directly to matplotlib
# fig = plt.figure()
# ax = fig.gca()
# for i,c in enumerate(channels):
#     ax.plot(time, deg[c])
# ax.set_xlim(time[[0, -1]])
# ax.legend(frameon=False, ncol=4, bbox_to_anchor=(0.5, 1), loc='lower center')
# ax.set_xlabel('Time')
# ax.set_ylabel('Degradation')
# plt.show()

# dem,edem,elogt,chisq,dn_reg=dn2dem_pos_nb(data,edata,trmatrix,tresp_logt,temperatures,dem_norm0=dem_norm,max_iter=50)

# plt.plot(np.log10(temperatures[:-1]),np.log10(dem))
# plt.show()


# time_test = time.Time('2012-01-01T00:00:00', scale='utc')

# td = time.TimeDelta(11,format='sec')
# q=Fido.search(
#     attrs.Time(time_test, time_test+td),
#     attrs.Instrument('AIA'),
#     attrs.Wavelength(channels[0]) | attrs.Wavelength(channels[1]) | attrs.Wavelength(channels[2]) | attrs.Wavelength(channels[3]) | attrs.Wavelength(channels[4]) | attrs.Wavelength(channels[5]),
# )
# print(q)
# files=Fido.fetch(q)
# # maps=sunpy.map.Map(files)
# # maps.peek(vmin=0)
# # print(files)

# maps = [Map(f) for f in files]
# maps = sorted(maps, key=lambda x: x.wavelength)
# maps = [aiaprep(m) for m in maps]
# maps = [Map(m.data/m.exposure_time.value, m.meta) for m in maps]


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
aia = Map(fits_files)

#read dimensions from the header
nx=int(aia[0].dimensions.x.value)
ny=int(aia[0].dimensions.y.value)
nf=len(file_str)

#normalise to dn/s
aia = [Map(m.data/m.exposure_time.value, m.meta) for m in aia]
#create data array
data=np.zeros([nx,ny,nf])
#convert from our list to an array of data
for j in np.arange(nf):
    data[:,:,j]=aia[j].data
data[data < 0]=0

#calculate our dem_norm guess
off=0.412
gauss_stdev=12
dem_norm0=np.zeros([nx,ny,nt])
dem_norm_temp=np.convolve(np.exp(-(np.arange(nt)+1-(nt-2)*(off+0.1))**2/gauss_stdev),np.ones(3)/3)[1:-1]
dem_norm0[:,:,:]=dem_norm_temp


serr_per=10.0
#errors in dn/px/s
npix=4096.**2/(nx*ny)
edata=np.zeros([nx,ny,nf])
gains=np.array([18.3,17.6,17.7,18.3,18.3,17.6])
dn2ph=gains*[94,131,171,193,211,335]/3397.0
rdnse=1.15*np.sqrt(npix)/npix
drknse=0.17
qntnse=0.288819*np.sqrt(npix)/npix
for j in np.arange(nf):
    etemp=np.sqrt(rdnse**2.+drknse**2.+qntnse**2.+(dn2ph[j]*abs(data[:,:,j]))/(npix*dn2ph[j]**2))
    esys=serr_per*data[:,:,j]/100.
    edata[:,:,j]=np.sqrt(etemp**2. + esys**2.)


# x1=300
# x2=500
# y1=600
# y2=800
x1=300
x2=500
y1=600
y2=601
pr = cProfile.Profile()
pr.enable()
dem,edem,elogt,chisq,dn_reg=dn2dem_pos(data[x1:x2,y1:y2,:],edata[x1:x2,y1:y2,:],trmatrix,tresp_logt,temperatures,dem_norm0=dem_norm0[x1:x2,y1:y2,:],max_iter=20)

# dem,edem,elogt,chisq,dn_reg=dn2dem_pos(data,edata,trmatrix,tresp_logt,temperatures,dem_norm0=dem_norm0,max_iter=20)
pr.disable()
s = StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())