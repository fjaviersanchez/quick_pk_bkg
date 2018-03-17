import numpy as np
import os
import fitsio
import skimage.transform
import matplotlib.pyplot as plt
import astropy.table
from scipy import fftpack
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--run-number',dest='rn',type=int,help='Run number to analyze')
parser.add_option('--show-plots',dest='show_plot',action='store_true',help='Show results')
(o, args) = parser.parse_args()
run_no = o.rn
test_number='%06d' %run_no # Number of the run to analyze
rfactor=16 # Rebinning factor
# I get the path to the e-images
data_path = ['/global/projecta/projectdirs/lsst/production/DC2/DC2-phoSim-test_WFD-r/output/%s/lsst_e_417009_f2_R22_S%s_E000.fits.gz' % (test_number,sensor) for sensor in ['00','01','02','10','11','12','20','21','22']]
# I read the e-images (using fitsio to improve speed)
data = [fitsio.read(dp) for dp in data_path]
images = [skimage.transform.resize(d,(d.shape[0]/rfactor,d.shape[1]/rfactor),preserve_range=True) for d in data]
# Check the background in one of the chips
if o.show_plot:
    plt.hist(images[0].flatten(),range=(400,1000),bins=100);
    plt.show()
# Merge the 9 images into 1 (I am not using the edges properly)
total_data = np.zeros((images[0].shape[0]*3,images[0].shape[1]*3))
for i in range(0,3):
    for j in range(0,3):
        total_data[images[0].shape[0]*i:images[0].shape[0]*(i+1),images[0].shape[1]*j:images[0].shape[1]*(j+1)]=images[3*i+j]
# FFT of the density contrast
F1 = fftpack.fft2((total_data/np.mean(total_data)-1))
 
# Now shift the quadrants around so that low spatial frequencies are in
# the center of the 2D fourier transformed image.
F2 = fftpack.fftshift( F1 )
 
# Calculate a 2D power spectrum
psd2D = np.abs( F2 )**2
if o.show_plot:
    plt.figure()
    plt.imshow( np.log10( total_data ), cmap='gray_r') 
    plt.show()
    plt.figure()
    plt.imshow( np.log10( psd2D )) 
    plt.show()

pix_scale=0.2/60*rfactor #pixel scale in arcmin
kx = 1./pix_scale*np.arange(-F2.shape[0]/2,F2.shape[0]/2)/(2.*np.pi)*1./F2.shape[0]
ky = 1./pix_scale*np.arange(-F2.shape[1]/2,F2.shape[1]/2)/(2.*np.pi)*1./F2.shape[1]
kxx, kyy = np.meshgrid(kx,ky)
rad = np.sqrt(kxx**2+kyy**2)
bins = 1./pix_scale*np.arange(0,F2.shape[0]/2)/(2.*np.pi)*1./F2.shape[0]
bin_space = bins[1]-bins[0]
ps1d = np.zeros(len(bins))
for i,b in enumerate(bins):
    ps1d[i]=np.mean(psd2D.T[(rad>b-0.5*bin_space) & (rad<b+0.5*bin_space)])/(F2.shape[0]*F2.shape[1])
if o.show_plot:
    plt.loglog(bins,ps1d)
    plt.ylim(1,3e2)
    plt.xlabel('$k$ [1/arcmin]')
    plt.ylabel('$P(k)$')
    plt.show()
tab_pk = astropy.table.Table([bins,ps1d],names=('k','Pk'))
tab_pk.write('Pk_%d.fits.gz' %run_no,overwrite=True)
