from __future__ import print_function, division
import numpy as np
import os
import fitsio
import matplotlib.pyplot as plt
import pandas as pd
import glob
from optparse import OptionParser

def compute_bkg(image):
    image = image.flatten()
    q95 = np.percentile(image,95)
    q5 = np.percentile(image,5)
    mask = (image>q5) & (image<q95)
    median_bkg = np.median(image[mask])
    mean_bkg = np.mean(image[mask])
    bkg_noise = np.std(image[mask])
    return mean_bkg, median_bkg, bkg_noise

def get_total_seeing(h):
    total_seeing = 0.
    for key in list(h.keys()):
        if ('SEE' in key) & (key !='SEED'):
            total_seeing = total_seeing + h[key]**2
    total_seeing = total_seeing**0.5
    return total_seeing

parser = OptionParser()
parser.add_option('--show-plots',dest='show_plot',action='store_true',help='Show results')
parser.add_option('--input-path',dest='input_path',type=str,help='Input directory')
parser.add_option('--output-path',dest='output_path',type=str,help='Output path')

(o, args) = parser.parse_args()

data_path = glob.glob(os.path.join(o.input_path,'lsst_e_*.fits.gz'))
parent_path = os.path.dirname(o.input_path)
instcat_path = glob.glob(os.path.join(parent_path,'instCat/phosim*.txt'))[0]
fh = pd.read_table(instcat_path, index_col=0, header=None, sep=' ').T
seeing = float(fh['seeing'].values[0])

psf_fwhm = []
mean_bkg = []
median_bkg = []
bkg_noise = []

for i,dp in enumerate(data_path):
    if i%10==0:
        print('Analyzed %d of %d images' % (i,len(data_path)))
    data, h = fitsio.read(dp,ext=0,header=True)
    psf_fwhm.append(get_total_seeing(h))
    aux1, aux2, aux3 = compute_bkg(data)
    mean_bkg.append(aux1)
    median_bkg.append(aux2)
    bkg_noise.append(aux3)

psf_fwhm = np.array(psf_fwhm)
mean_bkg = np.array(mean_bkg)
median_bkg = np.array(median_bkg)
bkg_noise = np.array(bkg_noise)

plt.figure()
plt.hist(psf_fwhm, label='PhoSim')
plt.plot(seeing*np.ones(3),np.linspace(0,len(data_path),3))
plt.xlabel('FWHM [arcsec]')
plt.figure()
plt.hist(mean_bkg, label='Mean', alpha=0.2)
plt.hist(median_bkg, label='Median', alpha=0.2)
plt.legend(loc='best')
plt.xlabel('Background counts [ADU]')
plt.figure()
plt.hist(bkg_noise)
plt.xlabel('Background noise [ADU]')
plt.show()
