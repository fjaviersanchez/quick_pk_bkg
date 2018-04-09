from __future__ import print_function, division
import numpy as np
import os
import fitsio
import matplotlib.pyplot as plt
import pandas as pd
import glob
from optparse import OptionParser
try:
    from lsst.sims.GalSimInterface import LSSTCameraWrapper
    from desc.imsim.skyModel import ESOSkyModel
    import desc.imsim
    from lsst.sims.photUtils import LSSTdefaults, PhotometricParameters
    from lsst.sims.utils import ObservationMetaData, radiansFromArcsec
    from lsst.sims.GalSimInterface import make_galsim_detector
    from lsst.sims.photUtils import BandpassDict
    from lsst.sims.GalSimInterface import GalSimInterpreter
    imsim_installed=True
except:
    imsim_installed=False
    print('imSim not installed, using OpSim DB to check background levels')

import sqlite3
from sqlite3 import Error
 
def create_connection(db_file):
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)
 
    return None
# I connect to the minion database to get additional information about the pointing

conn = create_connection('/global/cscratch1/sd/descpho/Pipeline-tasks/DBstaging/minion_1016_desc_dithered_v4.db')
'/global/cscratch1/sd/descpho/Pipeline-tasks/DBstaging/minion_1016_desc_dithered_v4.db'
bandpass_all = ['u','g','r','i','z','y']

def get_airmass_raw_seeing(conn,visit):
    cur = conn.cursor()
    cur.execute("SELECT airmass, filtSkyBrightness, finSeeing, rawSeeing, visitExpTime, fiveSigmaDepth FROM ObsHistory WHERE obsHistID==%d" %(visit))
    rows = cur.fetchall()
    return rows[0][0], rows[0][1], rows[0][2], rows[0][3], rows[0][4], rows[0][5]

def compute_sky_counts(mag,band,nsnap):
    # Data from https://github.com/lsst-pst/syseng_throughputs/blob/master/plots/table2
    if band == 'u':
        mag0 = 22.95
        counts0 = 50.2
    if band == 'g':
        mag0 = 22.24
        counts0 = 384.6
    if band == 'r':
        mag0 = 21.20
        counts0 = 796.2
    if band == 'i':
        mag0 = 20.47
        counts0 = 1108.1
    if band == 'z':
        mag0 = 19.60
        counts0 = 1687.9
    if band == 'y':
        mag0 = 18.63
        counts0 = 2140.8 
    return nsnap*counts0*10**(-0.4*(mag-mag0))

def compute_bkg(image):
    image = image.flatten()
    q95 = np.percentile(image,95)
    q5 = np.percentile(image,5)
    mask = (image>q5) & (image<q95)
    median_bkg = np.median(image[mask])
    mean_bkg = np.mean(image[mask])
    bkg_noise = np.std(image[mask])
    return mean_bkg, median_bkg, bkg_noise

def compare_to_imsim(phosim_commands):
    bandpass = bandpass_all[int(phosim_commands['filter'].values[0])]
    obs_md = ObservationMetaData(pointingRA=float(phosim_commands['rightascension'].values[0]),
                                pointingDec=float(phosim_commands['declination'].values[0]),
                                 mjd=float(phosim_commands['mjd'].values[0]),
                                 rotSkyPos=float(phosim_commands['rotskypos'].values[0]),
                                 bandpassName=bandpass,
                                 m5=LSSTdefaults().m5(bandpass),
                                 seeing=float(phosim_commands['seeing'].values[0]))
    noise_and_background = ESOSkyModel(obs_md, addNoise=True, addBackground=True)
    phot_params = PhotometricParameters(exptime=float(phosim_commands['vistime'].values[0]),
                                 nexp=int(phosim_commands['nsnap'].values[0]),
                                 gain=1,
                                 readnoise=0,
                                 darkcurrent=0,
                                 bandpass=bandpass)
    # We are going to check one sensor only
    detector_list = [make_galsim_detector(camera_wrapper, "R:2,2 S:1,1",
                                              phot_params, obs_md)]
    bp_dict = BandpassDict.loadTotalBandpassesFromFiles(bandpassNames=obs_md.bandpass)
    gs_interpreter = GalSimInterpreter(obs_metadata=obs_md,
                                       epoch=2000.0,
                                       detectors=detector_list,
                                       bandpassDict=bp_dict,
                                       noiseWrapper=noise_and_background,
                                       seed=1234)
    image = gs_interpreter.blankImage(detector=detector_list[0])
    image_2 = noise_and_background.addNoiseAndBackground(image,bandpass=obs_md.bandpass,
                                           m5=obs_md.m5,
                                           FWHMeff=obs_md.seeing,
                                       photParams=phot_params,
                                                     chipName=detector_list[0].name
                                      )
    return compute_bkg(image_2.array)

def get_total_seeing(h):
    total_seeing = 0.
    for key in list(h.keys()):
        if ('SEE' in key) & (key !='SEED'):
            total_seeing = total_seeing + h[key]**2
    return np.sqrt(total_seeing)

parser = OptionParser()

parser.add_option('--input-path',dest='input_path',type=str,help='Input directory')
parser.add_option('--output-path',dest='output_path',type=str,help='Output path')

(o, args) = parser.parse_args()

data_path = glob.glob(os.path.join(o.input_path,'lsst_e_*.fits.gz'))
parent_path = os.path.dirname(o.input_path)
instcat_path = glob.glob(os.path.join(parent_path,'instCat/phosim*.txt'))[0]
phosim_pars = pd.read_table(instcat_path, index_col=0, header=None, sep=' ').T
seeing = float(phosim_pars['seeing'].values[0])
bandpass = bandpass_all[int(phosim_pars['filter'].values[0])]
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

if imsim_installed:
   mean_bkg_imsim, median_bkg_imsim, bkg_noise_imsim = compare_with_imsim(phosim_pars)
else:
    sky_brightness = get_airmass_raw_seeing(conn,int(phosim_pars['obshistid'].values[0]))[1]
    print('OpSim skybrightness is: %f in %s band' % (sky_brightness,bandpass))
    mean_bkg_imsim = compute_sky_counts(sky_brightness,bandpass,int(phosim_pars['nsnap'].values[0]))
    print('The expected sky-count level is: %f' % mean_bkg_imsim)
    median_bkg_imsim = mean_bkg_imsim
    bkg_noise_imsim = np.sqrt(mean_bkg_imsim)
plt.figure()
plt.hist(psf_fwhm, label='PhoSim')
plt.plot(seeing*np.ones(3),np.linspace(0,len(data_path),3),label='OpSim')
plt.legend(loc='best')
plt.xlabel('FWHM [arcsec]')
plt.figure()
plt.hist(mean_bkg, label='Mean', alpha=0.2)
plt.plot(mean_bkg_imsim*np.ones(3),np.linspace(0,len(data_path),3),label='OpSim Mean')
plt.hist(median_bkg, label='Median', alpha=0.2)
plt.plot(median_bkg_imsim*np.ones(3),np.linspace(0,len(data_path),3),label='OpSim Median')
plt.legend(loc='best')
plt.xlabel('Background counts [ADU]')
plt.figure()
plt.hist(bkg_noise, label='PhoSim')
plt.plot(bkg_noise_imsim*np.ones(3),np.linspace(0,len(data_path),3),label='OpSim')
plt.xlabel('Background noise [ADU]')

# Additional checks
bkg_far = np.where(np.fabs(median_bkg/median_bkg_imsim-1.)>0.5)[0]
if len(bkg_far)>0:
    print(bkg_far)
    print('Images with background significantly different (50 percent) than expected: ')
    print(np.array(data_path)[bkg_far])
median_median_bkg = np.median(median_bkg)
diff = median_median_bkg/median_bkg_imsim-1
if np.fabs(diff)>0.2:
    print('Median background in the focal plane %.1f percent different than expected' %(100*np.fabs(diff)))
    print('Check observing conditions, please: ')
    print('Moon altitude (deg): ', phosim_pars['moonalt'].values)
    print('Distance to moon (deg): ', phosim_pars['dist2moon'].values)
    print('Sun altitude (deg): ', phosim_pars['sunalt'].values)
plt.show()
