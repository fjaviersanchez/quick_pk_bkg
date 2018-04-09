# quick_pk_bkg
Quick checks of background in LSST simulated images

Requirements:

`astropy.table`
`fitsio`
`skimage`
`scipy`
`numpy`
`pandas`

`quick_pk.py` computes the background's power-spectrum for a given raft. The usage is the following:

`python quick_pk.py --input-path $PATH_TO_LSST_IMAGES --output-path OUTPUT_NAME --raft-number RAFT_NUMBER`


`quick_bkg_check.py` computes the background of all e-images in the `--input-path` directory and compares either to the imSim's `ESOSkyModel` (if imSim is available) or to OpSim's value. It produces histograms of the background's mean, median andstandard deviation, adding a line for OpSim's/imSim's predicted values. It prints additional information if it detects that some images have very different background than predicted (in at least 50 percent) or if they differ a lot from the median (in at least 20 percent). Usage:
 
`python quick_bkg_check.py --input-path $PATH_TO_EIMAGES`
