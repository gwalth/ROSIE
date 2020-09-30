#!/usr/bin/env python

import argparse,os,sys


import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# usage:
#   cd /Users/gwalth/data/IMACS/REDUCTIONS/ROSIE/RIFUSLIT1_bkup
#   ~/python.linux/dev/ROSIE/ROSIE_cube.py  ift2005c_rsum.fits

def extract():

    return None

    

parser = argparse.ArgumentParser(description='Construction of 3D cube')

parser.add_argument('fits_file', metavar='file', type=str, nargs='?',
                    default=None,help='FITS file')
#parser.add_argument('-dir', metavar='file', type=str, nargs='?',
#                    default=None,help='directory of 2d spectra')

#args = parser.parse_args(namespace=)
args = parser.parse_args()
#print args
#print dir(args)
#print args.__dict__.keys()
fits_file = args.fits_file

# variables
ext = 0 
plot = 0


if not fits_file:
    print "Need at least one argument!"
    sys.exit()



A = []

i = 0
xsize_max = -1
ysize_max = -1


pf = fits.open(fits_file)
#print len(pf)
data = pf[ext].data
head = pf[ext].header

crpix1  = head['crpix1']     # starting pixel
crval1  = head['crval1']     # starting wavelength
cdelt1  = head['cdelt1']      # dispersion
dc_flag = head['dc-flag']    # Log-linear flag

nslits  = head['nslits']


objid  = [head["objid%i" % (i+1)] for i in range(nslits)]
csecta = [head["csect%ia" % (i+1)] for i in range(nslits)]
csectb = [head["csect%ib" % (i+1)] for i in range(nslits)]

ysize,xsize = data.shape

ysize_csect = np.max(np.array(csectb)-np.array(csecta))
N = len(objid)

print data.shape
A = np.zeros((xsize,N,ysize_csect))
i = 0
for obj,a,b in zip(objid,csecta,csectb):
    print obj,a,b
    slitlet = data[a-1:b-1,:]
    print slitlet.shape

    #A.append(slitlet)
    ysize,xsize = slitlet.shape
    A[:xsize,i,:ysize] = np.swapaxes(slitlet,0,1)[::-1,::-1]

    i += 1
    

    if plot:

        vmin = np.percentile(data,5)
        vmax = np.percentile(data,95)

        fig = plt.figure()
        p = fig.add_subplot(111)
        #p.imshow(data)
        p.imshow(slitlet, vmin=vmin, vmax=vmax, origin='lower',
                 interpolation="nearest", cmap=cm.gray_r,
                 extent=(w[0],w[-1],0,xsize),
                )
        #p.hist(data.flatten(),bins=100)
        #p.set_yscale("log")
        plt.show()


# construct header

new_head = head.copy()

# WCS
# need correction for pointing etc
# CD matrix from rotation angle etc

# new_wcs = wcs.WCS(naxis=2)
# new_wcs.wcs.crpix = [head["crpix1"],head["crpix2"]]
# new_wcs.wcs.cdelt = np.array([head["cdelt1"],head["cdelt2"]])
# new_wcs.wcs.crval = [head["crval1"],head["crval2"]]
# new_wcs.wcs.ctype = [head["ctype1"],head["ctype2"]]

new_crval1 = head['ra-d']
new_crval2 = head['dec-d']
new_crpix1 = ysize_max/2.
new_crpix2 = N/2.
new_cdelt1 = 0.2/3600.# IMACS f/2 pixel scale
#            0.207 mm? from SMF file
new_cdelt2 = 0.207 * 2.89856 / 3600.  # what is the slice width?

# rescale pixel shifts into mm in the focal plane
# for f/2 the scale is (1/2.89856) mm per arcsec
#  MMperPixel = 0.2 * 0.34500
#  PixelperMM = 1/(0.2 * 0.34500)

#new_head["naxis"] = 3
#new_head["naxis1"] = N
#new_head["naxis2"] = ysize_max
#new_head["naxis3"] = xsize_max

new_head["wcsaxes"] = 3
new_head["wcsname"] = 'ROSIEWCS'
new_head["radesys"] = 'FK5'


# all the slices
new_head["crval1"] = new_crval1
new_head["crpix1"] = new_crpix1
new_head["cdelt1"] = new_cdelt1
new_head["cd1_1"]  = new_cdelt1
new_head["ctype1"] = "RA---TAN"
new_head["cunit1"] = "deg"

# each slice
new_head["crval2"] = new_crval2
new_head["crpix2"] = new_crpix2
new_head["cdelt2"] = new_cdelt2
new_head["cd2_2"]  = new_cdelt2
new_head["ctype2"] = "DEC--TAN"
new_head["cunit2"] = "deg"

# dispersion direction
new_head["crval3"] = crval1
new_head["crpix3"] = crpix1
new_head["cdelt3"] = cdelt1
new_head["cd3_3"]  = cdelt1
new_head["ctype3"] = "wave"
new_head["cunit3"] = "Angstroms"
new_head['dc-flag'] = dc_flag    # Log-linear flag

   
output = fits_file.replace(".fits","") + "_cube.fits"
hdu = fits.PrimaryHDU(A, header=new_head)
hdu.writeto(output, clobber=True)
