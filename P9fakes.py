import numpy as np
import os, random, glob, shutil, time, subprocess, math
from astropy.io import fits
from astropy import wcs
import warnings
#warnings.filterwarnings('ignore', category=UserWarning, append=True) #Ignores UserWarnings otherwise Astropy spits out loads when it overwrites files
#warnings.filterwarnings('ignore', category=Warning, append=True)
from astropy.io.fits import getheader
import sip_to_pv as s2p
import sys
__author__ = "Chris Frohmaier"
__email__ = "c.frohmaier2soton.ac.uk"
__status__ = "In development"

def file_structure(): #This definition creates the file structure for results to be added into
	if not os.path.exists('Output_Images_V'+str(vnum)+''):
		os.makedirs('Output_Images_V'+str(vnum)+'')
	if not os.path.exists('Results_V'+str(vnum)+''):
		os.makedirs('Results_V'+str(vnum)+'')
	if not os.path.exists('Results_V'+str(vnum)+'/Catalog'):
		os.makedirs('Results_V'+str(vnum)+'/Catalog')
	if not os.path.exists('Results_V'+str(vnum)+'/Fake_Star_Catalog'):
		os.makedirs('Results_V'+str(vnum)+'/Fake_Star_Catalog')
	if not os.path.exists('Results_V'+str(vnum)+'/Fakes_added'):
		os.makedirs('Results_V'+str(vnum)+'/Fakes_added')
	if not os.path.exists('Results_V'+str(vnum)+'/Galaxies'):
		os.makedirs('Results_V'+str(vnum)+'/Galaxies')

def Sextract(in_image, w_image,zeropoint,seeing,saturation,gain): #Runs Sextractor and creates a catalog of all the stars
	"""
	Run SExtractor

	Make sure you have the default.conv and default.nwm etc in this directory, or in the param file, specify the location of these files.
	"""

	PTFname=in_image.split('/')[-1].split(".fits")[0]
	subprocess.call('sex -c Pipe_sexfile_Peter.sex '+in_image+' -PARAMETERS_NAME PTF_Transform_Param.param -FILTER_NAME default.conv -CATALOG_NAME Results_V'+str(vnum)+'/Catalog/'+PTFname+'_Catalog_V'+str(vnum)+'.cat -WEIGHT_IMAGE '+w_image+' -MAG_ZEROPOINT'+'	'+str(zeropoint)+' -SEEING_FWHM '+str(seeing)+' -SATUR_LEVEL '+str(saturation)+' -GAIN '+str(gain)+' -PHOT_FLUXFRAC 0.2,0.5,0.9 -VERBOSE_TYPE QUIET',shell=True)

def Enough_Objects(in_image): #Checks that sextractor has found a suitable number of objects
	"""Counts the number of objects in the source extractor catalog

	"""

	PTFname=in_image.split('/')[-1].split(".fits")[0]
	enough=True
	test=os.popen('wc -l Results_V'+str(vnum)+'/Catalog/'+PTFname+'_Catalog_V'+str(vnum)+'.cat').read()
	rows=test.split()
	if float(rows[0])<300: 
			return False

def SelectBright(in_image): #Selected the top 20 brightest stars in the catalog
	"""Will randomly select a bright star in the catalog to become our fake
	catP9in columns
	#   1 NUMBER          Running object number
	#   2 FLUX_AUTO       Flux within a Kron-like elliptical aperture     [count]
	#   3 MAG_AUTO        Kron-like elliptical aperture magnitude         [mag]
	#   4 X_IMAGE         Object position along x                         [pixel]
	#   5 Y_IMAGE         Object position along y                         [pixel]
	#   6 BACKGROUND      Background at centroid position                 [count]
	#   7 ELONGATION      A_IMAGE/B_IMAGE
	#   8 ELLIPTICITY     1 - B_IMAGE/A_IMAGE
	#   9 FLAGS           Extraction flags
	#  10 CLASS_STAR      S/G classifier output
	#  11 CXX_IMAGE       Cxx object ellipse parameter                    [pixel**(-2)]
	#  12 CYY_IMAGE       Cyy object ellipse parameter                    [pixel**(-2)]
	#  13 CXY_IMAGE       Cxy object ellipse parameter                    [pixel**(-2)]
	#  14 X2_IMAGE        Variance along x                                [pixel**2]
	#  15 Y2_IMAGE        Variance along y                                [pixel**2]
	#  16 XY_IMAGE        Covariance between x and y                      [pixel**2]
	#  17 FWHM_IMAGE      FWHM assuming a gaussian core                   [pixel]
	#  18 MAG_BEST        Best of MAG_AUTO and MAG_ISOCOR                 [mag]
	#  19 ALPHA_SKY       Right ascension of barycenter (native)          [deg]
	#  20 DELTA_SKY       Declination of barycenter (native)              [deg]
	#  21 FLUX_RADIUS     Fraction-of-light radii                         [pixel]
	"""
	PTFname=in_image.split('/')[-1].split(".fits")[0]
	catP9in=np.loadtxt('Results_V'+str(vnum)+'/Catalog/'+PTFname+'_Catalog_V'+str(vnum)+'.cat',skiprows=21, comments="#")
	source=catP9in[(catP9in[:,2]>15) & (catP9in[:,7]<0.3) & (catP9in[:,9]>0.5) & (catP9in[:,8]==0) & (catP9in[:,3]>100.) & (catP9in[:,3]<1948.) & (catP9in[:,4]>100.) & (catP9in[:,4]<3996.)]
	source=source[source[:,2].argsort()][:20]
	source=source[np.random.randint(len(source), size=1),:][0]
	return source

def scaleStar(source, fmag, zpt): #def Scaling(science_image ,xcord, ycord, mag_array, flux_array, background_array, zpt, fake_stars, CCD_Num, magnitude_best,alpha_sky, delta_sky):
	fakeFlux=10.0**((fmag-zpt)/(-2.5))
	scaling_factor=((fakeFlux)/source[1])
	print scaling_factor
	return scaling_factor

def addFake(inimage,sourcearr, hostlessx, hostlessy, scale_fac):
	science_data=inimage[0].data
	#---Old area to be scaled---
	sourcex=sourcearr[3]
	sourcey=sourcearr[4]
	back=sourcearr[5]
	startx=int(sourcex-10.0)
	starty=int(sourcey-10.0)
	finx=int(sourcex+10.0)
	finy=int(sourcey+10.0)

	#---New area to have flux added---
	Nstartx=hostlessx-10.0
	Nstarty=hostlessy-10.0
	Nfinx=hostlessx+10.0
	Nfiny=hostlessy+10.0

	fbox1=np.sum(science_data[hostlessy,hostlessx])
	fbox2=np.sum(science_data[hostlessy,hostlessx]) + np.sum(science_data[hostlessy-1.0,hostlessx]) + np.sum(science_data[hostlessy+1.0,hostlessx]) + np.sum(science_data[hostlessy,hostlessx-1.0]) + np.sum(science_data[hostlessy,hostlessx+1.0])
	fbox3=np.sum(science_data[hostlessy-1.0:hostlessy+2.0, hostlessx-1.0:hostlessx+2.0])
	fbox4=np.sum(science_data[hostlessy-1.0:hostlessy+2.0, hostlessx-1.0:hostlessx+2.0]) + np.sum(science_data[hostlessy-2.0,hostlessx]) + np.sum(science_data[hostlessy+2.0,hostlessx]) + np.sum(science_data[hostlessy, hostlessx-2.0]) + np.sum(science_data[hostlessy, hostlessx+2.0])
	fbox5=np.sum(science_data[hostlessy-2.0:hostlessy+3.0, hostlessx-2.0:hostlessx+3.0])
	fbox6=np.sum(science_data[hostlessy-5.0:hostlessy+6.0, hostlessx-5.0:hostlessx+6.0])

	newdata=np.ones((20,20)) #Preparing a blank gird for scaled objects
	newdata[0:20,0:20]=(((science_data[starty:finy,startx:finx]))-back)*scale_fac #inserting scaled object
	science_data[Nstarty:Nfiny, Nstartx:Nfinx]= (science_data[Nstarty:Nfiny, Nstartx:Nfinx]) + newdata
	
	inimage.flush()
	inimage.close()
global vnum
vnum=int(sys.argv[2])
file_structure()
image_list=np.genfromtxt(sys.argv[1], dtype=None)

for i in range(len(image_list)):
	oimage=image_list[i][0]
	wimage=image_list[i][1]
	Rafake=image_list[i][2]
	Decfake=image_list[i][3]
	fakemag=image_list[i][4]

	fname=oimage.split('/')[-1].split(".fits")[0]+"_P9fakes_V"+str(vnum)+".fits"


	s2p.sip_to_pv(oimage,'Output_Images_V'+str(vnum)+'/'+fname, preserve=True)
	print 'SIP Start'
	hdulist_multi_sci=fits.open('Output_Images_V'+str(vnum)+'/'+fname, mode='update')
	print 'SIP Finish'

	w=wcs.WCS(hdulist_multi_sci[0].header)

	#print oimage
	print "---------"
	print 'Ra, Dec (in): ', image_list[i][2],image_list[i][3]
	x, y = w.all_world2pix(Rafake,Decfake,1)

	print 'X, Y (in: ', (x,y) 

	ra, dec= w.all_pix2world(x,y,0)

	print 'Ra, Dec (out): ', (ra, dec)
	print "---------"

	zeropoint=float(hdulist_multi_sci[0].header['UB1_ZP'])
	seeing=float(hdulist_multi_sci[0].header['SEEING'])
	saturation=55000.0 #float(hdulist_multi_sci[0].header['SATURATE'])
	gain=float(hdulist_multi_sci[0].header['GAIN'])	
	Sextract(oimage, wimage,zeropoint,seeing,saturation,gain)	
	catsize=Enough_Objects(oimage)
	sourcestar=SelectBright(oimage)
	sfactor=scaleStar(sourcestar, fakemag, zeropoint)
	addFake(hdulist_multi_sci,sourcestar, int(round(x)), int(round(y)), sfactor)
	print "Done"
	