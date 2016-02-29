import numpy as np
import os, subprocess
from astropy.io import fits
from astropy import wcs
import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True) #Ignores UserWarnings otherwise Astropy spits out loads when it overwrites files
warnings.filterwarnings('ignore', category=Warning, append=True)
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
	
	

def Sextract(in_image, w_image,zeropoint,seeing,saturation,gain): #Runs Sextractor and creates a catalog of all the stars
	"""
	Run SExtractor

	Make sure you have the default.conv and default.nwm etc in this directory, or give path in the param file to specify the location of these files.
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
	return scaling_factor, fakeFlux

def addFake(inimage, oimage,sourcearr, hostlessx, hostlessy,fmag,fflux, scale_fac, head):
	fakesRecord=open('Results_V'+str(vnum)+'/Fake_Star_Catalog/'+sys.argv[1].split('.')[0]+'_fakesAdded_V'+str(vnum)+'.dat','a')
	PTFname=oimage.split('/')[-1].split(".fits")[0]
	science_data=inimage[0].data
	##---Old area to be scaled---
	sourcex=sourcearr[3]
	sourcey=sourcearr[4]
	back=sourcearr[5]
	startx=int(sourcex-10.0)
	starty=int(sourcey-10.0)
	finx=int(sourcex+10.0)
	finy=int(sourcey+10.0)
	##---New area to have flux added---
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
	"""Sourcex, Sourcey, a_source, dec_source, x_loc, y_loc, source_mag_auto, source_mag_best, flux_source, mag_fake, flux_fake, background, scaling_factor, PTFField, CCD, fbox1, fbox2, fbox3, fbox4, fbox5, fbox6, gain, readnoise, MOONILLF, MoonRA, MoonDec, AIRMASS, seeing, ELLIP, MEDSKY, SKYSIG), zeropoint, LMT_MG, MJD"""

	fakesRecord.write(str(oimage)+' '+str(sourcearr[3])+' '+str(sourcearr[4])+' '+str(sourcearr[18])+' '+str(sourcearr[19])+' '+str(float(hostlessx))+' '+str(float(hostlessy))+' '+str(sourcearr[2])+' '+str(sourcearr[17])+' '+str(sourcearr[1])+' '+str(fmag)+' '+str(fflux)+' '+str(sourcearr[5])+' '+str(scale_fac)+' '+str(int(head[4]))+' '+str(head[5])+' '+str(fbox1)+' '+str(fbox2)+' '+str(fbox3)+' '+str(fbox4)+' '+str(fbox5)+' '+str(fbox6)+' '+str(head[3])+' '+str(head[6])+' '+str(head[7])+' '+str(head[14])+' '+str(head[15])+' '+str(head[8])+' '+str(head[1])+' '+str(head[9])+' '+str(head[10])+' '+str(head[11])+' '+str(head[0])+' '+str(head[12])+' '+str(head[13])+'\n')
	



image_list=np.genfromtxt(sys.argv[1], dtype=None)


for i in range(len(image_list)):

	oimage=image_list[i][0]
	wimage=image_list[i][1]
	Rafake=image_list[i][2]
	Decfake=image_list[i][3]
	fakemag=image_list[i][4]
	global vnum
	vnum=int(image_list[i][5])
	file_structure()
	fname=oimage.split('/')[-1].split(".fits")[0]+"_P9fakes_V"+str(vnum)+".fits"
	PTFname=oimage.split('/')[-1].split(".fits")[0]
	s2p.sip_to_pv(oimage,'Output_Images_V'+str(vnum)+'/'+fname, preserve=True)
	hdulist_multi_sci=fits.open('Output_Images_V'+str(vnum)+'/'+fname, mode='update')
	w=wcs.WCS(hdulist_multi_sci[0].header)
	x, y = w.all_world2pix(Rafake,Decfake,1)
	ra, dec= w.all_pix2world(x,y,0)
	"""Header Values"""
	zeropoint=float(hdulist_multi_sci[0].header['UB1_ZP'])
	seeing=float(hdulist_multi_sci[0].header['SEEING'])
	saturation=55000.0 #float(hdulist_multi_sci[0].header['SATURATE'])
	gain=float(hdulist_multi_sci[0].header['GAIN'])	
	CCD_Num=float(hdulist_multi_sci[0].header['CCDID'])
	PTFFIELD=int(hdulist_multi_sci[0].header['PTFFIELD'])
	readnoise=float(hdulist_multi_sci[0].header['READNOI'])
	MOONILLF=float(hdulist_multi_sci[0].header['MOONILLF'])
	AIRMASS=float(hdulist_multi_sci[0].header['AIRMASS'])
	ELLIP=(hdulist_multi_sci[0].header['ELLIP'])
	if ELLIP=='NAN.0':
		bad_images=open('Results_V'+str(vnum)+'/Bad_Images_V'+str(vnum)+'.dat','a')
		bad_images.write(str(science_image[0])+str(science_image[1])+str('.fits')+' '+str('Reason: ELLIP has a NAN')+'\n')
		bad_images.close()
		#print science_image[0]+science_image[1], ' Has a NAN'
		continue
	else:
		ELLIP=float(hdulist_multi_sci[0].header['ELLIP'])
	MEDSKY=float(hdulist_multi_sci[0].header['MEDSKY'])
	SKYSIG=float(hdulist_multi_sci[0].header['SKYSIG'])
	LMT_MG=float(hdulist_multi_sci[0].header['LMT_MG'])
	MJD=float(hdulist_multi_sci[0].header['OBSMJD'])
	MoonRA=float(hdulist_multi_sci[0].header['MOONRA'])
	MoonDec=float(hdulist_multi_sci[0].header['MOONDEC'])

	head=[zeropoint, seeing, saturation, gain, CCD_Num, PTFFIELD, readnoise, MOONILLF,AIRMASS, ELLIP, MEDSKY, SKYSIG, LMT_MG, MJD, MoonRA, MoonDec]

	Sextract(oimage, wimage,zeropoint,seeing,saturation,gain)	
	catsize=Enough_Objects(oimage)
	if catsize==False:
		bad_images=open('Results_V'+str(vnum)+'/Bad_Images_V'+str(vnum)+'.dat','a')
		bad_images.write(str(science_image[0])+str(science_image[1])+str('.fits')+' '+str('Reason: Sextractor did not detect enough objects (<300)')+'\n')
		os.remove('Results_V'+str(vnum)+'/Catalog/'+science_image[1]+'_Catalog_V'+str(vnum)+'.cat')
		continue
	sourcestar=SelectBright(oimage)
	sfactor, fflux=scaleStar(sourcestar, fakemag, zeropoint)
	addFake(hdulist_multi_sci, oimage ,sourcestar, int(round(x)), int(round(y)),fakemag,fflux, sfactor, head)
	
