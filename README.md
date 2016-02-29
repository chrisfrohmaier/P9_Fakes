# P9_Fakes
Generating fakes for the P9 hunt.

Based on the old supernova code for the PTF efficiencies.
##Run
```python
python P9fakes.py inputlist 
```
##Requirements

- Take an image
- Specify an Ra, Dec and magnitude
- Lay a fake down at that location and magnitude

##How to use the code

Run the python script with a ascii file as an arguement.

The ascii file structure for the input must be:

`path/to/image path/to/weight Ra_fake Dec_fake Mag_fake version_number`

The version number exists if you want to use the same new image multiple times

**Note:** Ra and Dec in decimal
	
This code currently only drops one fake in per image.

Not currently multicore. Submit multiple instances with different input files?

The script creates new directories:
```
|-Output_Images_V?/
|--PTF20*_P9fakes_V?.fits (Images with the fakes)
|
|-Results_V?/
|--Catalog/
|---PTF20*.cat (SExtractor Catalog)
|--Fake_Star_Catalog/
|---{input_name}_fakesAdded_V?.dat
```

The `{input_name}_fakesAdded_V?.dat` is an ascii file detailing the features of our fake star, our source star which was copied and scaled and other image parameters.
```
Its structure is as follows:
Column Name
0 Paths
1 Sourcex
2 Sourcey
3 a_source
4 dec_source
5 x_loc
6 y_loc
7 source_mag_auto
8 source_mag_best
9 flux_source
10 mag_fake
11 flux_fake
12 background
13 scaling_factor
14 PTFField
15 CCD
16 fbox1
17 fbox2
18 fbox3
19 fbox4
20 fbox5
21 fbox6
22 gain
23 readnoise
24 MOONILLF
25 MoonRA
26 MoonDec
27 AIRMASS
28 seeing
29 ELLIP
30 MEDSKY
31 SKYSIG)
32 zeropoint
33 LMT_MG
34 MJD
```