# P9_Fakes
Generating fakes for the P9 hunt.

Based on the old supernova code for the PTF efficiencies.

##Requirements

-Take an image
-Specify an Ra, Dec and magnitude
-Lay a fake down at that location and magnitude

##How to use the code

Run the python script with a ascii file as an arguement.

The ascii file structure must be:

`path/to/image path/to/weight Ra_fake Dec_fake Mag_fake version_number`

The version number exists if you want to use the same new image multiple times

**Note:** Ra and Dec in decimal
	
This code currently only drops one fake in per image.

Not currently multicore. Submit multiple instances with different input files?