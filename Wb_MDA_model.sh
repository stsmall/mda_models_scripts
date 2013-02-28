#!/bin/sh

### time v. effectiveness







##values for 2008
#---------
# No MDA
ms 100 10000 -t 32 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/Nomda.txt

#30% effective
ms 100 10000 -t 32 -G 33373.0942 -eG .0000374063 0.0 -eN .0000374063 .16807 -eN 0.0000400781 .2401 -eN 0.00004275 .343 -eN 0.0000454219 .49 -eN .0000480938 .7 -eN .0000507657 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/30mda.txt

#50% effective
ms 100 10000 -t 32 -G 64855.8765 -eG .0000374063 0.0 -eN .0000374063 .03125 -eN 0.0000400781 .0625 -eN 0.00004275 .125 -eN 0.0000454219 .25 -eN .0000480938 .5 -eN .0000507657 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/50mda.txt
#90% effective
ms 100 10000 -t 32 -G 215446.558 -eG .0000374063 0.0 -eN .0000374063 .00001 -eN 0.0000400781 .0001 -eN 0.00004275 .001 -eN 0.0000454219 .01 -eN .0000480938 .1 -eN .0000507657 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/90mda.txt



#---------
##value for 1998

# No MDA
#ms 100 10 -t 15 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats 

#30% effective
ms 400 10000 -t 32 -eN 0.000178571429 .2401 -eN 0.000191326531 .343 -eN 0.000204081633 .49 -eN 0.000216836735 .7 -eN .000229 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/30mda_1998.txt

#50% effective
ms 400 10000 -t 32 -eN 0.000178571429 .0625 -eN 0.000191326531 .125 -eN 0.000204081633 .25 -eN 0.000216836735 .5 -eN .000229 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/50mda_1998.txt

#90% effective
ms 400 10000 -t 32 -eN 0.000178571429 .0001 -eN 0.000191326531 .001 -eN 0.000204081633 .01 -eN 0.000216836735 .1 -eN .000229 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/90mda_1998.txt

#--------
##values for 1997
# No MDA
#ms 100 10 -t 15 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats

#30% effective
ms 400 10000 -t 32 -eN 0.000191326531 .343 -eN 0.000204081633 .49 -eN 0.000216836735 .7 -eN .000229 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/30mda_1997.txt

#50% effective
ms 400 10000 -t 32 -eN 0.000191326531 .125 -eN 0.000204081633 .25 -eN 0.000216836735 .5 -eN .000229 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/50mda_1997.txt

#90% effective
ms 400 10000 -t 32 -eN 0.000191326531 .001 -eN 0.000204081633 .01 -eN 0.000216836735 .1 -eN .000229 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/90mda_1997.txt

#-------
##values for 1996
# No MDA
#ms 100 10 -t 15 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats 

#30% effective
ms 400 10000 -t 32 -eN 0.000204081633 .49 -eN 0.000216836735 .7 -eN .000229 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/30mda_1996.txt

#50% effective
ms 400 10000 -t 32 -eN 0.000204081633 .25 -eN 0.000216836735 .5 -eN .000229 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/50mda_1996.txt

#90% effective
ms 400 10000 -t 32  -eN 0.000204081633 .01 -eN 0.000216836735 .1 -eN .000229 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/90mda_1996.txt

#------
##values for 1995
# No MDA
#ms 100 1000 -t 15 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats

#30% effective
ms 400 10000 -t 32 -eN 0.000216836735 .7 -eN .000229 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/30mda_1995.txt

#50% effective
ms 400 10000 -t 32 -eN 0.000216836735 .5 -eN .000229 1.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/50mda_1995.txt

#90% effective
ms 400 10000 -t 32 -eN 0.000216836735 .1 -eN .000229 2.0 | /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/msdir/sample_stats >/Volumes/home/Users/stsmall/Desktop/90mda_1995.txt
