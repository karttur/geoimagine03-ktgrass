"""
Created 4 Apr 2021

GRASS
==========================================

Package belonging to Karttur´s GeoImagine Framework.

Author
------
Thomas Gumbricht (thomas.gumbricht@karttur.com)

"""

import os
import sys

from .version import __version__, VERSION, metadataD


os.environ['GRASSBIN']=r"/Applications/GRASS-7.8.app/Contents/MacOS/Grass"

os.environ['GISBASE'] = "/Applications/GRASS-7.8.app/Contents/Resources"
        
os.environ['GISRC']="/tmp/grass7-thomasgumbricht-61909/gisrc"
        
os.environ['PATH'] = "/Applications/GRASS-7.8.app/Contents/Resources/bin:/Applications/GRASS-7.8.app/Contents/Resources/scripts:/Users/thomasgumbricht/Library/GRASS/7.8/Modules/bin:/Users/thomasgumbricht/Library/GRASS/7.8/Modules/scripts:/usr/bin:/bin:/usr/sbin:/etc:/usr/lib"
        
os.environ['PYTHONPATH']="/Applications/GRASS-7.8.app/Contents/Resources/etc/python"
        
sys.path.append("/Applications/GRASS-7.8.app/Contents/Resources/etc/python")

from .grass02 import ProcessGRASS

from .grassdem import ProcessGRASSDEM
