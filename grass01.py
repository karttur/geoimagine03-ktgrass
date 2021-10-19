'''
Created on 2 Apr 2021

@author: thomasgumbricht

Core copied from https://grass.osgeo.org/grass78/manuals/libpython/script.html#module-script.setup
'''


import os
import sys
import subprocess

os.environ['GRASSBIN']=r"/Applications/GRASS-7.8.app/Contents/MacOS/Grass"

os.environ['GISBASE'] = "/Applications/GRASS-7.8.app/Contents/Resources"

os.environ['GISRC']="/tmp/grass7-thomasgumbricht-61909/gisrc"

os.environ['PATH'] = "/Applications/GRASS-7.8.app/Contents/Resources/bin:/Applications/GRASS-7.8.app/Contents/Resources/scripts:/Users/thomasgumbricht/Library/GRASS/7.8/Modules/bin:/Users/thomasgumbricht/Library/GRASS/7.8/Modules/scripts:/usr/bin:/bin:/usr/sbin:/etc:/usr/lib"

os.environ['PYTHONPATH']="/Applications/GRASS-7.8.app/Contents/Resources/etc/python"

sys.path.append("/Applications/GRASS-7.8.app/Contents/Resources/etc/python")

import geoimagine.grass.script as grass

grass.run_command("g.list", flags="f", type="rast,vect")  

#grass.run_command("r.in.gdal", flags='e', input='/Volumes/karttur/ease2n/ESA/tiles/dem/x04y07/0/dem_copdem_x04y07_0_v01-90m.tif',output='tileDEMx04y07')  

#grass.run_command("r.slope.aspect", flags="en",elevation="tileDEMx04y07",slope='slope',aspect='aspect',pcurvature='pcurve',tcurvature='tcurve')

#r.geomorphon

#r.param.scale
## method elev, slope, aspect, profc, planc, longc, crosc, minic, maxic, feature
## size 3, 5, 7, 9, 11
### 
#grass.run_command("r.param.scale", input='tileDEMx04y07', output='slope_2x3x3', method='slope',size=5)

grass.run_command("r.param.scale", input='tileDEMx04y07', output='elev3x3', method='elev',size=3)

grass.run_command("r.out.gdal", input="elev3x3", output="/Volumes/karttur/ease2n/ESA/tiles/dem/x04y07/0/ps-elev_copdem_x04y07_0_v01-90m-3x3.tif", format="GTiff")

grass.run_command("g.remove", flags="f", type="raster", name="elev3x3")

grass.run_command("r.param.scale", input='tileDEMx04y07', output='elev5x5', method='elev',size=5)

grass.run_command("r.out.gdal", input="elev5x5", output="/Volumes/karttur/ease2n/ESA/tiles/dem/x04y07/0/ps-elev_copdem_x04y07_0_v01-90m-5x5.tif", format="GTiff")

grass.run_command("g.remove", flags="f", type="raster", name="elev5x5")

grass.run_command("r.param.scale", input='tileDEMx04y07', output='elev7x7', method='elev',size=7)

grass.run_command("r.out.gdal", input="elev7x7", output="/Volumes/karttur/ease2n/ESA/tiles/dem/x04y07/0/ps-elev_copdem_x04y07_0_v01-90m-7x7.tif", format="GTiff")

grass.run_command("g.remove", flags="f", type="raster", name="elev7x7")

grass.run_command("r.param.scale", input='tileDEMx04y07', output='elev9x9', method='elev',size=9)

grass.run_command("r.out.gdal", input="elev9x9", output="/Volumes/karttur/ease2n/ESA/tiles/dem/x04y07/0/ps-elev_copdem_x04y07_0_v01-90m-9x9.tif", format="GTiff")

grass.run_command("g.remove", flags="f", type="raster", name="elev9x9")

'''
grass.run_command("r.param.scale", input='tileDEMx04y07', output='slope3x3', method='slope',size=3)

grass.run_command("r.out.gdal", input="slope3x3", output="/Volumes/karttur/ease2n/ESA/tiles/dem/x04y07/0/ps-slope_copdem_x04y07_0_v01-90m-3x3.tif", format="GTiff")

grass.run_command("g.remove", flags="f", type="raster", name="slope3x3")

grass.run_command("r.param.scale", input='tileDEMx04y07', output='slope5x5', method='slope',size=5)

grass.run_command("r.out.gdal", input="slope5x5", output="/Volumes/karttur/ease2n/ESA/tiles/dem/x04y07/0/ps-slope_copdem_x04y07_0_v01-90m-5x5.tif", format="GTiff")

grass.run_command("g.remove", flags="f", type="raster", name="slope5x5")

grass.run_command("r.param.scale", input='tileDEMx04y07', output='slope7x7', method='slope',size=7)

grass.run_command("r.out.gdal", input="slope7x7", output="/Volumes/karttur/ease2n/ESA/tiles/dem/x04y07/0/ps-slope_copdem_x04y07_0_v01-90m-7x7.tif", format="GTiff")

grass.run_command("g.remove", flags="f", type="raster", name="slope7x7")

grass.run_command("r.param.scale", input='tileDEMx04y07', output='slope9x9', method='slope',size=9)

grass.run_command("r.out.gdal", input="slope9x9", output="/Volumes/karttur/ease2n/ESA/tiles/dem/x04y07/0/ps-slope_copdem_x04y07_0_v01-90m-9x9.tif", format="GTiff")

grass.run_command("g.remove", flags="f", type="raster", name="slope9x9")
'''

#input='tileDEMx04y07' output=profc37 method=profc size=11

#r.param.scale input='tileDEMx04y07' output='slope_2x3x3' method='slope' size=5

#grass.run_command("r.geomorphon", elevation='tileDEMx04y07', forms='geomorph3x3', search=3, skip=0)

#grass.run_command("r.out.gdal", input="geomorph3x3", output="/Volumes/karttur/ease2n/ESA/tiles/dem/x04y07/0/geomorpth_copdem_x04y07_0_v01-90m-3x3.tif", format="GTiff")

#grass.run_command("r.geomorphon", elevation='tileDEMx04y07', forms='geomorph5x5', search=5, skip=0)

#grass.run_command("r.out.gdal", input="geomorph5x5", output="/Volumes/karttur/ease2n/ESA/tiles/dem/x04y07/0/geomorpth_copdem_x04y07_0_v01-90m-5x5.tif", format="GTiff")

#grass.run_command("g.remove", flags="f", type="raster", name="geomorph5x5")
