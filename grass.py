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

os.environ['PATH'] = "/Applications/GRASS-7.8.app/Contents/Resources/bin:/Applications/GRASS-7.8.app/Contents/Resources/scripts:/Users/thomasgumbricht/Library/GRASS/7.8/Modules/bin:/Users/thomasgumbricht/Library/GRASS/7.8/Modules/scripts:/usr/bin:/bin:/usr/sbin:/etc:/usr/lib"

os.environ['PYTHONPATH']="/Applications/GRASS-7.8.app/Contents/Resources/etc/python"

sys.path.append("/Applications/GRASS-7.8.app/Contents/Resources/etc/python")

# import (some) GRASS Python bindings
#from grass_session import Session

from grass_session import Session as gsession

from geoimagine.grass.script import core as gcore

import geoimagine.grass.script.setup as gsetup

class ProcessGRASSv03():
    
    def __init__(self):
        
        self.gisdb = os.path.join(os.path.expanduser("~"), "GRASSDATA")
        
        self.location = 'ease2n'
        
        self.mapset = 'PERMANENT'
        
        #with gsession(gisdb=self.gisdb, self.location,
        #      create_opts="EPSG:4326"):
            
        with Session(gisdb=self.gisdb, location=self.location, mapset="test",
              create_opts=""):
            
            print(gcore.parse_command("g.gisenv", flags="s"))
            


class ProcessGRASSv02():
    
    def __init__(self):
        
        gisbase = "/Applications/GRASS-7.8.app/Contents/Resources"
        
        grass_pydir = "/Applications/GRASS-7.8.app/Contents/Resources/etc/python"

        os.environ['GRASSBIN']=r"/Applications/GRASS-7.8.app/Contents/MacOS/Grass"

        os.environ['GISBASE'] = "/Applications/GRASS-7.8.app/Contents/Resources"
        
        os.environ['PATH'] = "/Applications/GRASS-7.8.app/Contents/Resources/bin:/Applications/GRASS-7.8.app/Contents/Resources/scripts:/Users/thomasgumbricht/Library/GRASS/7.8/Modules/bin:/Users/thomasgumbricht/Library/GRASS/7.8/Modules/scripts:/usr/bin:/bin:/usr/sbin:/etc:/usr/lib"
        
        os.environ['PYTHONPATH']="/Applications/GRASS-7.8.app/Contents/Resources/etc/python"
        
        sys.path.append("/Applications/GRASS-7.8.app/Contents/Resources/etc/python")

        grass7bin = r'/Applications/GRASS-7.8.app/Contents/MacOS/Grass'

        
        startcmd = [grass7bin, '--config', 'path']
        
        try:
            p = subprocess.Popen(startcmd, shell=False,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            out, err = p.communicate()
            
        except OSError as error:
            
            sys.exit("ERROR: Cannot find GRASS GIS start script"
                     " {cmd}: {error}".format(cmd=startcmd[0], error=error))
            
        if p.returncode != 0:
            
            sys.exit("ERROR: Issues running GRASS GIS start script"
                     " {cmd}: {error}"
                     .format(cmd=' '.join(startcmd), error=err))
            
        print ('out',out)
        print ('err',err)
        
        #gisbase = out.strip(os.linesep)
        
        #gisbase= os.path.join(os.path.expanduser("~"), "grassdata")
        
        os.environ['GISBASE'] = gisbase
        #grass_pydir = os.path.join(gisbase, "etc", "python")
        sys.path.append(grass_pydir)
        
        import geoimagine.grass.script.setup as gsetup
        
        import geoimagine.grass.script as gscript
        
        #import /Applications/GRASS-7.8.app/Contents/Resources/etc/python/grass.script.setup as gsetup
        
        gisdb = os.path.join(os.path.expanduser("~"), "Documents/grassdata")
        location = "nc_spm_08"
        mapset = "user1"
        
        rcfile = gsetup.init(gisbase, gisdb, location, mapset)
        
    def _PrintCurrent(self):
        
        print ('in _PrintCurrent')

        # example calls
        self.gscript.message('Current GRASS GIS 7 environment:')
        print ( self.gscript.gisenv() )
        
    def _ListRaster(self):

        gscript.message('Available raster maps:')
        for rast in gscript.list_strings(type='raster'):
            print ( rast )
            
    def _ListVectors(self):

        gscript.message('Available vector maps:')
        for vect in gscript.list_strings(type='vector'):
            print ( vect )

    def _Quit(self):
        # clean up at the end
        gsetup.cleanup()


class ProcessGRASSv01():
    
    def __init__(self):
        '''
        self.session = session
        self.verbose = verbose
        self.process = process
        '''
        
        # define GRASS Database
        # add your path to grassdata (GRASS GIS database) directory
        gisdb = os.path.join(os.path.expanduser("~"), "grassdata")
        # the following path is the default path on MS Windows
        # gisdb = os.path.join(os.path.expanduser("~"), "Documents/grassdata")
        
        if not os.path.exists(gisdb):
            
            gisdb = os.path.join('/Volumes','grassdata')

        # specify (existing) Location and Mapset
        location = self.process.pp.procesys.dstsystem
        mapset = "PERMANENT"

        # path to the GRASS GIS launch script
        # we assume that the GRASS GIS start script is available and on PATH
        # query GRASS itself for its GISBASE
        # (with fixes for specific platforms)
        # needs to be edited by the user
        
        grass7bin = 'grass78'
        if sys.platform.startswith('win'):
            # MS Windows
            grass7bin = r'C:\OSGeo4Win\grass78.bat'
            # uncomment when using standalone WinGRASS installer
            # grass7bin = r'C:\Program Files (x86)\GRASS GIS 7.8.0\grass78.bat'
            # this can be avoided if GRASS executable is added to PATH
        elif sys.platform == 'darwin':
            # Mac OS X
            # TODO: this have to be checked, maybe unix way is good enough
            grass7bin = '/Applications/GRASS-7.8.app/'

        # query GRASS GIS itself for its GISBASE
        startcmd = [grass7bin, '--config', 'path']
        try:
            p = subprocess.Popen(startcmd, shell=False,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
        except OSError as error:
            sys.exit("ERROR: Cannot find GRASS GIS start script"
                     " {cmd}: {error}".format(cmd=startcmd[0], error=error))
        if p.returncode != 0:
            sys.exit("ERROR: Issues running GRASS GIS start script"
                     " {cmd}: {error}"
                     .format(cmd=' '.join(startcmd), error=err))
        gisbase = out.strip(os.linesep)

        # set GISBASE environment variable
        os.environ['GISBASE'] = gisbase

        # define GRASS-Python environment
        grass_pydir = os.path.join(gisbase, "etc", "python")
        sys.path.append(grass_pydir)

        # import (some) GRASS Python bindings
        #import grass.script as gscript
        #import grass.script.setup as gsetup

        # launch session
        rcfile = gsetup.init(gisbase, gisdb, location, mapset)
        
    def _PrintCurrent(self):

        # example calls
        gscript.message('Current GRASS GIS 7 environment:')
        print ( gscript.gisenv() )
        
    def _ListRaster(self):

        gscript.message('Available raster maps:')
        for rast in gscript.list_strings(type='raster'):
            print ( rast )
            
    def _ListVectors(self):

        gscript.message('Available vector maps:')
        for vect in gscript.list_strings(type='vector'):
            print ( vect )

    def _Quit(self):
        # clean up at the end
        gsetup.cleanup()
        
        
if __name__ == "__main__":
    
    '''
    grass = ProcessGRASSv02()
    
    print ('printing')
    
    grass._PrintCurrent
    '''
    grass = ProcessGRASSv03()