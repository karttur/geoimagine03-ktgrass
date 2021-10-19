'''
Created on 13 Apr 2021

@author: thomasgumbricht
'''

import os
import sys
import subprocess
from shutil import copy


import geoimagine.grass.script as grass

import geoimagine.support.karttur_dt as mj_dt

import geoimagine.gis as ktgis

class ProcessGRASSDEM():
    '''class for binding to GRASS DEM commands''' 
      
    def __init__(self, pp, session):
        '''
        '''
        
        self.session = session
                
        self.pp = pp  
        
        self.verbose = self.pp.process.verbose 
        
        self.session._SetVerbosity(self.verbose) 
        
        
                
        # Direct to subprocess
        if self.pp.process.processid == 'GrassDemFillDirTiles':
            
            self._GrassDemFillDirTiles()
            
        elif self.pp.process.processid == 'GrassDemHydroDemTiles':
            
            self.GrassDemHydroDemTiles()
                        
        else:
        
            print (self.pp.process.processid)

            exit('Process not recognized in processGRASSDEM')

    def _GrassDemFillDirTiles(self):
        '''
        '''

        srcComp = self.pp.srcCompL[0]
        
        dstComp = self.pp.dstCompL[0]
        
        if self.pp.process.parameters.asscript:
            
            today = mj_dt.Today()
            
            scriptFP = os.path.join('/Volumes',self.pp.dstPath.volume, self.pp.procsys.dstsystem, self.pp.dstCompD[dstComp].source,'tiles','demfillscript')
            
            if not os.path.exists(scriptFP):
                
                os.makedirs(scriptFP)
                
            scriptFN = 'grassDEM_%(comp)s-%(today)s.sh' %{'comp':dstComp, 'today':today}
            
            scriptFPN = os.path.join(scriptFP,scriptFN)
            
            self.scriptF = open(scriptFPN,'w')
            
            writeln = '# Script created by Kartturs GeoImagine Framework for GRASS processing %s, created %s\n' %(dstComp, today)
    
            self.scriptF.write(writeln)
            
        for datum in self.pp.srcPeriod.datumL:
                           
            for locus in self.pp.srcLayerD:
                
                self.tempFP = os.path.join(scriptFP,locus)
                
                if not os.path.exists(self.tempFP):
                
                    os.makedirs(self.tempFP)
                                    
                srcLayer = self.pp.srcLayerD[locus][datum][srcComp]
                
                if not os.path.exists(srcLayer.FPN):
                    
                    infostr = '    SKIPPING - input layer missing for ProcessGRASS\n        %s' %(srcLayer.FPN)
        
                    print (infostr)
                    
                    continue
                
                self.dstLayer = self.pp.dstLayerD[locus][datum][dstComp]
                
                if self.dstLayer._Exists() and not self.pp.process.overwrite:
                
                    continue
                
                
                if self.pp.process.parameters.mosaic:
                    
                    srcLayerFPN = srcLayer.FPN.replace(srcLayer.comp.ext, '-full.vrt')
                    
                    if not os.path.exists(srcLayerFPN):
                        
                        exitstr = 'EXITING, expecting a virtual mosaic as input for GDALDEMTiles'
                        
                        exit(exitstr)
                          
                else:
                    
                    srcLayerFPN = srcLayer.FPN
                    
                # Copy the source file to the destination
                
                copy(srcLayer.FPN, self.dstLayer.FPN)
                    
                srcPpatialRef, srcMeta = ktgis.GetRasterMetaData(srcLayerFPN)
                
                self.cellAreakm2 = float(srcMeta.cellsize)**2/1000000.00
                    
                if self.pp.process.parameters.asscript:
                
                    self._FillDirTilingAsScript(locus,srcLayerFPN,srcMeta)
                    
                else:
                    
                    NOTYET
                    
        self.scriptF.close()
        
        print ('Script file:',scriptFPN)
                    
                        
    def _FillDirTilingAsScript(self,locus,srcLayerFPN,srcMeta):
        '''
        '''
        
        paramsD = {'locus':locus,'src':srcLayerFPN}
        '''
        cmd = '# Pit filling  DEM\n'
        cmd += '# Created by the kartturs GeoImagine Framework \n\n'
        
        cmd += '# The tiling speeds up the process hundredfold compared to using the entire DEM.\n'
        cmd += '# Just make sure you have the parameter "tileCellOverlap" set to capture the pit' 
        cmd += '# sizes you want to fill (i.e. somewhat larger than the "FillDirCell" parameter). \n\n'
                
        cmd += '# To run this script you must have GRASS GIS open in a Terminal window session.\n'
        cmd += '# Change the script to be executable by the command:\n'
        cmd += '# chmod 755 %(fpn)s\n' %{'fpn': self.GRASS1shFPN}
        cmd += '# Then execute the command form your GRASS terminal session:\n'
        cmd += '# GRASS 7.x.3x ("region"):~ > %(fpn)s\n\n'%{'fpn': self.GRASS1shFPN}
        
        cmd += '# The script uses the filled DEM from stage0-step1 "inland_comp_DEM" as input.\n'
        cmd += '# To use another DEM as input manually edit the 1st r.mapcalc before the looping.\n\n'
        
        cmd += '# Copy the DEM to use as baseline for burning the DEM fill vectors\n'
        cmd += 'cp %(src)s %(dst)s\n\n' %{'src':self.params.inlandCompDEMFPN_s0, 'dst':self.params.hydroFillDEMFPN_s0}
        
        cmd += '# Set the region to the full map extent\n'
        cmd += 'g.region raster=%(dem)s\n\n' %{'dem': self.params.grassDEM}
        '''
        
        cmd = '# create a new mapset\n'
        
        cmd += 'g.mapset -c %(locus)s\n\n' %paramsD
        
        cmd += '# import the source DEM\n'
        
        cmd += 'r.in.gdal input=%(src)s output=iniDEM  --overwrite\n\n' %paramsD    
        
                            
        cmd += '# Set region after DEM\n'
        
        cmd += 'g.region raster=iniDEM\n\n'  
        
        cmd += '# Reclass sea level to null to speed up by skipping al tiles that are only null\n'
        
        cmd += 'r.mapcalc "nullDEM = if((iniDEM == 0), null(), iniDEM )" --overwrite\n\n'

        cmd += '# There is a bug in r.fill.dir and you must reclass all no data to a very low elevation (-1000)\n'
        
        cmd += 'r.mapcalc "srcDEM = if(isnull(iniDEM), -1000, iniDEM )" --overwrite\n\n'
        
        if self.pp.process.parameters.peaksize:
            
            cmd += '#r.watershed to get upstream areas for the inverted height suppresion\n'
            
            cmd += 'r.watershed -a elevation=srcDEM accumulation=MFD_upstream_filldir threshold=2000--overwrite\n\n'

        cmd += '# Loop over smaller tiles to run r.fill.dir\n\n'
        self.scriptF.write(cmd)
        
        tilesize = self.pp.process.parameters.tilesize
        
        overlap = self.pp.process.parameters.tileoverlap
        
        west, south, east, north = srcMeta.bounds
        
        #west, south, east, north = int(west), int(south), int(east), int(north)
        
        # featureD is a dict that hold the id and geometry of each polygon
        featureD = {}
        
        for c in range(0,int(srcMeta.cols/tilesize)+1):
            

            west0 = west+c*tilesize*srcMeta.cellsize - overlap*srcMeta.cellsize

            rW = max(west, west0)
            
            rE = min(east,tilesize*srcMeta.cellsize+west+c*tilesize*srcMeta.cellsize+overlap*srcMeta.cellsize)
            
            for r in range(0,int(srcMeta.lins/tilesize)+1):
            
                rS = max(south,south+r*tilesize*srcMeta.cellsize-overlap*srcMeta.cellsize)
                
                rN = min(north,tilesize*srcMeta.cellsize+south+r*tilesize*srcMeta.cellsize+overlap*srcMeta.cellsize)
                
                # Reset region
                
                cmd = '\ng.region -a n=%(n)f s=%(s)f e=%(e)f w=%(w)f\n\n' %{'n':rN, 's':rS, 'e':rE, 'w':rW}
                
                cmd += 'data="$(r.stats -p input=nullDEM null_value=\'null\')"\n'
                
                cmd += 'if echo "$data" | grep -q "null\s100"; then\n'
                
                cmd += '    echo "Null tile - skipping"\n'
                
                cmd += 'else\n'
                
                cmd += '    echo "Valid tile"\n'
                
                self.scriptF.write(cmd)
                
                tile_id = '%(c)d_%(r)d' %{'c':c,'r':r}
                
                featureD[tile_id] = {'id':tile_id, 'xmin':rW, 'xmax':rE, 'ymin':rS, 'ymax':rN}
                
                if self.pp.process.parameters.pitsize:
                
                    if self.pp.process.parameters.pitsize == 1:
                        
                        self._FillDirPoint(tile_id)
                        
                    else:
                        
                        self._FillDirArea(tile_id)
                      
                if self.pp.process.parameters.peaksize:
                    
                    if self.pp.parameters.peaksize >= 1:
                    
                        self._InvertedFillDirPoint(tile_id)
                        
                # Remove all layers - takes too much space
                cmd = '    g.remove -f type=raster pattern="*_%(tid)s"\n' %{'tid':tile_id}
                
                cmd += '    g.remove -f type=vector pattern="*_%(tid)s"\n' %{'tid':tile_id}
                  
                # end the if statement       
                cmd += 'fi\n'
                                       
                
                
                self.scriptF.write(cmd)

        self.scriptF.write('\n')
                
        cmd = 'g.region raster=iniDEM\n\n'
        
        self.scriptF.write(cmd)
               
        # Save the tiles used for the filldir 
        #self.WriteStage0DS(self.params.fillDEMtiles_s0, None, featureD)
        
    def _FillDirPoint(self, tile_id):
        '''
        '''
  
        ptshpfn = 'inland-fill-pt-%(tid)s.shp' %{'tid':tile_id}
                 
        ptshpfpn = os.path.join(self.tempFP, ptshpfn)
           
        # For single cells only with the -f flag    
        cmd = '    r.fill.dir -f input=srcDEM output=hydro_cellfill_DEM_%(tid)s direction=hydro_cellfill_draindir_%(tid)s areas=hydro_cellfill_problems_%(tid)s --overwrite\n' %{'tid':tile_id}
        
        # Check differences
        cmd += '    r.mapcalc "DEM_cellfill_diff_%(tid)s = hydro_cellfill_DEM_%(tid)s - srcDEM" --overwrite\n' %{'tid':tile_id}
    
        # get only diff areas
        cmd += '    r.mapcalc "inland_fill_cell_%(tid)s = if(DEM_cellfill_diff_%(tid)s != 0, 1, null() )" --overwrite\n' %{'tid':tile_id}
        
        # Convert to vector (area and pt)
        cmd += '    r.to.vect input=inland_fill_cell_%(tid)s output=inland_fill_pt_%(tid)s type=point --overwrite\n' %{'tid':tile_id}
        
        # Extract DEM fill value for points
        cmd += '    v.db.addcolumn map=inland_fill_pt_%(tid)s columns="filldem DOUBLE PRECISION" \n' %{'tid':tile_id}



        # Extract data from hydro_cellfill_DEM                 
        cmd += '    v.what.rast map=inland_fill_pt_%(tid)s column=filldem raster=hydro_cellfill_DEM_%(tid)s\n' %{'tid':tile_id}

        # Export always
        cmd += '    v.out.ogr input=inland_fill_pt_%(tid)s type=point format=ESRI_Shapefile output=%(pthshpfpn)s --overwrite\n\n' %{'tid':tile_id,'pthshpfpn':ptshpfpn}
           
        if self.pp.process.parameters.pitquery:
            
            cmd += '    GDAL_rasterize -a filldem -where "%(query)s" %(ptshpfpn)s %(dst)s\n' %{'ptshpfpn':ptshpfpn,'query': self.pp.process.parameters.pitquery, 'dst':self.dstLayer.FPN}
        
        else:
            
            cmd += '    GDAL_rasterize -a filldem  %(ptshpfpn)s %(dst)s\n' %{'ptshpfpn':ptshpfpn, 'dst':self.dstLayer.FPN}
         
        self.scriptF.write(cmd)
                      
    def _FillDirArea(self, tile_id): 
        '''
        '''
        areashpfn = 'inland-fill-area-%(tid)s.shp' %{'tid':tile_id}
                
        areashpfpn = os.path.join(self.tempFP, areashpfn) 
        
        # For all fill areas  without the _-f_ flag            
        cmd = '    r.fill.dir input=srcDEM output=hydro_areafill_DEM_%(tid)s direction=hydro_areafill_draindir_%(tid)s areas=hydro_areafill_problems_%(tid)s --overwrite\n' %{'tid':tile_id}
    
        # Check differences     
        cmd += '    r.mapcalc "DEM_areafill_diff_%(tid)s = hydro_areafill_DEM_%(tid)s - srcDEM" --overwrite\n' %{'tid':tile_id}
        
        # get only diff areas 
        cmd += '    r.mapcalc "inland_fill_area_%(tid)s = if(DEM_areafill_diff_%(tid)s != 0, 1, null() )" --overwrite\n' %{'tid':tile_id}
         
        # Convert to vector (area)
        cmd += '    r.to.vect input=inland_fill_area_%(tid)s output=inland_fill_area_%(tid)s type=area --overwrite\n' %{'tid':tile_id}
        
        # Get area of fill
        cmd += '    v.to.db map=inland_fill_area_%(tid)s type=centroid option=area columns=area_km2 units=kilometers\n' %{'tid':tile_id}

        cmd += '    v.db.addcolumn map=inland_fill_area_%(tid)s columns="filldem DOUBLE PRECISION" \n' %{'tid':tile_id}

        # Extract data from                    
        cmd += '    v.what.rast type=centroid map=inland_fill_area_%(tid)s column=filldem raster=hydro_areafill_DEM_%(tid)s\n' %{'tid':tile_id}

        # Export always, this layer to be used in GDAL_rasterize - simplest way to correct data  
        cmd += '    v.out.ogr input=inland_fill_area_%(tid)s type=area format=ESRI_Shapefile output=%(areahshpfpn)s --overwrite\n' %{'tid':tile_id,'areahshpfpn':areashpfpn}
        
        if self.pp.process.parameters.pitquery:
            
            cmd += '    GDAL_rasterize -a filldem -where "%(query)s" %(areahshpfpn)s %(dst)s\n' %{'areahshpfpn':areashpfpn,'query': self.pp.process.parameters.pitquery, 'dst':self.dstLayer.FPN}
        
        else:
            
            cmd += '    GDAL_rasterize -a filldem -where "area_km2 <= %(area)f" %(areahshpfpn)s %(dst)s\n' %{'areahshpfpn':areashpfpn,'area': self.cellAreakm2*self.pp.process.parameters.pitsize, 'dst':self.dstLayer.FPN}

        self.scriptF.write(cmd)
                
    def _InvertedFillDirPoint(self, tile_id):
        '''
        '''
  
        ptshpfn = 'inverted-fill-pt-%(tid)s.shp' %{'tid':tile_id}
                
        ptshpfpn = os.path.join(self.stage0datafp, ptshpfn)
 
        
        cmd = '    r.mapcalc "invertedDEM_%(tid)s = 10000-srcDEM" --overwrite\n' %{'tid':tile_id}
                
        cmd += '    r.fill.dir -f input=invertedDEM_%(tid)s output=inverted_cellfill_DEM_%(tid)s direction=inverted_cellfill_draindir_%(tid)s areas=inverted_cellfill_problems_%(tid)s --overwrite\n' %{'tid':tile_id}
    
        # Check differences
        cmd += '    r.mapcalc "inverted_cellfill_diff_%(tid)s = inverted_cellfill_DEM_%(tid)s - invertedDEM_%(tid)s" --overwrite\n' %{'tid':tile_id}
        
        # get only diff areas
        cmd += '    r.mapcalc "inverted_fill_cell_%(tid)s = if(inverted_cellfill_diff_%(tid)s != 0, 1, null() )" --overwrite\n' %{'tid':tile_id}
        
        # Convert to vector (area and pt)
        cmd += '    r.to.vect input=inverted_fill_cell_%(tid)s output=inverted_fill_pt_%(tid)s type=point --overwrite\n' %{'tid':tile_id}

        #Extract DEM fill value for points
        cmd += '    v.db.addcolumn map=inverted_fill_pt_%(tid)s\
         columns="invdem DOUBLE PRECISION, filldem DOUBLE PRECISION, nbupmax DOUBLE PRECISION, nbupq3 DOUBLE PRECISION, nbdemq1 DOUBLE PRECISION, nbdemmin DOUBLE PRECISION" \n' %{'tid':tile_id}

        # Extract data from hydro_cellfill_DEM             
        cmd += '    v.what.rast map=inverted_fill_pt_%(tid)s column=invdem raster=inverted_cellfill_DEM_%(tid)s\n' %{'tid':tile_id}
        
        cmd += '    v.db.update map=inverted_fill_pt_%(tid)s column=filldem query_column="10000-(invdem)"\n' %{'tid':tile_id}

        cmd += '    r.neighbors input=inland_comp_DEM selection=inverted_fill_cell_%(tid)s\
         output=dem_neighbor_quart1_%(tid)s,dem_neighbor_min_%(tid)s\
          method=quart1,minimum --overwrite\n' %{'tid':tile_id}
          
        cmd += '    r.neighbors input=MFD_upstream_filldir selection=inverted_fill_cell_%(tid)s\
         output=inverted_neighbor_max_%(tid)s,inverted_neighbor_quart3_%(tid)s\
          method=maximum,quart3 --overwrite\n' %{'tid':tile_id}
        
        cmd += '    v.what.rast map=inverted_fill_pt_%(tid)s column=nbupmax raster=inverted_neighbor_max_%(tid)s\n' %{'tid':tile_id}
                
        cmd += '    v.what.rast map=inverted_fill_pt_%(tid)s column=nbupq3 raster=inverted_neighbor_quart3_%(tid)s\n' %{'tid':tile_id}

        cmd += '    v.what.rast map=inverted_fill_pt_%(tid)s column=nbdemq1 raster=dem_neighbor_quart1_%(tid)s\n' %{'tid':tile_id}

        cmd += '    v.what.rast map=inverted_fill_pt_%(tid)s column=nbdemmin raster=dem_neighbor_min_%(tid)s\n' %{'tid':tile_id}

        cmd += '    v.out.ogr input=inverted_fill_pt_%(tid)s type=point format=ESRI_Shapefile output=%(pthshpfpn)s --overwrite\n\n' %{'tid':tile_id,'pthshpfpn':ptshpfpn}
     
        #cmd = 'GDAL_rasterize -a filldem %(src_datasource)s %(dst)s\n'  %{'src_datasource':ptshpfpn, 'dst':self.params.invertedFillDEMFPN_s0}           
    
        #self.GDALshF.write(cmd)
        
        self.GRASS1shF.write(cmd)
        
    def _InvertedFillDirArea(self, tile_id):
        '''
        '''
  
        areashpfn = 'inverted-fill-area-%(tid)s.shp' %{'tid':tile_id}
                
        areashpfpn = os.path.join(self.stage0datafp, areashpfn)
 
        
        cmd = '    r.mapcalc "invertedDEM_%(tid)s = 10000-srcDEM" --overwrite\n' %{'tid':tile_id}
                
        cmd += '    r.fill.dir input=invertedDEM_%(tid)s output=inverted_areafill_DEM_%(tid)s direction=inverted_areafill_draindir_%(tid)s areas=inverted_areafill_problems_%(tid)s --overwrite\n' %{'tid':tile_id}
    
        # Check differences
        cmd += '    r.mapcalc "inverted_areafill_diff_%(tid)s = inverted_areafill_DEM_%(tid)s - invertedDEM_%(tid)s" --overwrite\n' %{'tid':tile_id}
        
        # get only diff areas
        cmd += '    r.mapcalc "inverted_fill_area_%(tid)s = if(inverted_areafill_diff_%(tid)s != 0, 1, null() )" --overwrite\n' %{'tid':tile_id}
        
        # Convert to vector (area and pt)
        cmd += '    r.to.vect input=inverted_fill_area_%(tid)s output=inverted_fill_area_%(tid)s type=area --overwrite\n' %{'tid':tile_id}

        #Extract DEM fill value for points
        cmd += '    v.db.addcolumn map=inverted_fill_area_%(tid)s\
         columns="invdem DOUBLE PRECISION, filldem DOUBLE PRECISION, nbupmax DOUBLE PRECISION, nbupq3 DOUBLE PRECISION, nbdemq1 DOUBLE PRECISION, nbdemmin DOUBLE PRECISION" \n' %{'tid':tile_id}

        # Extract data from hydro_cellfill_DEM             
        cmd += '    v.what.rast map=inverted_fill_area_%(tid)s column=invdem raster=inverted_areafill_DEM_%(tid)s\n' %{'tid':tile_id}
        
        cmd += '    v.db.update map=inverted_fill_area_%(tid)s column=filldem query_column="10000-(invdem)"\n' %{'tid':tile_id}

        cmd += '    r.neighbors input=inland_comp_DEM \
         output=dem_neighbor_quart1_%(tid)s,dem_neighbor_min_%(tid)s\
          method=quart1,minimum --overwrite\n' %{'tid':tile_id}
          
        cmd += '    r.neighbors input=MFD_upstream_filldir \
         output=inverted_neighbor_max_%(tid)s,inverted_neighbor_quart3_%(tid)s\
          method=maximum,quart3 --overwrite\n' %{'tid':tile_id}
        
        cmd += '    v.what.rast map=inverted_fill_area_%(tid)s column=nbupmax raster=inverted_neighbor_max_%(tid)s\n' %{'tid':tile_id}
                
        cmd += '    v.what.rast map=inverted_fill_area_%(tid)s column=nbupq3 raster=inverted_neighbor_quart3_%(tid)s\n' %{'tid':tile_id}

        cmd += '    v.what.rast map=inverted_fill_area_%(tid)s column=nbdemq1 raster=dem_neighbor_quart1_%(tid)s\n' %{'tid':tile_id}

        cmd += '    v.what.rast map=inverted_fill_area_%(tid)s column=nbdemmin raster=dem_neighbor_min_%(tid)s\n' %{'tid':tile_id}

        cmd += '    v.out.ogr input=inverted_fill_area_%(tid)s type=point format=ESRI_Shapefile output=%(areahshpfpn)s --overwrite\n\n' %{'tid':tile_id,'areahshpfpn':areashpfpn}
     


    def _GrassDemHydroDemTiles(self):
        '''
        '''

        srcComp = self.pp.srcCompL[0]
        
        dstComp = self.pp.dstCompL[0]
        
        if self.pp.process.parameters.asscript:
            
            today = mj_dt.Today()
            
            scriptFP = os.path.join('/Volumes',self.pp.dstPath.volume, self.pp.procsys.dstsystem, self.pp.dstCompD[dstComp].source,'tiles','demfillscript')
            
            if not os.path.exists(scriptFP):
                
                os.makedirs(scriptFP)
                
            scriptFN = 'grassDEM-hydrodem_%(comp)s-%(today)s.sh' %{'comp':dstComp, 'today':today}
            
            scriptFPN = os.path.join(scriptFP,scriptFN)
            
            self.scriptF = open(scriptFPN,'w')
            
            writeln = '# Script created by Kartturs GeoImagine Framework for GRASS processing %s, created %s\n' %(dstComp, today)
    
            self.scriptF.write(writeln)
            
        for datum in self.pp.srcPeriod.datumL:
                           
            for locus in self.pp.srcLayerD:
                
                self.tempFP = os.path.join(scriptFP,locus)
                
                if not os.path.exists(self.tempFP):
                
                    os.makedirs(self.tempFP)
                                    
                srcLayer = self.pp.srcLayerD[locus][datum][srcComp]
                
                if not os.path.exists(srcLayer.FPN):
                    
                    infostr = '    SKIPPING - input layer missing for ProcessGRASS\n        %s' %(srcLayer.FPN)
        
                    print (infostr)
                    
                    continue
                
                self.dstLayer = self.pp.dstLayerD[locus][datum][dstComp]
                
                if self.dstLayer._Exists() and not self.pp.process.overwrite:
                
                    continue
                
                if self.pp.process.parameters.mosaic:
                    
                    srcLayerFPN = srcLayer.FPN.replace(srcLayer.comp.ext, '-full.vrt')
                    
                    if not os.path.exists(srcLayerFPN):
                        
                        exitstr = 'EXITING, expecting a virtual mosaic as input for GDALDEMTiles'
                        
                        exit(exitstr)
                          
                else:
                    
                    srcLayerFPN = srcLayer.FPN
                    
                # Copy the source file to the destination
                
                copy(srcLayer.FPN, self.dstLayer.FPN)
                    
                srcPpatialRef, srcMeta = ktgis.GetRasterMetaData(srcLayerFPN)
                
                #self.cellAreakm2 = float(srcMeta.cellsize)**2/1000000.00
                    
                if self.pp.process.parameters.asscript:
                
                    self._HydroDEMTilingAsScript(locus,srcLayerFPN,srcMeta)
                    
                else:
                    
                    NOTYET
                    
        self.scriptF.close()
        
        print ('Script file:',scriptFPN)
                    
                        
    def _HydroDEMTilingAsScript(self,locus,srcLayerFPN,srcMeta):
        '''
        '''
        
        paramsD = {'locus':locus,'src':srcLayerFPN}

        
        cmd = '# create a new mapset\n'
        
        cmd += 'g.mapset -c %(locus)s\n\n' %paramsD
        
        cmd += '# import the source DEM\n'
        
        cmd += 'r.in.gdal input=%(src)s output=iniDEM  --overwrite\n\n' %paramsD    
        
                            
        cmd += '# Set region after DEM\n'
        
        cmd += 'g.region raster=iniDEM\n\n'  
        
        cmd += '# Reclass sea level to null to speed up by skipping all tiles that are only null\n'
        
        cmd += 'r.mapcalc "nullDEM = if((iniDEM == 0), null(), iniDEM )" --overwrite\n\n'

        cmd += '# There is a bug in r.fill.dir and you must reclass all no data to a very low elevation (-1000)\n'
        
        cmd += 'r.mapcalc "srcDEM = if(isnull(iniDEM), -1000, iniDEM )" --overwrite\n\n'
        
    
        cmd += '# Loop over smaller tiles to run r.hydrodem\n\n'
        self.scriptF.write(cmd)
        
        tilesize = self.pp.process.parameters.tilesize
        
        overlap = self.pp.process.parameters.tileoverlap
        
        west, south, east, north = srcMeta.bounds
                
        # featureD is a dict that hold the id and geometry of each polygon
        featureD = {}
        
        for c in range(0,int(srcMeta.cols/tilesize)+1):
            
            west0 = west+c*tilesize*srcMeta.cellsize - overlap*srcMeta.cellsize

            rW = max(west, west0)
            
            rE = min(east,tilesize*srcMeta.cellsize+west+c*tilesize*srcMeta.cellsize+overlap*srcMeta.cellsize)
            
            for r in range(0,int(srcMeta.lins/tilesize)+1):
            
                rS = max(south,south+r*tilesize*srcMeta.cellsize-overlap*srcMeta.cellsize)
                
                rN = min(north,tilesize*srcMeta.cellsize+south+r*tilesize*srcMeta.cellsize+overlap*srcMeta.cellsize)
                
                # Reset region
                
                cmd = '\ng.region -a n=%(n)f s=%(s)f e=%(e)f w=%(w)f\n\n' %{'n':rN, 's':rS, 'e':rE, 'w':rW}
                
                cmd += 'data="$(r.stats -p input=nullDEM null_value=\'null\')"\n'
                
                cmd += 'if echo "$data" | grep -q "null\s100"; then\n'
                
                cmd += '    echo "Null tile - skipping"\n'
                
                cmd += 'else\n'
                
                cmd += '    echo "Valid tile"\n'
                
                self.scriptF.write(cmd)
                
                tile_id = '%(c)d_%(r)d' %{'c':c,'r':r}
                
                featureD[tile_id] = {'id':tile_id, 'xmin':rW, 'xmax':rE, 'ymin':rS, 'ymax':rN}
                
                # get the central tile
                
                rW = c*tilesize*srcMeta.cellsize
            
                rE = tilesize*srcMeta.cellsize+west+c*tilesize*srcMeta.cellsize
                                
                rS = south+r*tilesize*srcMeta.cellsize
                        
                rN = tilesize*srcMeta.cellsize+south+r*tilesize*srcMeta.cellsize
                
                self._HydroDEM(tile_id,rW,rE,rS,rN)
                       
                # Remove all layers - takes too much space
                cmd = '    g.remove -f type=raster pattern="*_%(tid)s"\n' %{'tid':tile_id}
                                  
                # end the if statement       
                cmd += 'fi\n'
                                       
                self.scriptF.write(cmd)

        self.scriptF.write('\n')
                
        cmd = 'g.region raster=iniDEM\n\n'
        
        self.scriptF.write(cmd)
               
        # Save the tiles used for the filldir 
        #self.WriteStage0DS(self.params.fillDEMtiles_s0, None, featureD)
        
    def _HydroDEM(self, tile_id, rW, rE, rS, rN):
        '''
        '''
  
        ptshpfn = 'hydrodem-%(tid)s.shp' %{'tid':tile_id}
                 
        ptshpfpn = os.path.join(self.tempFP, ptshpfn)
           
        # For single cells only with the -f flag    
        cmd = '    r.hydrodem -f input=srcDEM output=hydro_cellfill_DEM_%(tid)s direction=hydro_cellfill_draindir_%(tid)s areas=hydro_cellfill_problems_%(tid)s --overwrite\n' %{'tid':tile_id}
        
        #Reset region to the core tile and export
        
        cmd += '\ng.region -a n=%(n)f s=%(s)f e=%(e)f w=%(w)f\n\n' %{'n':rN, 's':rS, 'e':rE, 'w':rW}
        
        # Export 
        cmd += '    r.out.gdal input=inland_fill_pt_%(tid)s type=point format=ESRI_Shapefile output=%(pthshpfpn)s --overwrite\n\n' %{'tid':tile_id,'pthshpfpn':ptshpfpn}
                   
        self.scriptF.write(cmd)
                      
                        
    