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
            
            self._GrassDemHydroDemTiles()
            
        elif self.pp.process.processid == 'GrassDemReInterpolateTiles':
            
            self._GrassDemReInterpolateTiles()
                        
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
                
            scriptFN = 'grassDEM-filldir_%(comp)s-%(today)s.sh' %{'comp':dstComp, 'today':today}
            
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
                        
                        exitstr = 'EXITING, expecting a virtual mosaic as input for GDALDEMTiles\n    %s' %(srcLayerFPN)
                            
                        print (exitstr)
                            
                        SNULLE
                          
                else:
                    
                    srcLayerFPN = srcLayer.FPN
                    
                # Copy the source file to the destination
                if os.path.exists(self.dstLayer.FPN) and self.pp.process.parameters.superimpose:
                    
                    pass
                
                else:
                
                    copy(srcLayer.FPN, self.dstLayer.FPN)
                    
                srcSpatialRef, srcMeta = ktgis.GetRasterMetaData(srcLayerFPN)
                
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
        
        self._InitiateDemCorrection(paramsD)
    
        tilesize = self.pp.process.parameters.tilesize
        
        overlap = self.pp.process.parameters.tileoverlap
        
        west, south, east, north = srcMeta.bounds
                
        # featureD is a dict that hold the id and geometry of each polygon
        featureD = {}
        
        for c in range(0,int(srcMeta.cols/tilesize)+1):
            
            rW = max(west, west+c*tilesize*srcMeta.cellsize - overlap*srcMeta.cellsize)
            
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
                
                rWt = west+c*tilesize*srcMeta.cellsize
            
                rEt = tilesize*srcMeta.cellsize+west+c*tilesize*srcMeta.cellsize
                                
                rSt = south+r*tilesize*srcMeta.cellsize
                        
                rNt = tilesize*srcMeta.cellsize+south+r*tilesize*srcMeta.cellsize
                                
                if self.pp.process.parameters.pitsize:
                                    
                    if self.pp.process.parameters.pitsize == 1:
                        
                        self._FillDirPoint(tile_id, rWt, rEt, rSt, rNt)
                        
                    else:
                        
                        self._FillDirArea(tile_id, rWt, rEt, rSt, rNt)
                      
                if self.pp.process.parameters.peaks:
                                                
                    self._InvertedFillDirPoint(tile_id, rWt, rEt, rSt, rNt)
                                 
                # Remove all layers 
                cmd = '    g.remove -f type=raster pattern="*_%(tid)s"\n' %{'tid':tile_id}
                
                cmd += '    g.remove -f type=vector pattern="*_%(tid)s"\n' %{'tid':tile_id}
                  
                # end the if statement       
                cmd += 'fi\n'
                                       
                self.scriptF.write(cmd)

        self.scriptF.write('\n')
                
        cmd = 'g.region raster=iniDEM\n\n'
        
        self.scriptF.write(cmd)
               
    def _InitiateDemCorrection(self,paramsD):
        '''
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
        
        cmd += '#r.watershed to get upstream areas and streams for queries\n'
            
        cmd += 'r.watershed -a elevation=nullDEM accumulation=MFD_updrain stream=MFD_stream length_slope=MFD_length_slope threshold=2000--overwrite\n\n'

        cmd += '# Loop over smaller tiles to run r.hydrodem\n\n'
        
        self.scriptF.write(cmd)
        
    def _InnerLoop(self,paramsD,rNt, rSt, rEt,rWt, single=False, peaks=False, hydrodem=False):
        '''
        '''
        
        #Reset region to the core tile and export
        
        cmd = '    g.region -a n=%(n)f s=%(s)f e=%(e)f w=%(w)f\n\n' %{'n':rNt, 's':rSt, 'e':rEt, 'w':rWt}
        
        # Check differences
        if peaks:
            
            cmd += '    r.mapcalc "diff_%(tid)s = %(output)s - invertedDEM_%(tid)s" --overwrite\n' %paramsD
        
        else:
                 
            cmd += '    r.mapcalc "diff_%(tid)s = %(output)s - srcDEM" --overwrite\n' %paramsD
        
        # get only diff areas 
        if paramsD['diffthreshold'] == 0:
            
            cmd += '    r.mapcalc "fillarea_%(tid)s = if( diff_%(tid)s != 0, 1, null() )" --overwrite\n' %paramsD

        else:
            
            cmd += '    r.mapcalc "fillarea_%(tid)s = if( abs(diff_%(tid)s) >= %(diffthreshold)f, 1, null() )" --overwrite\n' %paramsD
         
        if not single:
            
            # Convert to vector (area)
            cmd += '    r.to.vect input=fillarea_%(tid)s output=fillarea_%(tid)s type=area --overwrite\n' %paramsD
            
            # Get area of fill
            cmd += '    v.to.db map=fillarea_%(tid)s type=centroid option=area columns=area_km2 units=kilometers\n' %paramsD
    
            cmd += '    v.db.addcolumn map=fillarea_%(tid)s columns="filldem DOUBLE PRECISION, area_cells INT" \n' %paramsD
    
            # Extract data from raster                   
            cmd += '    v.what.rast type=centroid map=fillarea_%(tid)s column=filldem raster=%(output)s\n' %paramsD
    
            # Calculate cell area from km2                    
            cmd += '    v.db.update fillarea_%(tid)s column=area_cells qcol="area_km2/%(cellarea)f"\n' %paramsD
            
            if hydrodem:
                
                cmd += '    v.rast.stats map=fillarea_%(tid)s raster=MFD_updrain column_prefix=u method=maximum\n' %paramsD
            
            
            # Export always, this layer to be used in GDAL_rasterize - simplest way to correct data  
            cmd += '    v.out.ogr input=fillarea_%(tid)s type=area format=ESRI_Shapefile output=%(areashpfpn)s --overwrite\n' %paramsD
                
            if hydrodem:
                
                if paramsD['query']:
            
                    cmd += '    GDAL_rasterize -a %(burnattribute)s -where "%(query)s" %(areashpfpn)s %(dst)s\n' %paramsD
        
                else:
                    
                    cmd += '    GDAL_rasterize -a %(burnattribute)s %(areashpfpn)s %(dst)s\n' %paramsD
                
                self.scriptF.write(cmd)
                
                return
  
        # Convert to vector (point)
        cmd += '    r.to.vect input=fillarea_%(tid)s output=fillpt_%(tid)s type=point --overwrite\n' %paramsD

        if peaks:
            
            # add database columns
            cmd += '    v.db.addcolumn map=fillpt_%(tid)s\
             columns="area_cells INT, area_km2 DOUBLE PRECISION, invdem DOUBLE PRECISION,\
             filldem DOUBLE PRECISION,\
             nbupmax DOUBLE PRECISION,\
             nbupq3 DOUBLE PRECISION,\
             nbdemq1 DOUBLE PRECISION"\n' %paramsD
              
               
            # Extract data from filldir             
            cmd += '    v.what.rast map=fillpt_%(tid)s column=invdem raster=%(output)s\n' %paramsD
            
            cmd += '    v.db.update map=fillpt_%(tid)s column=filldem query_column="10000-(invdem)"\n' %paramsD
    
            # neighbor analysis on DEM only around diff 
            cmd += '    r.neighbors input=srcDEM selection=diff_%(tid)s\
             output=dem_neighbor_quart1_%(tid)s,dem_neighbor_min_%(tid)s\
              method=quart1,minimum --overwrite\n' %paramsD
     
            cmd += '    r.neighbors input=MFD_updrain selection=diff_%(tid)s\
             output=inverted_neighbor_max_%(tid)s,inverted_neighbor_quart3_%(tid)s\
              method=maximum,quart3 --overwrite\n' %paramsD
              

            
            cmd += '    v.what.rast map=fillpt_%(tid)s column=nbupmax raster=inverted_neighbor_max_%(tid)s\n' %paramsD
                    
            cmd += '    v.what.rast map=fillpt_%(tid)s column=nbupq3 raster=inverted_neighbor_quart3_%(tid)s\n' %paramsD
    
            cmd += '    v.what.rast map=fillpt_%(tid)s column=nbdemq1 raster=dem_neighbor_quart1_%(tid)s\n' %paramsD
    
            #cmd += '    v.what.rast map=fillpt_%(tid)s column=nbdemmin raster=dem_neighbor_min_%(tid)s\n' %paramsD

            #cmd += '    v.what.rast map=fillpt_%(tid)s column=nbstream raster=inverted_stream_%(tid)s\n' %paramsD

        else:
            # add database columns
            cmd += '    v.db.addcolumn map=fillpt_%(tid)s\
             columns="area_cells INT, area_km2 DOUBLE PRECISION, filldem DOUBLE PRECISION,\
              updrain DOUBLE PRECISION" \n' %paramsD
        
            # Extract the data
            if not single: 
                
                # Convert vector area to raster
                cmd += '    v.to.rast input=fillarea_%(tid)s output=area_km2_%(tid)s use=attr attribute_column=area_km2 memory=%(mem)d --overwrite\n' %paramsD
                
                
                cmd += '    v.what.rast map=fillpt_%(tid)s column=area_km2 raster=area_km2_%(tid)s\n' %paramsD
        
            cmd += '    v.what.rast map=fillpt_%(tid)s column=filldem raster=%(output)s\n' %paramsD
            
            if paramsD['query']:
            
                cmd += '    v.what.rast map=fillpt_%(tid)s column=updrain raster=MFD_updrain\n' %paramsD
    
                #cmd += '    v.what.rast map=fillpt_%(tid)s column=stream raster=MFD_stream\n' %paramsD
    
                #cmd += '    v.what.rast map=fillpt_%(tid)s column=slopelen raster=MFD_length_slope\n' %paramsD
    
                # Calculate cell area from km2 
                if single:
                                       
                    cmd += '    v.db.update fillpt_%(tid)s column=area_cells value=1\n' %paramsD
                
                else:
                    
                    cmd += '    v.db.update fillpt_%(tid)s column=area_cells qcol="area_km2/%(cellarea)f"\n' %paramsD
    
        # Export  
        cmd += '    v.out.ogr input=fillpt_%(tid)s type=point format=ESRI_Shapefile output=%(ptshpfpn)s --overwrite\n' %paramsD
    
        # Export 
        #cmd += '    r.out.gdal input=%(output)s format=GTiff output=%(geotif)s --overwrite\n\n' %paramsD
            
        if paramsD['query']:
            
            cmd += '    GDAL_rasterize -a %(burnattribute)s -where "%(query)s" %(ptshpfpn)s %(dst)s\n' %paramsD

        else:
            
            cmd += '    GDAL_rasterize -a %(burnattribute)s %(ptshpfpn)s %(dst)s\n' %paramsD
        
        self.scriptF.write(cmd)
                    
    def _FillDirPoint(self, tile_id, rWt, rEt, rSt, rNt):
        '''
        '''
        
        grasslayer = 'filldir_%(tid)s' %{'tid':tile_id}
        
        geotiffn = 'filldir_pt_%(tid)s.tif' %{'tid':tile_id}
               
        geotiffpn = os.path.join(self.tempFP, geotiffn)
        
        ptshpfn = 'filldir-pt-%(tid)s.shp' %{'tid':tile_id}
                
        ptshpfpn = os.path.join(self.tempFP, ptshpfn)
        
        paramsD = {'flags':"", 'output': grasslayer, 'geotif':geotiffpn, 
                   'diffthreshold':0,
                   'tid':tile_id,
                   'ptshpfpn':ptshpfpn,
                   'cellarea': self.cellAreakm2,
                   'query': self.pp.process.parameters.pitquery,
                   'burnattribute': self.pp.process.parameters.pitburnattribute,
                   'dst': self.dstLayer.FPN,
                   'mem':self.pp.process.parameters.memory}
        
        # For all fill areas  with the -f flag            
        cmd = '    r.fill.dir -f input=srcDEM output=%(output)s direction=hydro_ptfill_draindir_%(tid)s areas=hydro_ptfill_problems_%(tid)s --overwrite\n' %paramsD
    
        self.scriptF.write(cmd)
        
        self._InnerLoop(paramsD,rNt, rSt, rEt,rWt, single=True)
                       
    def _FillDirArea(self, tile_id, rWt, rEt, rSt, rNt):
        '''
        '''
        
        grasslayer = 'filldir_%(tid)s' %{'tid':tile_id}
        
        geotiffn = 'filldir_area_%(tid)s.tif' %{'tid':tile_id}
               
        geotiffpn = os.path.join(self.tempFP, geotiffn)
        
        areashpfn = 'filldir-area-%(tid)s.shp' %{'tid':tile_id}
                
        areashpfpn = os.path.join(self.tempFP, areashpfn) 
        
        ptshpfn = 'filldir-area-pt-%(tid)s.shp' %{'tid':tile_id}
                
        ptshpfpn = os.path.join(self.tempFP, ptshpfn)
        
        paramsD = {'flags':"", 'output': grasslayer, 'geotif':geotiffpn, 
                   'diffthreshold':0,
                   'tid':tile_id,
                   'areashpfpn':areashpfpn,
                   'ptshpfpn':ptshpfpn,
                   'cellarea': self.cellAreakm2,
                   'query': self.pp.process.parameters.pitquery,
                   'burnattribute': self.pp.process.parameters.pitburnattribute,
                   'dst': self.dstLayer.FPN,
                   'mem':self.pp.process.parameters.memory}
        
        # For all fill areas  without the _-f_ flag            
        cmd = '    r.fill.dir input=srcDEM output=%(output)s direction=hydro_areafill_draindir_%(tid)s areas=hydro_areafill_problems_%(tid)s --overwrite\n' %paramsD
    
        self.scriptF.write(cmd)
        
        self._InnerLoop(paramsD,rNt, rSt, rEt,rWt, single=False)
                
    def _InvertedFillDirPoint(self, tile_id, rWt, rEt, rSt, rNt):
        '''
        '''
  
        grasslayer = 'flattenpeak_%(tid)s' %{'tid':tile_id}
        
        geotiffn = 'flattenpeak_pt_%(tid)s.tif' %{'tid':tile_id}
               
        geotiffpn = os.path.join(self.tempFP, geotiffn)
        
        ptshpfn = 'flattenpeak-pt-%(tid)s.shp' %{'tid':tile_id}
                
        ptshpfpn = os.path.join(self.tempFP, ptshpfn)
        
        paramsD = {'flags':"", 'output': grasslayer, 'geotif':geotiffpn, 
                   'diffthreshold':0,
                   'tid':tile_id,
                   'ptshpfpn':ptshpfpn,
                   'cellarea': self.cellAreakm2,
                   'query': self.pp.process.parameters.peakquery,
                   'burnattribute': self.pp.process.parameters.peakburnattribute,
                   'dst': self.dstLayer.FPN,
                   'mem':self.pp.process.parameters.memory}
                
        # r.fill.dir with -f flag  
        
        cmd = '    r.mapcalc "invertedDEM_%(tid)s = 10000-srcDEM" --overwrite\n' %{'tid':tile_id}

        cmd += '    r.fill.dir -f input=invertedDEM_%(tid)s output=%(output)s direction=hydro_pkflat_draindir_%(tid)s areas=hydro_pkflat_problems_%(tid)s --overwrite\n' %paramsD
    
        self.scriptF.write(cmd)
        
        self._InnerLoop(paramsD,rNt, rSt, rEt,rWt, single=True, peaks=True)
   
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
                        
                        exitstr = 'EXITING, expecting a virtual mosaic as input for GDALDEMTiles\n    %s' %(srcLayerFPN)
                            
                        print (exitstr)
                            
                        SNULLE
                          
                else:
                    
                    srcLayerFPN = srcLayer.FPN

                # Copy the source file to the destination
                if os.path.exists(self.dstLayer.FPN) and self.pp.process.parameters.superimpose:
                    
                    pass
                
                else:
                
                    copy(srcLayer.FPN, self.dstLayer.FPN)
                    
                # Get the metadata of the source file
                srcSpatialRef, srcMeta = ktgis.GetRasterMetaData(srcLayerFPN)
                
                self.cellAreakm2 = float(srcMeta.cellsize)**2/1000000.00
                    
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
        
        self._InitiateDemCorrection(paramsD)
        
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
                
                # get the central tile coordinates
                
                rWt = west+c*tilesize*srcMeta.cellsize
            
                rEt = tilesize*srcMeta.cellsize+west+c*tilesize*srcMeta.cellsize
                                
                rSt = south+r*tilesize*srcMeta.cellsize
                        
                rNt = tilesize*srcMeta.cellsize+south+r*tilesize*srcMeta.cellsize
                
                self._HydroDEM(tile_id,rWt,rEt,rSt,rNt)
                       
                # Remove all layers - takes too much space
                cmd = '    g.remove -f type=raster pattern="*_%(tid)s"\n' %{'tid':tile_id}
                
                cmd += '    g.remove -f type=vector pattern="*_%(tid)s"\n' %{'tid':tile_id}
                                  
                # end the if statement       
                cmd += 'fi\n'
                                       
                self.scriptF.write(cmd)

        self.scriptF.write('\n')
                
        cmd = 'g.region raster=iniDEM\n\n'
        
        self.scriptF.write(cmd)
                           
    def _HydroDEM(self, tile_id, rWt, rEt, rSt, rNt):
        '''
        '''
  
        grasslayer = 'hydrodem_%(tid)s' %{'tid':tile_id}
        
        areashpfn = 'hydrodem_area_%(tid)s.shp' %{'tid':tile_id}
        
        areashpfpn = os.path.join(self.tempFP, areashpfn)
        
        ptshpfn = 'hydrodem_pt_%(tid)s.shp' %{'tid':tile_id}
        
        ptshpfpn = os.path.join(self.tempFP, ptshpfn)
        
        geotiffn = 'hydrodem_%(tid)s.tif' %{'tid':tile_id}
               
        geotiffpn = os.path.join(self.tempFP, geotiffn)
           
        if self.pp.process.parameters.flags:
            
            if self.pp.process.parameters.flags[0] == '-':
                
                flags = self.pp.process.parameters.flags
                
            else:
                
                flags = '-%s' %( self.pp.process.parameters.flags )      
        else:
            
            flags = ""
                                        
        paramsD = {'flags':flags, 'output': grasslayer, 'geotif':geotiffpn, 
                   'mod':self.pp.process.parameters.mod,
                   'size':self.pp.process.parameters.size,
                   'diffthreshold':self.pp.process.parameters.diffthreshold,
                   'tid':tile_id,
                   'areashpfpn':areashpfpn,
                   'ptshpfpn':ptshpfpn,
                   'cellarea': self.cellAreakm2,
                   'query': self.pp.process.parameters.query,
                   'burnattribute': 'filldem',
                   'dst': self.dstLayer.FPN,
                   'mem':self.pp.process.parameters.memory}
           
        cmd = '    r.hydrodem %(flags)s input=srcDEM output=%(output)s size=%(size)d mod=%(mod)d memory=%(mem)d --overwrite\n' %paramsD
        
        self.scriptF.write(cmd)
        
        self._InnerLoop(paramsD,rNt, rSt, rEt,rWt, False, False, True)
                  
    def _GrassDemReInterpolateTiles(self):
        '''
        '''

        srcComp = self.pp.srcCompL[0]
        
        dstComp = self.pp.dstCompL[0]
        
        if self.pp.process.parameters.asscript:
            
            today = mj_dt.Today()
            
            scriptFP = os.path.join('/Volumes',self.pp.dstPath.volume, self.pp.procsys.dstsystem, self.pp.dstCompD[dstComp].source,'tiles','demfillscript')
            
            if not os.path.exists(scriptFP):
                
                os.makedirs(scriptFP)
                
            scriptFN = 'grassDEM-resample_%(comp)s-%(today)s.sh' %{'comp':dstComp, 'today':today}
            
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
                        
                        exitstr = 'EXITING, expecting a virtual mosaic as input for GDALDEMTiles\n    %s' %(srcLayerFPN)
                            
                        print (exitstr)
                            
                        SNULLE
                          
                else:
                    
                    srcLayerFPN = srcLayer.FPN
                
                srcSpatialRef, srcMeta = ktgis.GetRasterMetaData(srcLayerFPN)

                if self.pp.process.parameters.ew_res == 0:
                    
                    ew_res = srcMeta.cellsize
                    
                else:
                    
                    ew_res = self.pp.process.parameters.ew_res
                    
                    
                if self.pp.process.parameters.ns_res == 0:
                    
                    ns_res = srcMeta.cellsize
                    
                else:
                    
                    ns_res = self.pp.process.parameters.ns_res
                    
                # g.extension extension=r.resamp.rst
                       
                paramsD = {'locus':locus,'src':srcLayerFPN,
                    'flags':self.pp.process.parameters.flags,     
                    'output': 'resampleDEM', 
                   'ew_res':ew_res,
                   'ns_res':ns_res,
                   'overlap': self.pp.process.parameters.overlap,
                   'tension': self.pp.process.parameters.tension,
                   'dst':self.dstLayer.FPN}

                cmd = '# create a new mapset\n'
        
                cmd += 'g.mapset -c %(locus)s\n\n' %paramsD
                
                cmd += '# import the source DEM\n'
                
                cmd += 'r.in.gdal input=%(src)s output=srcDEM  --overwrite\n\n' %paramsD    
                                  
                cmd += '# Set region after DEM\n'
                
                cmd += 'g.region raster=srcDEM\n\n'  
                              
                cmd += 'r.resamp.rst %(flags)s input=srcDEM ew_res=%(ew_res)f ns_res=%(ns_res)f elevation=%(output)s\
                  overlap=%(overlap)d tension=%(tension)f  --overwrite\n\n' %paramsD
                  
                cmd += 'r.out.gdal input=%(output)s format=GTiff output=%(dst)s --overwrite\n\n' %paramsD
         
                self.scriptF.write(cmd)
                  
        self.scriptF.close()
        
        print ('Script file:',scriptFPN)
        