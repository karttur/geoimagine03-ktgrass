'''
Created on 2 Apr 2021

@author: thomasgumbricht

Core copied from https://grass.osgeo.org/grass78/manuals/libpython/script.html#module-script.setup
'''


import os
import sys
import subprocess



import geoimagine.grass.script as grass

import geoimagine.support.karttur_dt as mj_dt

class ProcessGRASS():
    '''class for binding to GRASS commands''' 
      
    def __init__(self, pp, session):
        '''
        '''
        
        self.session = session
                
        self.pp = pp  
        
        self.verbose = self.pp.process.verbose 
        
        self.session._SetVerbosity(self.verbose) 
        
        
                
        # Direct to subprocess
        if self.pp.process.processid == 'Grass1to1Tiles':
            
            self._Grass1to1Tiles()
            
        elif self.pp.process.processid == 'GrassOnetoManyTiles':
            
            self._GrassOnetoManyTiles()
                        
        else:
        
            print (self.pp.process.processid)

            exit('Process not recognized in processGRASS')

    def _Grass1to1Tiles(self):
        '''
        '''

        srcComp = self.pp.srcCompL[0]
        
        dstComp = self.pp.dstCompL[0]
        
        if self.pp.process.parameters.asscript:
            
            today = mj_dt.Today()
            
            scriptFP = os.path.join('/Volumes',self.pp.dstPath.volume, self.pp.procsys.dstsystem, self.pp.dstCompD[dstComp].source,'region','reprojectscript')
            
            if not os.path.exists(scriptFP):
                
                os.makedirs(scriptFP)
                
            scriptFN = 'grass_%(comp)s-%(today)s.sh' %{'comp':dstComp, 'today':today}
            
            scriptFPN = os.path.join(scriptFP,scriptFN)
            
            self.scriptF = open(scriptFPN,'w')
            
            writeln = '# Script created by Kartturs GeoImagine Framework for GRASS processing %s, created %s\n' %(dstComp, today)
    
            self.scriptF.write(writeln)
            
        for datum in self.pp.srcPeriod.datumL:
                           
            for locus in self.pp.srcLayerD:
                                    
                srcLayer = self.pp.srcLayerD[locus][datum][srcComp]
                
                if not os.path.exists(srcLayer.FPN):
                    
                    infostr = '    SKIPPING - input layer missing for ProcessGRASS\n        %s' %(srcLayer.FPN)
        
                    print (infostr)
                    
                    continue
                
                dstLayer = self.pp.dstLayerD[locus][datum][dstComp]
                
                if dstLayer._Exists() and not self.pp.process.overwrite:
                
                    continue
                
                if self.pp.process.parameters.mosaic:
                    
                    srcLayerFPN = srcLayer.FPN.replace(srcLayer.comp.ext, '-full.vrt')
                    
                    if not os.path.exists(srcLayerFPN):
                        
                        exitstr = 'EXITING, expecting a virtual mosaic as input for GDALDEMTiles'
                        
                        exit(exitstr)
                          
                else:
                    
                    srcLayerFPN = srcLayer.FPN
                    
                cmdReplaceD = {'srcFPN':srcLayerFPN,
                                'srcFN':os.path.splitext(self.pp.srcLayerD[locus][datum][srcComp].FN)[0].replace('-',''),
                                'dstFPN':self.pp.dstLayerD[locus][datum][dstComp].FPN}
                                
                                
                # set mapset to region
                qwargs = {'flags':'c','mapset':locus}
                
                if not self.pp.process.dryrun:
                            
                    if self.pp.process.parameters.asscript:
                        
                        self._WriteGrassCmd('g.mapset',qwargs)
                    
                    else:
                        
                        grass.run_command('g.mapset', **qwargs)
                           
                for param in self.pp.process.parameters.subparameter:
                                       
                    cmdD = dict ( list ( param.__dict__.items() ) )
                     
                    for cmd in cmdD:
                        
                        kwargs = dict ( list ( cmdD[cmd].__dict__.items() ) )
                        
                        # Always set overwrite to True
                        
                        kwargs['overwrite'] = True
                                      
                        for key in kwargs:
                            
                            if kwargs[key] in cmdReplaceD:
                                
                                kwargs[key] = cmdReplaceD[kwargs[key]]
                                
                                
                        if cmd in ['r.out.gdal'] and self.pp.process.parameters.mosaic: 
                            
                            self._SetOutRegion()
                                
                        if self.verbose:   
                            
                            print ('        ',cmd,kwargs)
                            
                        if not self.pp.process.dryrun:
                            
                            if self.pp.process.parameters.asscript:
                                
                                self._WriteGrassCmd(cmd,kwargs)
                            
                            else:
                                           
                                grass.run_command(cmd, **kwargs) 
                        
                        if cmd in ['r.out.gdal'] and self.pp.process.parameters.mosaic: 
                            
                            self._ReSetRegion(self.pp.process.parameters.regionlayer)
                        
                        '''        
                        if cmd == 'r.out.gdal' and self.pp.process.parameters.mosaic: 
                            
                            # set the regions to the original src file (i.e. the core tile)
                            qwargs = {'flags':'ec','input':self.pp.srcLayerD[locus][datum][srcComp].FPN,'output':'dummy'}
                            
                            grass.run_command('r.out.gdal', **qwargs)
                            
                        print ('kwargs',kwargs)
                        
                        grass.run_command(cmd, **kwargs) 
                        '''
        if self.pp.process.parameters.asscript:
            
            self.scriptF.close()
            
            infostr = 'To render the commands, run the script file:\n    %s' %(scriptFPN)
            
            print (infostr)
                        
    def _GrassOnetoManyTiles(self):
        '''
        '''

        srcComp = self.pp.srcCompL[0]
        
        cmpstr = '%s-etc' %(self.pp.dstCompL[0])
        
        if self.pp.process.parameters.asscript:
            
            today = mj_dt.Today()
            
            scriptFP = os.path.join('/Volumes',self.pp.dstPath.volume, self.pp.procsys.dstsystem, self.pp.dstCompD[self.pp.dstCompL[0]].source,'region','grasscript')
            
            if not os.path.exists(scriptFP):
                
                os.makedirs(scriptFP)
                
            scriptFN = 'grass_%(comp)s-%(today)s.sh' %{'comp':cmpstr, 'today':today}
            
            scriptFPN = os.path.join(scriptFP,scriptFN)
            
            self.scriptF = open(scriptFPN,'w')
            
            writeln = '# Script created by Kartturs GeoImagine Framework for GRASS processing %s, created %s\n' %(cmpstr, today)
    
            self.scriptF.write(writeln)
            
        for datum in self.pp.srcPeriod.datumL:
                           
            for locus in self.pp.srcLayerD:
                                    
                srcLayer = self.pp.srcLayerD[locus][datum][srcComp]
                
                if not os.path.exists(srcLayer.FPN):
                    
                    infostr = '    SKIPPING - input layer missing for ProcessGRASS\n        %s' %(srcLayer.FPN)
        
                    print (infostr)
                    
                    continue
                
                dstLayerD = {}
                
                skipLayerCmdL = []
                
                for dstComp in self.pp.dstLayerD[locus][datum]:
                    
                    dstLayer = self.pp.dstLayerD[locus][datum][dstComp]
                
                    if dstLayer._Exists() and not self.pp.process.overwrite:
                        
                        skipLayerCmdL.append(dstComp)
                        
                        skipLayerCmdL.append(dstComp.replace('out',''))
                        
                        if self.verbose:
                            
                            infostr = '        Output layer already exists:\n        %s ' %(dstLayer.FPN)
                            
                            print (infostr)
                    
                        continue
                    
                    if self.verbose:
                            
                        infostr = '        Adding output composition (%s) as layer:\n        %s ' %(dstComp,dstLayer.FPN)
                            
                        print (infostr)   
                    
                    dstLayerD[dstComp] = self.pp.dstLayerD[locus][datum][dstComp]
                    
                if len(dstLayerD) == 0:
                    
                    if self.verbose:
                            
                        infostr = '        SKIPPING - all output layers exists: %s' %(locus)
                            
                        print (infostr)
                    
                    continue
                                    
                if self.pp.process.parameters.mosaic:
                    
                    srcLayerFPN = srcLayer.FPN.replace(srcLayer.comp.ext, '-full.vrt')
                    
                    if not os.path.exists(srcLayerFPN):
                        
                        exitstr = 'EXITING, expecting a virtual mosaic as input for GDALDEMTiles'
                        
                        exit(exitstr)
                          
                else:
                    
                    srcLayerFPN = srcLayer.FPN
                    
                cmdReplaceD = {'srcFPN':srcLayerFPN,
                                'srcFN':os.path.splitext(self.pp.srcLayerD[locus][datum][srcComp].FN)[0].replace('-','')}
                
                for cl in dstLayerD:
                                         
                    cmdReplaceD[cl] = dstLayerD[cl].FPN
                                        
                # set mapset to region
                qwargs = {'flags':'c','mapset':locus}
                
                if not self.pp.process.dryrun:
                            
                    if self.pp.process.parameters.asscript:
                        
                        self._WriteGrassCmd('g.mapset',qwargs)
                    
                    else:
                        
                        grass.run_command('g.mapset', **qwargs)
                    
                if self.pp.process.parameters.mosaic:
                    
                    # Import the original tile to use as export termpalte
                    kwargs = {'input':self.pp.srcLayerD[locus][datum][srcComp].FPN,'output':'originalTile','overwrite':True}
                    
                    if self.verbose:   
                            
                        print ('        ','r.in.gdal',kwargs)
                        
                    if not self.pp.process.dryrun:
                        
                        if self.pp.process.parameters.asscript:
                            
                            self._WriteGrassCmd('r.in.gdal',kwargs)
                        
                        else:
                    
                            grass.run_command('r.in.gdal', **kwargs)
                                
                srcLayer = False 
                           
                for param in self.pp.process.parameters.subparameter:
                    
                    cmdD = dict ( list ( param.__dict__.items() ) )
                                        
                    for cmd in cmdD:
                        
                        kwargs = dict ( list ( cmdD[cmd].__dict__.items() ) )
                        
                        # if this process includes an input or output that is already done, skip
                        cont = False 
                        
                        for key in kwargs:
                            
                            if kwargs[key] in skipLayerCmdL:
                                
                                cont = True
                                
                        if cont:
                            
                            continue

                        # Always set overwrite to True
                        
                        if not cmd == 'null':
                        
                            kwargs['overwrite'] = True
                                                              
                        for key in kwargs:
                            
                            if kwargs[key] in cmdReplaceD:
                                
                                kwargs[key] = cmdReplaceD[kwargs[key]]
                                
                                if cmd in ["v.out.ogr","v.in.ogr"]:
                                    
                                    ext = os.path.splitext( kwargs['output'] )[1]
                                    
                                    print (kwargs['output']  )
                                    
                                    kwargs["output"] = kwargs["output"].replace(ext,'.shp')
                                           
                        if cmd == 'r.in.gdal'  and not srcLayer:
                            
                            srcLayer = kwargs['output']
                                
                        if cmd in ['r.out.gdal'] and self.pp.process.parameters.mosaic: 
                            
                            self._SetOutRegion()
                                
                        if self.verbose:   
                            
                            print ('        ',cmd,kwargs)
                            
                        if not self.pp.process.dryrun:
                            
                            if self.pp.process.parameters.asscript:
                                
                                self._WriteGrassCmd(cmd,kwargs)
                            
                            else:
                                           
                                grass.run_command(cmd, **kwargs) 
                        
                        if cmd in ['r.out.gdal'] and self.pp.process.parameters.mosaic: 
                            
                            self._ReSetRegion(self.pp.process.parameters.regionlayer)

                '''            
                if self.pp.process.parameters.mosaic:
                    
                    # Remove the original tile
                    kwargs = {'flags':'f', 'type':'raster', 'name': 'originalTile'}
                    
                    if self.verbose:   
                            
                        print ('        ','g.remove',kwargs)
                    
                    if not self.pp.process.dryrun:
                        
                        if self.pp.process.parameters.asscript:
                            
                            self._WriteGrassCmd('g.remove',kwargs)
                        
                        else:
                    
                            grass.run_command('g.remove', **kwargs)
                '''
                                              
        if self.pp.process.parameters.asscript:
            
            self.scriptF.close()
            
            infostr = 'To render the commands, run the script file:\n    %s' %(scriptFPN)
            
            print (infostr)
            
            
    def _SetOutRegion(self):
        '''
        '''
        # set the regions to the original src file (i.e. the core tile)
        qwargs = {'raster':'originalTile'}
        
        if self.verbose:   
        
                print ('        ','g.region',qwargs)
        
        if not self.pp.process.dryrun:
            
            if self.pp.process.parameters.asscript:
                
                self._WriteGrassCmd('g.region',qwargs)
            
            else:
                
                grass.run_command('g.region', **qwargs)
                
    def _ReSetRegion(self,regionlayer):
        '''
        '''
                            
        # reset the regions to the virtual, enlarged, region)
        qwargs = {'raster':regionlayer}
        
        if self.verbose:   
        
            print ('        ','g.region',qwargs)
        
        if not self.pp.process.dryrun:
            
            if self.pp.process.parameters.asscript:
            
                self._WriteGrassCmd('g.region',qwargs)
            
            else:
            
                grass.run_command('g.region', **qwargs)
                        
    def _WriteGrassCmd(self,cmd,kwargs):
        
        if cmd == 'null':
            
            cmd = ''
                
        for key, value in kwargs.items():
                        
            if key == 'flags':
                
                cmd += ' -%s' %(value)
            
            elif isinstance(value, bool) and value: 
            
                cmd += ' --%s' %(key)
                
            elif key == 'null':
                
                cmd += '%s' %(value)
                
            else:
                
                cmd += ' %s=%s' %(key,value)
                
        cmd += '\n'
        
        self.scriptF.write(cmd)
            
        
                    