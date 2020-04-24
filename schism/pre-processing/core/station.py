#!/usr/bin/env python3.7
import copy
from scipy.interpolate import griddata
import numpy as np
import os

class Station(object):

    def __init__(self,fileout,hgrid,outputs, logger=None):
        """ Constructor"""

        if logger:
            self.logger = logger
            self.logger.info("Processing atmospheric forcing")

        self.hgrid=copy.deepcopy(hgrid.mesh)
        self.hgrid.longitude=hgrid.longitude
        self.hgrid.latitude=hgrid.latitude

        self.outputs=outputs
        self.fileout=fileout


        self.X=[]
        self.Y=[]
        self.Z=[]
        self.name=[]


    def add_station(self,station,name=None):
        if 'node' in station:
           node=station['node']
           x=self.hgrid.x[node-1]
           y=self.hgrid.y[node-1]
           h=self.hgrid.values[node-1]

        elif 'xy' in station:
           x=station['xy'][0]
           y=station['xy'][1]
           h= griddata((self.hgrid.x,self.hgrid.y),self.hgrid.values,station['xy'], method='linear')[0]

        elif 'lonlat' in station:
           node= griddata((self.hgrid.longitude,self.hgrid.latitude),np.arange(0,len(self.hgrid.latitude)),station['lonlat'], method='nearest')[0]
           x=self.hgrid.x[node]
           y=self.hgrid.y[node]
           h=self.hgrid.values[node]


        for depth in station['depth']:
            self.X.append(x)
            self.Y.append(y)
            self.Z.append(-1*(h*(depth/100.)))
            if name:
                self.name.append(name)


    def write_station(self):
        f=open(self.fileout,'w')

        f.write(" %.f %.f %.f %.f %.f %.f %.f %.f %.f\n" % \
        (self.outputs['elev'],self.outputs['air pressure'],\
         self.outputs['windx'],self.outputs['windy'],\
         self.outputs['temp'],self.outputs['sal'],\
         self.outputs['u'],self.outputs['v'],self.outputs['w']))

        f.write("%.f\n" % len(self.X))
        name=''
        for i in range(0,len(self.X)):
            if len(self.name)>0:
               name=self.name[i]

            f.write("%.f %.2f %.2f %.2f "  % (i+1,self.X[i],self.Y[i],self.Z[i]))
            f.write(u"!")  
            f.write(" %s\n" % name)

        f.close()
        self.logger.info("	%s exported" % self.fileout)