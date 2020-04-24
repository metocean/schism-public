import sys,os
from base_io import BaseIO
from pyschism.mesh import Hgrid



import argparse
# -------------------------------------------- BASE MODEL CLASSES -----------------------------
WGS84 = 4326 # Spatial Reference System
#-------------------------------------------- GRID METHODS -----------------------------
def transform_proj(x, y, in_espg=2193, out_espg=4326):
    ''' Performs cartographic transformations (converts from
        longitude,latitude to native map projection x,y coordinates and
        vice versa) using proj (http://trac.osgeo.org/proj/).
    '''
    in_proj = Proj(init='epsg:%s'%in_espg)
    out_proj = Proj(init='epsg:%s'%out_espg)
    x2, y2 = transform(in_proj, out_proj, x, y)
    return x2, y2

class HorizontalGrid(BaseIO):
    """A class that manages I/O of HGRID files
    - Transformation of grid files
    - Make available all grid data & metadata
    """

    def __init__(self,path_gr3,format='gr3',epsg=2193, logger=None):
        """ Constructor"""
        super(HorizontalGrid, self).__init__()

        if logger:
            self.logger = logger

        self.path_gr3=path_gr3
            
        # ----------------------------------------------------------- grid properties -----------       
        self.hgrid = self.load(epsg)
        self.mesh = self.hgrid

        self.nopen = len(self.mesh.boundaries[None]) # number of open boundaries segments    
        self.nnode = [self.mesh.boundaries[None][bnd]['indexes'] for bnd in self.mesh.boundaries[None]]
        self.h = self.mesh.values
        self.epsg=epsg

        if epsg == WGS84:
            self.longitude = self.mesh.x
            self.latitude = self.mesh.y
        else:
            self.path_ll=self.path_gr3[:-3]+'ll'
            self.longitude = self.mesh.get_x(crs="EPSG:%i"%WGS84)
            self.latitude = self.mesh.get_y(crs="EPSG:%i"%WGS84)

            if not os.path.isfile(self.path_ll):
                self.mesh.transform_to(WGS84)
                self.logger.info('\tCreating: %s' % (self.path_gr3[:-3]+'ll'))
                self.mesh.write(self.path_gr3[:-3]+'ll')
                self.mesh.transform_to(epsg)


 

 
    def __repr__(self):
        rep = "\n%s\n" % self.__class__
        rep += "\ninfo: %s\n\t" % self.info.__dict__
        rep += "\nGrid types:\n \t path_ll: '%s' \n\t path_gr3: '%s'\n" % (self.path_ll, self.path_gr3)
        return rep

    def get_properties(self):
        '''Docstring'''
        pass

    def load(self,epsg):
        '''Docstring'''
        return Hgrid.open(self.path_gr3,crs="EPSG:%i" % epsg)

    def close(self):
        '''Docstring'''
        self.hgrid.close()

    def copy_to_root(self,rootdir):
        if not os.path.isfile(os.path.join(rootdir,'hgrid.gr3') ):
            os.system('cp %s %s' % (self.path_gr3, os.path.join(rootdir,'hgrid.gr3') ))

        if self.epsg != WGS84 and not os.path.isfile(os.path.join(rootdir,'hgrid.ll') ):    
            os.system('cp %s %s' % (self.path_ll, os.path.join(rootdir,'hgrid.ll') ))

def run_as_script():
    parser = argparse.ArgumentParser(description="SCHISM grid reader")
    parser.add_argument('-g', '--hgrid', type=str,
                        help='grid file')
    parser.add_argument('-f', '--format', type=str,default='gr3',
                        help='grid format')
    parser.add_argument('-e', '--epsg', type=int, default=2193,
                        help="grid coordinate reference")

    args = parser.parse_args()
    
    
    gr=HorizontalGrid(args.hgrid,format=args.format,epsg=args.epsg,logger=None)

if __name__ == '__main__':
    run_as_script()
