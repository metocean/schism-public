#!/usr/bin/env python3
# prepare schism input
import os,sys
from os.path import join
sys.path.append(join(os.path.dirname(__file__),'core')) 
import glob
import argparse
import yaml
import json
import logging
import dateutil
from matplotlib.dates import num2date,date2num
from hgrid import HorizontalGrid
from vgrid import VerticalGrid
from param import ModelConfig
from bctide import  BCinputs
from datasourcing import download_data
from openbnd import  OpenBoundaries
from meteo import Meteo
from ic import InitialConditions
from atmospheric import get_convergence
from station import Station

logging.basicConfig(filename=None,
                    filemode='w',
                    format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='[%Y-%m-%d %H:%M:%S]',
                    level=10)

class SCHISM():
    ''' 
    This is the wrapper for preparing and setting up SCHISM.
    '''

    def __init__(self, rootdir, hydro_config,vgrid_config,obc,param_tmp,exec_bin,
                 hgrid_file,timing,input_files,forcings,ic=None,meteo=None,stations=None,
                 indir=None,
                 logdir=None,
                 errors=None,
                 epsg=2193,
                 **kwargs):

        # ----------------------------------------------------- Mandatory arguments -----------
        self.rootdir = rootdir
        self.hydro_config  = hydro_config
        self.vgrid_config  = vgrid_config
        self.input_files = input_files
        self.forcings= forcings
        self.meteo=meteo
        self.obc  = obc
        self.timing= timing
        self.ic=ic
        self.stations=stations

        # ----------------------------------------------------- Run paramterss -----------
        self.outdir = join(self.rootdir, 'outputs')
        self.logdir = logdir or join(self.rootdir, 'log')
        self.indir = indir or join(self.rootdir, 'in')

        #------------------------------------------------------- Model environment ---------
        self.param_tmp = param_tmp
        self.exe_dirs = exec_bin
        self.hfile= hgrid_file
        self.epsg=epsg



        #----------------------------------------------------------------- Logging ---------
        self.logger = logging
        self.logger.level = 20
        self.errorsfname = errors

    def _set_environment(self):
        # creates implementation logs directory
        if not os.path.isdir(self.logdir):
            os.makedirs(self.logdir)

        # creates hotfile directory
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)

        # Open error files
        errorsfname = self.errorsfname or join(self.logdir, 'errors.log')
        self.errors = open(errorsfname, 'a')

    def prepare_run(self):
        self.logger.info('----------------------------------------------------------')
        self.logger.info('\tRunning: %s' % (self.exe_dirs))
        self.logger.info('\ttemplate: %s' % (self.param_tmp))
        self.logger.info('\trootdir: %s' % (self.rootdir))
        self.logger.info('----------------------------------------------------------')
        self._set_environment()

        #----------------------- Load horizontal grid objects ----------
        self.logger.info('\tReading: %s' % (self.hfile))
        self.hgrid = HorizontalGrid(self.hfile,format='gr3',epsg=self.epsg, logger=self.logger)
        self.hgrid.copy_to_root(self.rootdir)



        #----------------------- Load vertical grid objects ----------
        vgrid_reader = VerticalGrid(logger=self.logger)
        if 'vfile' not in self.vgrid_config:
            self.vgrid_config['vfile']=join(self.rootdir,'vgrid.in')

        if not os.path.isfile(self.vgrid_config['vfile']):
            self.logger.info('\tCreating: %s' % (self.vgrid_config['vfile']))
            vgrid_reader.write(**self.vgrid_config)
   
        self.vgrid=vgrid_reader.load(self.vgrid_config['vfile'])
        


        lat0=sum(self.hgrid.latitude)/len(self.hgrid.latitude)
        t0 = dateutil.parser.parse(self.timing["time start"])
        t1 = dateutil.parser.parse(self.timing["time end"])

        #  #----------------------- Download Initial and Boundary fields ---------- 
        dwnl=download_data(t0,t1,logger=self.logger)
        for file in self.input_files.keys():     
          dwnl.get_input_data(self.input_files[file])


        self.logger.info('----------------------------------------------------------')
        #----------------------- Write command file (param.in) ------------------
        cfg = ModelConfig(hydro=self.hydro_config,t0=t0,logger=self.logger)
        cfg.make_config(self.param_tmp,join(self.rootdir,'param.nml'),'hydro')

        # #----------------------- Set Boundary Conditions (bctides.in) -----------
        # # Store boundary arrays in each obc bctype object (Ex: self.obc['btype']['7']['iettype'])

        bcinput = BCinputs(obc=self.obc,nnode=self.hgrid.nnode, lat0=lat0,t0=t0, logger=self.logger)
        bcinput.make_bctides(join(self.rootdir,'bctides.in'))

       #  # ------------------- Create Ocean boundary forcing -----------------
        for key in self.forcings.keys():
          if not os.path.isfile(join(self.rootdir,key)):
            Obf = OpenBoundaries(obc=self.forcings[key],hgrid=self.hgrid,vgrid=self.vgrid,t0=t0,t1=t1, logger=self.logger)
            if 'tidal' in self.forcings[key]:
              Obf.add_tide(self.forcings[key]['tidal'])

            if 'residual' in self.forcings[key]:
              Obf.add_res(self.forcings[key]['residual'])

            Obf.make_boundary(join(self.rootdir,key),self.forcings[key].get('dt',3600))


       #  # ------------------- Create Oceanic initial conditions ---------------------
        if self.ic:
           for key in self.ic:
               filename=os.path.join(self.rootdir,key)
               if not os.path.isfile(filename):
                 ic = InitialConditions(filename,hgrid=self.hgrid,t0=t0,value=self.ic[key].get('value',None),\
                  ncfile=self.ic[key].get('filename',None),\
                  var=self.ic[key].get('var',None),\
                  shapefile=self.ic[key].get('shapefile',None),logger=self.logger)



       # #------------------------- Check/Prepare for hotstart --------------------
       #  self.hot = HotStart(nest=self.nest, logger=self.logger)
       #  self.hot.set_hotstart()


        # ------------------- Create Atmospheric forcing --------------------
        if self.meteo and not os.path.isfile(os.path.join(self.rootdir,'sflux','sflux_inputs.txt')):

          atm=Meteo(os.path.join(self.rootdir,'sflux'),t0=t0,t1=t1,dt=self.meteo.pop('dt'),logger=self.logger)
          for key in range(1,len(self.meteo.keys())+1):
            filename=self.meteo[key].pop('filename')
            atm.add_dataset(filename,self.meteo[key])

          atm.make_meteo()

        if self.meteo and not os.path.isfile(os.path.join(self.rootdir,'windrot_geo2proj.gr3')):
        ## Create the GR3 files
          filename = os.path.join(self.rootdir,'windrot_geo2proj.gr3')
          ic = InitialConditions(filename,hgrid=self.hgrid,t0=t0,value=get_convergence(self.hgrid.latitude,self.hgrid.longitude,self.epsg),\
            logger=self.logger)


        if self.stations and not os.path.isfile(os.path.join(self.rootdir,'station.in')):
          outputs=self.stations.pop('outputs')
          st = Station(os.path.join(self.rootdir,'station.in'),self.hgrid,outputs, logger=self.logger)
          for key in self.stations:
            st.add_station(self.stations[key],name=key)

          st.write_station()

       #  # ------------------- Create Wave boundary forcing -----------------
       #  if self.nest.wave:
       #      wave = Waves(nest=self.nest, logger=self.logger)
       #      wave.make_waves()
        

def load_action(yfile):
    with open(yfile, 'r') as f:
        return yaml.load(f)

def cycle_sim(action='tests/model.tinyapp.yml', **kwargs):
    """ Simulate scheduler call 
    """
    config = load_action(action)
    config.update(kwargs)
    model = SCHISM(**config)
    model.prepare_run()

#    model.set_mpiexec(ncores, hosts='localhost')
#    model.run()
#    return model

def run_as_script():
    parser = argparse.ArgumentParser(description="SCHISM Wrapper")
    parser.add_argument('-a', '--action', type=str,
                        default='tests/model.tinyapp.yml',
                        help='Model action')
    parser.add_argument('-k', '--kwargs', type=json.loads, default={},
                        help="""kwargs dictionary - e.g. '{"spinup":0}'""")
    args = parser.parse_args()
    
    
    cycle_sim(action=args.action, **args.kwargs)

if __name__ == '__main__':
    run_as_script()
