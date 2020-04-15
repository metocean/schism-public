#!/usr/bin/env python2.7
from mako.template import Template
from mako.runtime import Context
from io import StringIO
from mako import exceptions
import ast

class ModelConfig(object):
    
    def __init__(self, hydro,t0,
        wwm=None,
        sed=None,
        eco=None,
        icm=None,
        cosine=None,
        fib=None,
        ice=None,
        sed2d=None,
        marsh=None,
        logger=None):
        '''Constructor'''

        if logger:
            self.logger = logger
            

        self.config={}
        self.userconfig={}
        self.userconfig['hydro'] = hydro
        self.config['wwm'] = wwm
        self.config['sed'] = sed
        self.config['eco'] = eco
        self.config['icm'] = icm
        self.config['cosine'] = cosine
        self.config['fib'] = fib
        self.config['ice'] = ice
        self.config['sed2d'] = sed2d
        self.config['marsh'] = marsh
        

        for module in self.userconfig:
            # order matters since some params are overwritten
            self.config[module]=self._get_default_params(module)
            self._replace_default_params(module)
            self._set_params_relations(module)


        self.config['hydro']['start_year']=t0.year
        self.config['hydro']['start_month']=t0.month
        self.config['hydro']['start_day']=t0.day
        self.config['hydro']['start_hour']=t0.hour




    def _define_logic(self):
        '''Defines the logic and interdependency between 
           model parameters (in param.in)'''
        pass

    def _check_logic(self):
        '''Checks if logic and interdependency between
           model parameters are correct (in param.in)'''
        pass
    
    def _set_logic(self):
        '''Translates the logic and interdependency between 
           model paramters into parctical values (in param.in)'''
        pass

    def _flow_logic(self):
        '''Streamline methods/modules to be used according to the 
           logic and interdependency of model paramters. This is 
           here just as a reminder - should probably live elsewhere'''
        pass

    def _get_default_params(self,module):
        ''' Docstring'''
        default={}
        if module is 'hydro':
            default['ibcc']=0
            default['ibtp']=1
            default['rnday']=30
            default['dt']=100
            default['ntracer_gen']=0
            default['ntracer_age']=0
            default['nspool']=36
            default['ihfskip']=864
            default['start_year']=2000
            default['start_month']=1
            default['start_day']=1
            default['start_hour']=0
            default['ihot']=0
            default['nramp']=1
            default['dramp']=1
            default['nrampbc']=0
            default['drampbc']=1
            default['indvel']=1
            default['ihorcon']=0
            default['hvis_coef0']=0.025
            default['ishapiro']=1
            default['shapiro0']=0.5
            default['niter_shap']=1
            default['thetai']=0.6
            default['imm']=0
            default['ibdef']=10
            default['sfea0']=45
            default['iunder_deep']=0
            default['h1_bcc']=50
            default['h2_bcc']=100
            default['hw_depth']=1e6
            default['hw_ratio']=0.5
            default['ihydraulics']=0
            default['if_source']=0
            default['ihdif']=0
            default['nchi']=0
            default['dzb_min']=0.5
            default['dzb_decay']=0
            default['hmin_man']=1
            default['ic_elev']=0
            default['nramp_elev']=0
            default['gen_wsett']=0
            default['ihhat']=1
            default['inunfl']=0
            default['h0']=0.01
            default['nadv']=1
            default['dtb_max']=30
            default['dtb_min']=10
            default['inter_mom']=0
            default['kr_co']=1
            default['itr_met']=3
            default['h_tvd']=5
            default['eps1_tvd_imp']=1e-4
            default['eps2_tvd_imp']=1e-14
            default['ip_weno']=2
            default['courant_weno']=0.5
            default['nquad']=2
            default['ntd_weno']=1
            default['epsilon1']=1e-3
            default['epsilon2']=1e-10
            default['i_prtnftl_weno']=0
            default['nws']=0
            default['iwind_form']=-1
            default['impose_net_flux']=0
            default['ihconsv']=0
            default['isconsv']=0
            default['itur']=3
            default['dfv0']=1e-2
            default['dfh0']=1e-4
            default['mid']="'KL'"
            default['stab']="'KC'"
            default['xlsc']=0.1
            default['step_nu_tr']=86400
            default['s1_mxnbt']=0.5
            default['s2_mxnbt']=3.5
            default['flag_fib']=1
            default['slr_rate']=120
            default['isav']=0
            default['sav_cd']=1.13
            default['nstep_ice']=1
            default['elev']=1
            default['bott']=0 #bottom stress
            default['wind']=0 #wind velocity vector [m/s] {wind_speed}
            default['wist']=0 #wind stress vector [m^2/s/s] {wind_stress}
            default['dav']=0 #depth-averaged vel vector [m/s] {dahv}
            default['vert']=0 #vertical velocity [m/s] {vertical_velocity}
            default['temp']=0 #water temperature [C] {temp}
            default['salt']=0 #water salinity [PSU] {salt}
            default['conc']=0 #water density [kg/m^3] {water_density}
            default['tdff']=0 #eddy diffusivity [m^2/s] {diffusivity}
            default['vdff']=0 #eddy viscosity [m^2/s] {viscosity}
            default['kine']=0 #turbulent kinetic energy {TKE}
            default['mixl']=0 #turbulent mixing length [m] {mixing_length}
            default['hvel']=0 #horizontal vel vector [m/s] {hvel}
            default['GEN_1']="'0'" #1st tracer {GEN_1}
            default['AGE_1']="'0'" #

            default['flag_ic_gen']=0 #EN (user defined module)
            default['flag_ic_age']=0 #Age
      
            default['flag_ic_timor']=0 #TIMOR
            default['flag_ic_fabm']=0 #FABM
            default['inu_tr_t']=0 #T
            default['inu_tr_s']=0 #S
            default['inu_tr_gen']=0 #GEN
            default['inu_tr_age']=0 #Age
            
            default['inu_tr_timor']=0 #TIMOR 
            default['inu_tr_fabm']=0 #FABM 
            default['icou_elfe_wwm']=0 #
            default['iwbl']=0 #
            default['hmin_radstress']=1 #
            default['turbinj']=0.15 #']=0 #
            default['shorewafo']=0 #

            default['sed_class']=0 #
            default['inu_tr_sed3d']=0 #SED3D
            default['flag_ic_sed3d']=0 #SED3D
            default['flag_ic_icm']=0 #ICM
            default['inu_tr_icm']=0 #ICM

            default['inu_tr_cosine']=0 #CoSINE
            default['flag_ic_cosine']=0 #CoSINE
            default['eco_class']=0 #
            default['inu_tr_ecosim']=0 #EcoSim 
            default['flag_ic_ecosim']=0 #EcoSim

            default['msc2']=24
            default['mdc2']=30
            default['flag_ic_fib']=0 #FI
            default['inu_tr_fib']=0


            default['nhot']=1
            default['nhot_write']=(24*3600)/default['dt']

            default['iout_sta']=0
            default['nspool_sta'] = (3600)/default['dt'] 


        if module is 'wwm':
            pass

        
        if module is 'sed':
            pass


        if module is 'eco':
            pass


        if module is 'icm':
            pass


        if module is 'cosine':
            pass

        if module is 'fib':
            pass
            

        if module is 'ice':
            pass

        if module is 'sed2d':
            pass

        if module is 'marsh':
            pass

 

        return default
    def _replace_default_params(self,module):
        for key in self.config[module]:
            if key in self.userconfig[module]:
                self.config[module][key]=self.userconfig[module][key]

        
    def _set_default_wwm(self):
        # WWM - waves params
        if 'coupling' in self.nest.physics.__dict__:
                wwm_params = parsectl(join(STATIC_DIR, 'wwm_default_input_params.yaml')).__dict__
                wwm_params = dict((k.lower(), v) for k, v in wwm_params.items())
        
                wwm_output = parsectl(join(STATIC_DIR, 'wwm_default_output_params.yaml')).__dict__
                wwm_output = dict((k.lower(), v) for k, v in wwm_output.items())
                
                config_dict['wwm_params'] = dict(container=Container(), defaults=wwm_params)
                config_dict['wwm_output'] = dict(container=Container(), defaults=wwm_output)
        
        for cfg in config_dict.keys():
            if cfg in self.nest:
                exec("params = self.nest.%s.__dict__"%cfg)
                params = dict((k.lower(), v) for k, v in params.items())
                for pkey in config_dict[cfg]['defaults']:
                    if pkey in params:
                        params[pkey] = params.get(pkey, config_dict[cfg]['defaults'][pkey])
                    else:
                        params[pkey] = config_dict[cfg]['defaults'][pkey]   
                config_dict[cfg]['container'].__dict__ = params
            else:   
                config_dict[cfg]['container'].__dict__ = config_dict[cfg]['defaults']
            exec("self.nest.%s = config_dict[cfg]['container']"%cfg) 
        
    def _set_params_relations(self,module):
        ''' Docstring'''

        self._add_more_tracers()
        self._update_timings()
        self._choose_diffusion()

    def _choose_diffusion(self):
        mode=self.userconfig['hydro'].get('mode','diffusion 1')
        inter_mom=self.userconfig['hydro'].get('inter_mom',0)
        if mode=='diffusion 1':
            self.config['hydro']['inter_mom']=0
            self.config['hydro']['ishapiro']=1
            self.config['hydro']['ihorcon']=0
            self.config['hydro']['indvel']=0
        elif mode=='diffusion 2':
            self.config['hydro']['inter_mom']=inter_mom
            self.config['hydro']['ishapiro']=0
            self.config['hydro']['ihorcon']=0
            self.config['hydro']['indvel']=1
        elif mode=='dispersion':
            self.config['hydro']['inter_mom']=0
            self.config['hydro']['ishapiro']=0
            self.config['hydro']['ihorcon']=0
            self.config['hydro']['indvel']=0
        else:
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print('!!!!!!!!Mode must be diffusion 1 or diffusion 2 or dispersion (eddying)!!!!!!!!!!!!!')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            import sys;sys.exit(-1)


    def _update_timings(self):
        dt=self.config['hydro']['dt']

        self.config['hydro']['nhot_write']=(self.userconfig['hydro']['hotstart dt']*3600)/dt
        self.config['hydro']['nspool_sta'] = (self.userconfig['hydro']['station dt']*3600)/dt
        self.config['hydro']['nspool_sta'] = (self.userconfig['hydro']['station dt']*3600)/dt
        self.config['hydro']['nspool'] = (self.userconfig['hydro']['output dt']*3600)/dt
        self.config['hydro']['ihfskip'] = (self.userconfig['hydro']['file length']*3600)/dt


    def _add_more_tracers(self):
        if self.config['hydro']['ntracer_gen']>0:
            self.config['hydro']['GEN_1']=''
            for n in range(0,self.config['hydro']['ntracer_gen']):
                self.config['hydro']['GEN_1']=self.config['hydro']['GEN_1']+'  iof_gen(%s) = 1 !N!' % str(n+1)
        else:
            self.config['hydro']['GEN_1']='!   iof_gen(1) = 0'

        self.config['hydro']['GEN_1']='"'+self.config['hydro']['GEN_1']+'"'


        if self.config['hydro']['ntracer_age']>0:
            self.config['hydro']['AGE_1']=''
            for n in range(0,self.config['hydro']['ntracer_age']):
                self.config['hydro']['AGE_1']=self.config['hydro']['AGE_1']+'  iof_age(%s) = 1 !N!' % str(n+1)
        else:
            self.config['hydro']['AGE_1']='!  iof_age(1) = 0 !N!!   iof_age(2) = 0'

        self.config['hydro']['AGE_1']='"'+self.config['hydro']['AGE_1']+'"'



    def make_config(self, cmdfile,paramfile,module):
        '''Docstring''' 

        self.logger.info("\tWriting:%s" %paramfile)
        cmdblank = Template(filename=cmdfile,input_encoding='utf-8',output_encoding='utf-8')
        buf = StringIO()
        stri='Context(buf,'
        
        for key in self.config[module]:
            stri=stri+key+'='+str(self.config[module][key])+','

        stri=stri[0:-1]+')'
        ctx = eval(stri.encode('utf-8'))
        

        #cmdblank.render_context(ctx)
        try:
            cmdblank.render_context(ctx)
        except:
            self.logger.info(exceptions.text_error_template().render())

	
        a=open(paramfile,'w')
        a.write(buf.getvalue())
        a.close()



