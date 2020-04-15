from matplotlib.dates import date2num
import numpy as np
from ttide.t_getconsts import t_getconsts
from ttide.t_vuf import t_vuf


TRACERS_MODULE=['GEN','AGE','SED3D','EcoSim','ICM','CoSiNE','FIB','TIMOR']

def Calculate(lat, t0, cons):

    const = t_getconsts(np.array([]))
    Const= [con.decode('UTF-8') for con in const[0]['name']] 


    consindex = [Const.index(con.ljust(4)) for con in cons]

    # V: astronomical phase, U: nodal phase modulation, F: nodal amplitude correction
    v,u,f = t_vuf('nodal', np.array(t0), consindex, lat)
    tear = 360.*(v+u)
    tfreq = (2*np.pi)*const[0]['freq'][consindex]/3600.
    talpha = const[0]['name'][consindex]
    return talpha, tfreq, tear, f

def develop_bnd(obc):
    # check if we have a loop in the dictionnary
    if 'cons' in obc:
        obc.pop('cons')

    btype=dict()
    
    for nk in obc.keys():
        if type(nk) is not int:
            [aa,bb]=nk.split('-')
            for nkk in range(int(aa)-1,int(bb)):
                btype[int(nkk)]=obc[nk]
                
        else:
    
            btype[int(nk)-1]=obc[nk]

    return btype

class BCinputs(object):
    """
    Class that manages input file generation for a model simulation.
    """

    def __init__(self,obc,nnode,lat0,t0, logger=None):
        '''Docstring'''  

        if logger:
            self.logger = logger

        self.lat0 = lat0
        self.t0= t0
        self.t0dec=date2num(t0)
        self.obc=obc
        self.nnode=nnode



    
    def _write_iettype(self,btypes,iettype):
            

            ## TAKE CARE OF iettype
            if iettype==2: #this boundary is forced by a constant elevation
                ethconst =btypes['iettype']['const']
                self.bctides.write("%.2f\n" % (ethconst)) 

            elif iettype==3 or iettype==5: #forced in frequency domain
                print('This option is not implemented yet')
                pass

    def _write_ifltype(self,btypes,ifltype):

            ## TAKE CARE OF ifltype
        if ifltype==2: #forced in frequency domain
            const=btypes['ifltype']['const']
            self.bctides.write("%.4f\n" % (const))
        elif ifltype==3 or ifltype==5: #this boundary is forced by a constant elevation
            print("not implemented yet")
        elif ifltype==-4: #time history of velocity 
            inflow = btypes['ifltype']['inflow']
            outflow = btypes['ifltype']['outflow']
            self.bctides.write("%.2f %.2f\n" % (inflow,outflow))
        elif ifltype==-1: #Flanther type radiation b.c.
            eta_m0 = btypes['ifltype']['eta_m0']
            qthcon = btypes['ifltype']['qthcon']
            self.bctides.write("%.2f %.2f\n" % (eta_m0,qthcon))


    def _write_tracers(self,btypes,module,flag):

        if flag==1 or flag==3 or flag==4: #
            tobc = btypes[module]['tobc']
            self.bctides.write("%.2f\n" % (tobc))
        elif flag==2: 
            const = btypes[module]['const']
            tobc = btypes[module]['tobc']
            self.bctides.write("%.2f\n" % (const))
            self.bctides.write("%.2f\n" % (tobc))

    def _write_bctides(self,filename, talpha=[], tfreq=[], tear=[], tnf=[]):
        '''
           ----------------------------------- Elevation b.c. section ------
            iettype(j): (DEFAULT: 5)

                0: elevations are not specified for this boundary (in this case the velocity must be specified).

                1: time history of elevation on this boundary
                  > no input in bctides.in; time history of elevation is read in from elev.th (ASCII);

                2: forced by a constant elevation
                  > ethconst (constant elevation value for this segment)

                3: forced by tides
                  > for k=1, nbfr
                  >   alpha(k) !tidal constituent name
                  >   for i=1, nond(j) !loop over all open boundary nodes on this segment
                  >       emo((j,i,k) efa (j,i,k) !amplitude and phase for each node on this open boundary
                  >   end for i
                  > end for k

                4: space- and time-varying input
                  > no input in this file; time history of elevation is read in from elev2D.th(binary)

                5: combination of (3) and (4):
                   > time history of elevation is read in from elev2D.th, and then added to tidal B.C. specified below
                   > for k = 1, nbfr
                   >     alpha(k) !tidal constituent name
                   >     for i = 1, nond(j) !loop over all open boundary nodes on this segment
                   >         emo((j, i, k) efa(j, i, k) !amplitude and phase for each node on this open boundary
                   >     end for i
                   > end for k

            ----------------------------------- Velocity b.c. section -----------------------
            ifltype(j): (DEFAULT: 5)

                0: velocity not specified (no input needed)

                1: time history of discharge on this boundary
                  > no input in this file; time history of discharge is read in from flux.th (ASCII)

                2: this boundary is forced by a constant discharge
                  > vthconst (constant discharge; note that a negative number means inflow)

                3: vel. (not discharge!) is forced in frequency domain
                  > for k=1, nbfr
                  >    vmo(j,k) vfa(j,k) (uniform amplitude and phase along each boundary segment)
                  >    end for;

                4 or -4: 3D input
                        > time history of velocity (not discharge!) is read in from uv3D.th (binary)

                5: combination of (4) and tides
                  > time history of velocity (not discharge!) is read in from uv3D.th (binary) and then added to tidal velocity specified below
                  >    for k=1, nbfr
                  >       alpha(k) (tidal constituent name)
                  >       for i=1, nond(j) (loop over all open boundary nodes on this segment)
                  >           umo(j,i,k) ufa(j,i,k) vmo(j,i,k) vfa(j,i,k) !amplitude and phase for (u,v) at each node on this open boundary
                  >       end for i
                  >    end for k

                -1: Flanther type radiation b.c. (iettype must be 0 in this case)
                   > eta_mean (mean elevation below)
                   > for i=1,nond(j) (loop over all nodes)
                   >     eta_m0(i) !mean elev at each node
                   > end for i
                   > vn_mean (mean normal velocity)
                   > for i=1,nond(j)
                   >    qthcon(1:Nz,i,j) (mean normal velocity at the node; at all levels)
                   > end for i

            ----------------------------- Temperature b.c. section ---------------------------
            itetype: (DEFAULT: ?)

                0: temperature not specified (no input needed)

                1: time history of temperature on this boundary
                  > tobc !nudging factor (between 0 and 1 with 1 being strongest nudging) for inflow; time history of temperature will be read in from TEM_1.th (ASCII)

                2: this boundary is forced by a constant temperature
                  > tthconst (constant temperature on this segment)
                  > tobc (nudging factor (between 0 and 1) for inflow)

                3: initial temperature profile for inflow
                  > tobc (nudging factor (between 0 and 1) for inflow)

                4: 3D input
                  > tobc (nudging factor (between 0 and 1); time history of temperature is read in from TEM_3D.th (binary)

            ---------------------------- Salintiy B.C. section ------------------------------
            isatype: (DEFAULT: ?)

                the structure is similar to temperature

            ---------------------------- Tracers B.C. section -------------------------------
            If any tracer module is invoked, you also need the corresponding B.C. part for each tracer module

            isatype: (DEFAULT: 0)

                the structure is similar to temperature and salinity

            end for !j: open boundary segment
            *** Notes: the tidal amplitudes and phases can be generated using utility scripts shown on the web.'''


        self.bctides = open(filename, 'w')
        # 48-character start time info string (only used for visualization with xmvis6)
        self.bctides.write(self.t0.strftime('%d/%m/%Y %H:%M:%S\n'))
        # ntip tip_dp > # of constituents used in earth tidal potential; cut-off depth for applying tidal potential (i.e., it is not calculated when depth < tip_dp).
        self.bctides.write("0 1000\n")
        # nbfr > # of tidal boundary forcing frequencies
        self.bctides.write("%d\n"%len(talpha))
        for k in range(0, len(talpha)):
            # tidal constituent name
        	self.bctides.write("%s\n"%(talpha[k])) 
            # angular frequency (rad/s), nodal factor, earth equilibrium argument (deg)
        	self.bctides.write("%.6f %.6f %.6f\n" % (tfreq[k], tnf[k], np.mod(tear[k],360)))

        


        #number of open boundary segments
        nopen=len(self.nnode)
        self.bctides.write('%s\n'%nopen)


        btypes=develop_bnd(self.obc)



        
        for k in range(0,nopen):
    # open boundary option
            iettype = btypes[k]['iettype']['value']
            ifltype = btypes[k]['ifltype']['value']
            itetype = btypes[k]['itetype']['value']
            isatype = btypes[k]['isatype']['value']

            
            self.bctides.write("%.f %.f %.f %.f %.f" % (len(self.nnode[k]),iettype,ifltype,itetype,isatype))
            ## add the tracers
            for modules in TRACERS_MODULE:
                if modules in btypes[k]:
                    self.bctides.write(" %.f" % (btypes[k][modules]['value']))
                else:
                    self.bctides.write(" %.f" % (0))
            
            self.bctides.write("\n" )


            self._write_iettype(btypes[k],iettype)
            self._write_ifltype(btypes[k],ifltype)


            self._write_tracers(btypes[k],'itetype',itetype) # temperature
            self._write_tracers(btypes[k],'isatype',isatype) # salinity

            for modules in TRACERS_MODULE:
                if modules in btypes[k]:
                    self._write_tracers(btypes[k],modules,btypes[k][modules]['value']) # salinity
                            
                        
        self.bctides.flush()
        self.bctides.close()


    def make_bctides(self,filename):
        '''Docstring'''
        
        if self.logger:
            self.logger.info("  Writing %s" %filename)


        talpha, tfreq, tear, tnf = Calculate(self.lat0, self.t0dec, self.obc['cons'].split(' '))
        self._write_bctides(filename,talpha, tfreq, tear, tnf)

