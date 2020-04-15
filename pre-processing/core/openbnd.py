from bctide import develop_bnd
import numpy as np
from matplotlib.dates import num2date,date2num
from tidal_tools import extract_HC,get_tide
from res_tools import get_file_interpolator,vertical_extrapolation
from filetype import create_ncTH
import cdms2
from vcmq import fill2d,grid2xy,griddata,create_time,create_depth,create_axis,MV2,N
import scipy.io
from scipy.interpolate import griddata,interp1d




class OpenBoundaries(object):
    def __init__(self,obc,hgrid,vgrid,t0,t1,z0=0.001,logger=None):
        '''Docstring'''  

        if logger:
            self.logger = logger

        self.obc=obc
        self.hgrid=hgrid
        self.vgrid=vgrid
        self.t0=t0
        self.t1=t1

        self.bnd_nodes=[self.hgrid.mesh.boundaries[None][bnd]['indexes'] for bnd in self.hgrid.mesh.boundaries[None]]
        for x in range(len(self.bnd_nodes)):
            self.bnd_nodes[x]=[int(y)-1 for y in self.bnd_nodes[x]]

        self.tidal=False
        self.residual=False
        self.i23d=2
        self.ivs=1
        self.z0=z0
        self.lat0=np.mean(self.hgrid.latitude)

        bnd=obc.get('bnd',None)
        if bnd: 
            self.llat,self.llon,self.zz=self.get_all_open_bnd(bnd)
            self.zz=np.array(self.zz)
        else:
            self.llat,self.llon,self.zz=self.get_all_nodes()
            self.zz=np.array(self.zz)

    def get_all_nodes(self):
        Lat_ocean=[]
        Lon_ocean=[]
        Zlayer_ocean=[] 
        for node_i in range(0,len(self.hgrid.longitude)):
            Lat_ocean.append(self.hgrid.latitude[node_i])
            Lon_ocean.append(self.hgrid.longitude[node_i])
            Zlayer_ocean.append(self.vgrid.sigma_to_zlayer(node_i,self.hgrid.h[node_i]*-1,0.,0.1))

        return Lat_ocean,Lon_ocean,Zlayer_ocean 

    def get_all_open_bnd(self,bnd):
        if  type(bnd[0])==str:
            bnd=[int(x) for x in bnd[0].split()]

        Lat_ocean=[]
        Lon_ocean=[]
        Zlayer_ocean=[]

        for nn in bnd:
            ocean_boundary = self.bnd_nodes[nn-1] 
            Lat_oce=[]
            Lon_oce=[]
            Zlayer_oce=[]
            for node_i in ocean_boundary:
                Lat_ocean.append(self.hgrid.latitude[node_i])
                Lon_ocean.append(self.hgrid.longitude[node_i])
                Zlayer_ocean.append(self.vgrid.sigma_to_zlayer(node_i,self.hgrid.h[node_i]*-1,0.,0.1))

        return Lat_ocean,Lon_ocean,Zlayer_ocean        


    def add_res(self,res):

       self.res_file=cdms2.open(res['filename'])
       self.res_vars=res['vars']
       self.residual=True
       if len(self.res_vars)>1:
          self.ivs=2

    def add_tide(self,tidal):

        self.HC,self.tfreq,self.constidx=extract_HC(tidal['filename'],tidal['vars'],self.llon,self.llat, logger=self.logger)
        self.tidal=True
        if len(self.HC.keys())>1:
          self.ivs=2



    def create_Dthnc(self,fileout,TimeSeries):     
        if '2D' in fileout:
            self.i23d=2
        else:
            self.i23d=3


         
        # create file
        if self.i23d==3:
            Nlev=self.zz.shape[1]
        else:
            Nlev=1

        time_Series,nc=create_ncTH(fileout,len(self.llon),Nlev,self.ivs,np.round((TimeSeries-TimeSeries[0])*24*3600))



        for n in range(0,len(TimeSeries)):
            tin=create_time(np.ones(len(self.llon)*Nlev)*(TimeSeries[n]+1),units='days since 1-1-1')

            total=np.zeros(shape=(self.ivs,len(self.llon),Nlev))

            # get tide
            if self.tidal:
                var=self.HC.keys()

                for i,v in enumerate(sorted(var)):
                    # horizontal interpolation
                    tmp=get_tide(self.constidx,self.tfreq,self.HC[v],np.array(TimeSeries[n]),self.lat0)
        
                    if self.i23d>2:# vertical interpolation
                        tmp=vertical_extrapolation(tmp,self.zz,z0=self.z0)
                        
                    total[i,:,:]=total[i,:,:]+tmp


            if self.residual:
                var=self.res_vars

                for i,v in enumerate(sorted(var)):
                    arri=self.res_file[v][:]
                    if self.i23d >2:
                       dep=create_depth(arri.getAxis(1)[:])
                       extra = create_axis(N.arange(1), id='member')
                       arri2=np.tile(arri,[1,1,1,1,1])
                       arri3 = MV2.array(arri2,axes=[extra,arri.getAxis(0), dep, arri.getAxis(2), arri.getAxis(3)], copy=False,fill_value=1e20)
                      
                       zi=arri.getAxis(1)[:]
                       if np.mean(zi)>0:
                          zi=zi*-1
                       tb=grid2xy(arri3,xo=np.tile(self.llon,[Nlev,1]).T.flatten(),yo=np.tile(self.llat,[Nlev,1]).T.flatten(),zo=self.zz.flatten(),method='linear',to=tin,zi=zi)

                    else:
                       tb=grid2xy(arri,xo=self.llon,yo=self.llat,method='linear',to=tin)

                    if np.any(tb.mask==True):
                        bad=tb.mask==True
                        if len(bad.shape)>1:
                            bad=bad[0,:]
                        tin_bad=create_time(np.ones(len(bad))*(TimeSeries[n]+1),units='days since 1-1-1')
                        
                        if self.i23d >2:
                            llon=np.tile(self.llon,[Nlev,1]).T.flatten()
                            llat=np.tile(self.llat,[Nlev,1]).T.flatten()
                            zz=self.zz.flatten()
                            zi=arri.getAxis(1)[:]
                            if np.mean(zi)>0:
                                zi=zi*-1

                            tb[0,bad]=grid2xy(arri3,xo=llon[bad],yo=llat[bad],zo=zz[bad],method='nearest',to=tin_bad,zi=zi)

                        else:
                            tb[bad]=grid2xy(arri,xo=np.array(self.llon)[bad].tolist(),yo=np.array(self.llat)[bad].tolist(),method='nearest',to=tin_bad)


                    if np.any(tb.mask==True):
                        print('probleme')


                    total[i,:,:]=total[i,:,:]+np.reshape(tb,(len(self.llon),Nlev))



            total=np.transpose(total,(1,2,0))
            

            if np.isnan(total).any():
                import pdb;pdb.set_trace()  
    
            if n % 100 == 0:
                self.logger.info('For timestep=%.f, max=%.4f, min=%.4f , max abs diff=%.4f' % (TimeSeries[n],total.max(),total.min(),abs(np.diff(total,n=1,axis=0)).max()))
            
            time_Series[n,:,:,:]=total
  
        nc.close()

    def create_th(self,fileout,TimeSeries,options):
        dt=(TimeSeries[1]-TimeSeries[0])*24.
        DtSeries=np.arange(0,(len(TimeSeries))*dt*3600,dt*3600)



        #### check if we have a loop in the dictionnary
        Opt={}
        for nk in options.keys():
            if type(nk) is not int:
                [aa,bb]=nk.split('-')
                for nkk in range(int(aa),int(bb)+1):
                    Opt[int(nkk)]=options[nk]                   
            else:
                Opt[int(nk)]=options[nk]
        options=Opt
        Y=np.ones(shape=(len(TimeSeries),len(options)+1))
        Y[:,0]=DtSeries
        fmt='%.f    '
        
        


        for n in range(0,len(options)):
            if 'filename' in options[n+1]:
                mat = scipy.io.loadmat(options[n+1]['filename'])
                x=mat[options[n+1]['X']]-366
                y=mat[options[n+1]['Y']].flatten(1)
                
                if np.mean(y)>=0 and fileout[-5:]=='ux.th':
                    y=y*-1

                if 'fac' in options[n+1]:
                    fac=options[n+1]['fac']
                    y=y+(y*(fac/100.0))
            
                y2 = interp1d(x.flatten(1), y, kind='linear')
                Y[:,n+1]=y2(TimeSeries)
                fmt=fmt+'%.4f   '
            elif type(options[n+1]['Y'])==type(list()):
                monthly_time = daterange(num2date(TimeSeries[0])-relativedelta(months=1), num2date(TimeSeries[-1])+relativedelta(months=1),delta=1,typ='months')
                data=options[n+1]['Y']
                starting_month=monthly_time[0].month-1
                ending_month=monthly_time[-1].month
                monthly_data=data[starting_month:ending_month]
                # for y in range(1,len(np.unique([x.year for x in monthly_time]))-1):
                #   monthly_data+=data
                
                # monthly_data+=data[:ending_month]

                y2 = interp1d(date2num(monthly_time).flatten(1), monthly_data, kind='cubic')
                Y[:,n+1]=y2(TimeSeries)
                fmt=fmt+'%.4f   '

            else:
                Y[:,n+1]=Y[:,n+1]*options[n+1]['Y']
                fmt=fmt+'%.4f   '

        fmt=fmt[:-1]
        np.savetxt(fileout, Y,fmt=fmt)

    def make_boundary(self,filename,dt=3600):
        if self.logger:
            self.logger.info("  Writing %s" %filename)
        TimeSeries=np.arange(date2num(self.t0),date2num(self.t1)+1,dt/(24.*3600.))
       
        if filename.endswith('.th.nc') or filename.endswith('_nu.nc'):
            self.create_Dthnc(filename,TimeSeries)
        elif filename.endswith('.th'):
            self.obc.pop('dt')
            self.create_th(filename,TimeSeries,self.obc)







