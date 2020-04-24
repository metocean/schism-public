#!/usr/bin/env python2.7
import sys,os

import netCDF4
import numpy as np

from ttide.t_tide import t_tide
from ttide.t_getconsts import t_getconsts
from ttide.t_vuf import t_vuf
from ttide import t_utils as tu

from matplotlib.dates import datestr2num,date2num,num2date
from datetime import timedelta
from pyproj import Proj, transform

import interpz

outputEPSG=4326

PREFIX='schout_'
SUFIX=''



def Calculate(mean_lat,Tstart,Cons):
    const=t_getconsts(np.array([]))
    itemindex=[]
    for k in range(0,len(Cons)):
        itemindex .append(np.where(const[0]['name']==Cons[k].ljust(4))[0][0])

    # Earth equilibrium at the beginning of the run
    [v,u,f]=t_vuf('nodal',Tstart,itemindex, mean_lat);
    tear=360.*(v+u);

    # Nodel factor at the middl of the run
    [v,u,f]=t_vuf('nodal',Tstart,itemindex, mean_lat);
    tnf=f;

    tfreq=2.*np.pi*const[0]['freq'][itemindex]/3600.;
    talpha=const[0]['name'][itemindex];
    return talpha,tfreq,tear,tnf


def create_filetide_var(filetide,t0,vari,cons,X,Y,depth,nv,do_variance,lev):
    FillValue=-99999
    ncons=len(cons);
    nconstr=4;

    ncout = netCDF4.Dataset(filetide,  'w')

    # dimensions
    dcons = ncout.createDimension('cons', ncons)
    dconstr  = ncout.createDimension('constr', nconstr)
    dnode   = ncout.createDimension('nSCHISM_hgrid_node', len(depth))
    dnface   = ncout.createDimension('nMaxSCHISM_hgrid_face_nodes', 4)
    dface   = ncout.createDimension('nSCHISM_hgrid_face', nv.shape[0])

    if lev is not None:
        dlev   = ncout.createDimension('lev', len(lev))


    # write variable

    vcons = ncout.createVariable('cons','S1',dimensions=('cons','constr'))
    vcons.long_name='Tide constituents'
    vcons.standard_name='tide constituents'

    vX=ncout.createVariable('Easting',np.float64,('nSCHISM_hgrid_node',))
    vY=ncout.createVariable('Northing',np.float64,('nSCHISM_hgrid_node',))
    vD=ncout.createVariable('Depth',np.float64,('nSCHISM_hgrid_node',))

    vN=ncout.createVariable('ele',np.float64,('nSCHISM_hgrid_face','nMaxSCHISM_hgrid_face_nodes',))

    vfreq=ncout.createVariable('freq',np.float64,('cons',))
    vfreq.long_name='freq'

    if lev is not None:
        vlev=ncout.createVariable('lev',np.float64,('lev',))


    all_mean={}
    all_exp={}
    all_amp={}
    all_pha={}
    for nvar in vari:
        if nvar is 'u' or nvar is 'v':
            all_amp[nvar]=ncout.createVariable(nvar+'_amp',np.float64,('cons','lev','nSCHISM_hgrid_node'))
            all_amp[nvar].fill_value=FillValue

            all_pha[nvar]=ncout.createVariable(nvar+'_pha',np.float64,('cons','lev','nSCHISM_hgrid_node',))
            all_pha[nvar].fill_value=FillValue
            all_pha[nvar].units='radians'

            if do_variance is True:
                all_mean[nvar]= ncout.createVariable(nvar+'_mean',np.float64,('lev','nSCHISM_hgrid_node',))
                all_mean[nvar].fill_value=FillValue

                all_exp[nvar]= ncout.createVariable(nvar+'_var',np.float64,('lev','nSCHISM_hgrid_node',))
                all_exp[nvar].fill_value=FillValue

        else:
            all_amp[nvar]=ncout.createVariable(nvar+'_amp',np.float64,('cons','nSCHISM_hgrid_node',))
            all_amp[nvar].fill_value=FillValue

            all_pha[nvar]=ncout.createVariable(nvar+'_pha',np.float64,('cons','nSCHISM_hgrid_node',))
            all_pha[nvar].fill_value=FillValue
            all_pha[nvar].units='radians'

            if do_variance is True:
                all_mean[nvar]= ncout.createVariable(nvar+'_mean',np.float64,('nSCHISM_hgrid_node',))
                all_mean[nvar].fill_value=FillValue

                all_exp[nvar]= ncout.createVariable(nvar+'_var',np.float64,('nSCHISM_hgrid_node',))
                all_exp[nvar].fill_value=FillValue




    for j in range(0,ncons):
        for k in range(0,len(cons[j])):
            vcons[j,k]=cons[j][k]


        
  	# write V0, freq 

    const, sat, cshallow = t_getconsts(np.array([t0]))
    freq=np.zeros(shape=(len(cons),1))
    i=0
    for ic in const['name']:
        if ic in cons:
            IDX=np.nonzero(ic==cons)[0][0]
            freq[IDX]=2*np.pi*24*const['freq'][i]       	
        i=i+1
	
    vfreq[:]=freq 
    # for nvar in vari:
    #     all_amp[nvar][:]=np.zeros([len(depth),ncons])
    #     all_pha[nvar][:]=np.zeros([len(depth),ncons])

    vX[:]=X
    vY[:]=Y
    vD[:]=depth
    vN[:]=nv


    if lev is not None:
        vlev[:]=np.abs(lev)

    return ncout,all_amp,all_pha,all_mean,all_exp

def initialize(hafreq,nodes,vari,do_variance,levs):
    mnharf=len(hafreq) 
    if hafreq[0]==0:
         nz=0
         nf=1
    else:
         nz=1
         nf=0

    nfreq=len(hafreq)-nf
    mm=2*nfreq+nf
    ha=np.zeros((mm,mm))
    icall=0
    globvar={}
    globAV={}
    globVA={}
    for nv in vari:
        if nv is 'u' or nv is 'v':
            globvar[nv]=np.zeros((mm,nodes,len(levs)))
            if do_variance is True:
                globAV[nv]=np.zeros((nodes,len(levs)))
                globVA[nv]=np.zeros((nodes,len(levs)))           
        else:
            globvar[nv]=np.zeros((mm,nodes))
            if do_variance is True:
                globAV[nv]=np.zeros((nodes))
                globVA[nv]=np.zeros((nodes))

    return icall,ha,nfreq,nf,nz,globvar,mm,mnharf,globAV,globVA

def LSQUPDEG(nf,nz,gloelv,elv,nfreq,hafreq,timeud):

    if nz==0:
        gloelv[0,:]=gloelv[0,:]+elv


    for i in range(1,nfreq+1):
        i1=(2*i-nz)-1
        i2=(i1+1)
        tf1=hafreq[i+nf-1]*timeud
        ctf1 = np.cos(tf1)
        stf1 = np.sin(tf1)
        gloelv[i1,:]=gloelv[i1,:]+elv*ctf1
        gloelv[i2,:]=gloelv[i2,:]+elv*stf1

    return gloelv



def LSQUPDLHS(nf,nfreq,icall,ha,hafreq,time):

#     Take care of the steady constituent if included in the analysis
    icall+=1
    if nf==1:
        ha[0,0]=icall
        for j in range(1,nfreq+1):
            tf1=hafreq[j+nf-1]*time
            ha[0,(2*j)-1]   = ha[0,(2*j)-1] + np.cos(tf1)
            ha[0,(2*j+1)-1] = ha[0,(2*j+1)-1] + np.sin(tf1)


#   Take care of the other constituents

    for i in range(1,nfreq+1):
        for j in range(1,nfreq+1):
            i1=2*i-(1-nf)
            i2=i1+1
            j1=2*j-(1-nf)
            j2=j1+1
            tf1=hafreq[i+nf-1]*time
            tf2=hafreq[j+nf-1]*time
            ha[i1-1,j1-1] = ha[i1-1,j1-1] + np.cos(tf1)*np.cos(tf2)
            ha[i1-1,j2-1] = ha[i1-1,j2-1] + np.cos(tf1)*np.sin(tf2)
            ha[i2-1,j2-1] = ha[i2-1,j2-1] + np.sin(tf1)*np.sin(tf2)
            if i2-1<=j1-1:
                ha[i2-1,j1-1] = ha[i2-1,j1-1] + np.sin(tf1)*np.cos(tf2)


     
    return icall,ha

def fulsol(idecom,mnharf,mm,ha,hap):
        

    if idecom==0:

        for j in range(0,mm):
            for i in range(j,mm):
                ha[i,j]=ha[j,i]

         
#     Decomposition of matrix a
		
        for ir in range(0,mm):
            ire=ir+1
            for j in range(ire,mm):
                ha[ir,j]=ha[ir,j]/ha[ir,ir]
            # if (ire>mm):
            #     return ha
            for j in range(ire,mm):
                for k in range(ire,mm):
                    ha[k,j]=ha[k,j]-ha[k,ir]*ha[ir,j]
            for j in range(ire,mm):
                ha[j,ir]=0.0

        return ha
    
      

#...  solve for y by forward substitution for l*y=p
    c=np.zeros((2*mnharf,1))
    y=np.zeros((2*mnharf,1))
    hax=np.zeros((2*mnharf,1))

    n  = ha.shape[0]
    c  = np.zeros((n),float)
    y  = np.zeros((n),float)
    for i in range(0,n):           # forward substitution -
         s = hap[i]
         for j in range(0,i):
              s = s - ha[i,j]*y[j]
         y[i] = s

    for i in range(n-1,-1,-1):     # backward substitution  - 
         s = y[i]
         for j in range(i+1,n):
              s = s - ha[i,j]*hax[j]
         hax[i] = s/ha[i,i]




    return hax

def LSQSOLEG(mnharf,nf,mm,ha,gloelv,haff,haface):

    convrd=180./np.pi

    phasee=np.zeros((mnharf,gloelv.shape[1]))
    emag=np.zeros((mnharf,gloelv.shape[1]))

#***** AT each node transfer each load vector to p, solve and write output
    nnodes=gloelv.shape[1]
    for n in range(0,nnodes):
        hax=fulsol(1,mnharf,mm,ha,gloelv[:,n])
        if n % 10000 == 0:
                print '      nodes %i of %i' % (n,nnodes)
#        Compute amplitude and phase for each frequency making sure that the
#!        phase is between 0 and 360 deg.  Then write output.

       
        if nf==1:
            emag[0,n]=hax[0]/haff[0]
            phasee[0,n]=0

        for i in range(2,hax.shape[0]/2+1):
                i1=(2*i-1-nf)-1
                i2=(i1+1)
                emag[i-1,n]=np.sqrt(hax[i1]*hax[i1]+hax[i2]*hax[i2])/haff[i-1]

                if hax[i1]==0 and hax[i2]==0:
                	phasee[i-1,n]=0
                else:
                	phasee[i-1,n] = np.arctan2(hax[i2],hax[i1])
                	phasee[i-1,n] =convrd*phasee[i-1,n] +haface[i-1]
                if phasee[i-1,n] <0:
                	phasee[i-1,n]=phasee[i-1,n]+360.

                if phasee[i-1,n]>=360:
                	phasee[i-1,n]=phasee[i-1,n]-360.
      

    return emag,phasee


def get_variance(elav,elva,hafreq,emag,phasee,ntsteps,dt):

    eav = 0.
    esq = 0.
    TIMEBEG=0.

    for it in range(0,ntsteps):
        if it % 1000 == 0:
            print '      timesteps %i of %i' % (it,ntsteps)
        time=TIMEBEG+dt*it
        rse=np.zeros((emag.shape[1]))
        for nf in range(0,hafreq.shape[0]):
            ftime=hafreq[nf]*time
            rse=rse+emag[nf,:]*np.cos(ftime-phasee[nf])
        eav=eav+rse
        esq=esq+rse**2

         
    eav=eav/ntsteps
    esq=esq/ntsteps-eav*eav

    eavdif=np.zeros((emag.shape[1]))
    gd=elav!=0
    eavdif[~gd]=1.0
    eavdif[gd]=eav[gd]/elav[gd]

    evadif=np.zeros((emag.shape[1]))
    gd=elva!=0
    evadif[~gd]=1.0
    evadif[gd]=esq[gd]/elva[gd]

    return eavdif,evadif

def read_initial_netcdf_file(file0,EPSG,is_3d):
    print 'READING INITIAL FILE: %s' % file0

    nc=netCDF4.Dataset(file0)

    unit=nc.variables['time'].units
    depth=nc.variables['depth'][:]
    X=nc.variables['SCHISM_hgrid_node_x'][:]
    Y=nc.variables['SCHISM_hgrid_node_y'][:]
    nv=nc.variables['SCHISM_hgrid_face_nodes'][:]

    t0 = netCDF4.num2date(nc.variables['time'][0],unit)
    tend = netCDF4.num2date(nc.variables['time'][-1],unit)

    inProj = Proj(init='epsg:'+str(EPSG))
    outProj = Proj(init='epsg:'+str(outputEPSG))

    mean_lon,mean_lat=transform(inProj,outProj,float(np.mean(X[:])),float(np.mean(Y[:]))) 
    dt=nc.variables['time'][2]-nc.variables['time'][1]

    Tstart=date2num(t0-timedelta(seconds=dt))

    nc.close()

    return mean_lat,Tstart,X,Y,depth,nv,t0,dt
def vertical_interpolation(zcor,e,lev):
    E=np.zeros((e.shape[0],len(lev)))

    if lev[0]==0.:
        E[:,0]=e[:,-1]
        lev2=lev[1:]
    else:
        lev2=lev
    
    e2=e.data
    e2[e.data==e.fill_value]=0
    z2=zcor.data;
    z2[z2>1e36]=-9999
    z2[np.isnan(z2)]=-9999
    surf=np.zeros((zcor.shape[0],1))
    surf[:,0]=z2[:,-1]
    z2=z2-surf
      
    EE=interpz.interpz1d(e2,z2,lev2,np=e2.shape[0],nzin=e2.shape[1],nzout=len(lev2),kz=1, null_value=-9.9e15)
    
    for n in range(0,len(lev2)):
        E[:,n+1]=EE[:,n]



    return E
def check_rayleigh(ndays,ray):
    tmptuple = tu.constituents(ray / (ndays*24.),\
                               np.array([]), np.array([]),\
                               np.array([]), np.array([]),\
                               np.array([]))
    nameu, fu, ju, namei, fi, jinf, jref = tmptuple

    return [x.replace(' ','') for x in nameu]

def process(fileout,dirout,INDstart,INDend,vari,Cons,EPSG,do_variance,lev,ray,ndays):


    gd_cons=check_rayleigh(ndays,ray)


    for con in Cons[1:]:
        if con not in gd_cons:
            raise argparse.ArgumentTypeError('%s violate raylight criterion:\n\
             Change criterion value or length of timesries or specify another constituent' % con)
        else:
            'Using %s'% con



    if 'u' in vari:
        is_3d=True
    else:
        is_3d=False



    mean_lat,Tstart,X,Y,depth,nv,t0,dt=read_initial_netcdf_file(os.path.join(dirout,PREFIX+str(INDstart)+SUFIX+'.nc'),EPSG,is_3d)

    [talpha,tfreq,tear,tnf]=Calculate(mean_lat,np.array([Tstart]),Cons)


    icall,ha,nfreq,nf,nz,globvar,mm,mnharf,globAV,globVA=initialize(tfreq,len(X),vari,do_variance,lev)
  
    
    ncout,all_amp,all_pha,all_mean,all_exp = create_filetide_var(fileout,t0.toordinal(), vari,talpha,X,Y,depth,nv,do_variance,lev)

    print 'LOADING DATA'

    time_in_sec=0

    for ifile in range(INDstart,INDend+1):
        file1=os.path.join(dirout,PREFIX+str(ifile)+SUFIX+'.nc')
        nc=netCDF4.Dataset(file1)
        nt=len(nc.variables['time'])
        print '   READING FILE: %s' % file1

        for n in range(0,nt):

            time_in_sec=time_in_sec+dt # Hopefully the timnestep never change!!!
            
            #...  ..UPDATE THE LHS MATRIX
            icall,ha=LSQUPDLHS(nf,nfreq,icall,ha,tfreq,time_in_sec)

            if is_3d:
                zcor=nc.variables['zcor'][n,:,:]

            for nv in vari:
                if nv == 'u' or nv == 'v':
                    if nv=='u':
                        e=nc.variables['hvel'][n,:,:,0]
                    elif nv=='v':
                        e=nc.variables['hvel'][n,:,:,1]

                    E=vertical_interpolation(zcor,e,lev)
                    for ilev in range(0,len(lev)):
                        #....IF DESIRED UPDATE GLOBAL ELEVATION LOAD VECTOR
                        globvar[nv][:,:,ilev]=LSQUPDEG(nf,nz,globvar[nv][:,:,ilev],E[:,ilev],nfreq,tfreq,time_in_sec)
                        if do_variance is True:
                            globAV[nv][:,ilev]=globAV[nv][:,ilev]+E[:,ilev]
                            globVA[nv][:,ilev]=globVA[nv][:,ilev]+E[:,ilev]**2


                else:                   
                    if nv=='e':
                        E=nc.variables['elev'][n,:]
                    elif nv=='um':
                        E=nc.variables['dahv'][n,:,0]
                    elif nv=='vm':
                        E=nc.variables['dahv'][n,:,1]

                    #....IF DESIRED UPDATE GLOBAL ELEVATION LOAD VECTOR
                    globvar[nv]=LSQUPDEG(nf,nz,globvar[nv],E,nfreq,tfreq,time_in_sec)

                    if do_variance is True:
                        globAV[nv]=globAV[nv]+E
                        globVA[nv]=globVA[nv]+E**2

        nc.close()



    #......Fill out and decompose the LHS harmonic analaysis matrix

    hax=fulsol(0,mnharf,mm,ha,[])
    if do_variance is True:
        for nv in vari:
            globAV[nv]=globAV[nv]/icall
            globVA[nv]=globVA[nv]/icall-globAV[nv]**2


    #......Solve the harmonic analysis problem and write the output
    #
    for nv in vari:
        print 'SOLVING %s' % nv
        if nv is 'u' or nv is 'v':
            for ilev in range(0,len(lev)):
                emag,phasee=LSQSOLEG(mnharf,nf,mm,hax,globvar[nv][:,:,ilev],tnf,tear)
                all_amp[nv][:,ilev,:]=emag
                all_pha[nv][:,ilev,:]=phasee           
                if do_variance is True:
                    print 'resynthesized time series for %s' % nv
                    e_mean,e_exp=get_variance(globAV[nv][:,ilev],globVA[nv][:,ilev],tfreq,emag,phasee,icall,dt)
                    all_exp[nv][ilev,:]=e_exp
                    all_mean[nv][ilev,:]=e_mean

        else:          
            emag,phasee=LSQSOLEG(mnharf,nf,mm,hax,globvar[nv],tnf,tear)
            all_amp[nv][:]=emag
            all_pha[nv][:]=phasee
            
            if do_variance is True:
                print 'resynthesized time series for %s' % nv
                e_mean,e_exp=get_variance(globAV[nv],globVA[nv],tfreq,emag,phasee,icall,dt)
                all_exp[nv][:]=e_exp
                all_mean[nv][:]=e_mean

    ncout.close()
	

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(prog='tidal_analysis.py', usage='%(prog)s fileout dirout INDstart INDend params [options]')
    ## main arguments
    parser.add_argument('fileout', type=str,help='name of the output tidal file')
    parser.add_argument('dirout', type=str,help='Folder with all the SCHSIM files')
    parser.add_argument('INDstart', type=int,help='First file to take')
    parser.add_argument('INDend', type=int,help='Last file to take')
    parser.add_argument('params', type=str,nargs='+',help='Parameters (i.e e um vm')

    ## options  
    parser.add_argument('--EPSG', type=int,help='EPSG number (2193 for NZTM)')
    parser.add_argument('--VAR', type=str2bool, help='calculate variance')
    parser.add_argument('--CONS', type=str,nargs='+',help='Constituents (i.e e um vm u v')
    parser.add_argument('--LEV', type=float,nargs='+',help='outputs level 0 -1 -2 ... ')
    parser.add_argument('--RAY', default=1, type=float,help='Rayleight criterion default=1 ')
    parser.add_argument('--NDAYS',default=220, type=float,help='Length of the timeseries in days')

    args = parser.parse_args()
    
    if args.EPSG is None:
        args.EPSG=2193
    if args.VAR is None:
        args.VAR=True
    if args.CONS is None:
        args.CONS=['Z0','M2','S2','N2','SK3','K1','M3','MSF','MS4','M4','L2','O1']


    if args.CONS[0]!='Z0':
        raise argparse.ArgumentTypeError('The first constituent must be Z0')

    if args.LEV is None and 'u' in args.params:
        raise argparse.ArgumentTypeError('Needs output levels')


    process(args.fileout,\
        args.dirout,\
        args.INDstart,args.INDend,\
        args.params,\
        args.CONS,\
        args.EPSG,args.VAR,\
        args.LEV,\
        args.RAY,args.NDAYS)

