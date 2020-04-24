import os
import numpy as np
#import regrid
from netCDF4 import Dataset
import cdms2
#from scipy.interpolate import griddata,interp1d
from vcmq import fill2d,grid2xy,griddata,create_time
from vacumm.misc.grid.regridding import GridData
from vacumm.misc.grid.kriging import krig

def vertical_extrapolation(u,dep,lev=[],z0=[]):
	zu=np.zeros(shape=(dep.shape))

	if u.shape[1]<2: ### 2D dataset

		if dep.shape[1]>2:
			for nlev in range(0,dep.shape[1]):
				if z0!=-99:
					### 2D dataset > log progile for 2D data
					dz=abs(dep[:,0])-abs(dep[:,nlev])
					good=np.nonzero(dz!=0)[0]
					
					if len(good)>1:
						
						zu[good,nlev]=u[good].T*(np.log(dz[good]/z0))/(np.log(abs(dep[good,0])/z0)-1)
				else:
					zu[:,nlev]=u
				

		else:
			zu[:,1]=u
			zu[:,0]=u # not too sure for this one maybe should be 0
	else: # 3D datasets
		
		for n in range(0,dep.shape[0]):
			not_nan = ~np.isnan(u[:,n])
			
			iz=zinterp1d(-lev[not_nan], dep[n,:], u[not_nan,n] ).squeeze()
		
			if any(iz>1e10):
				mask=np.zeros(shape=u[:,n].shape)
				mask[u[:,n]>1e10]=1
				var=np.ma.array(u[:,n],mask=mask)
				varout=interpmiss(var.data,var.mask,4)
				iz=zinterp1d(-lev[:], dep[n,:], varout).squeeze()
				
			zu[n,:] = iz
	bad=dep==-999
	zu[bad]=-999	
	return zu
	
def missing_interp(lon,lat,arr):
	lon,lat=np.meshgrid(lon,lat)
	tempf = fill2d(arr, method='carg')
	lon=np.ma.masked_where(tempf[0,:,:]==tempf.fill_value, lon)
	lat=np.ma.masked_where(tempf[0,:,:]==tempf.fill_value, lat)
	

	return lon,lat,tempf

def get_file_interpolator(sourcefile,var_res,lonr,latr):
	sourcedata = Dataset(sourcefile,'r')
	sourcedata_dms2 = cdms2.open(sourcefile)
	f_out={}
	LONR={}
	LATR={}
	time0=[t.torelative('days since 1-1-1').value for t in sourcedata_dms2['time'].asRelativeTime()]
	tres=np.array(time0)



	for varname in var_res:
		#	arri = griddata(lat.ravel(),lon.ravel(), tempf[0,:,:].compressed(), tempf.getGrid(), method='nearest')
		lon,lat,arri=missing_interp(sourcedata_dms2[varname].getAxis(2)[:],sourcedata_dms2[varname].getAxis(1)[:],sourcedata_dms2[varname][:])
		import pdb;pdb.set_trace()
		tt=create_time(time0,units='days since 1-1-1')
		tb=grid2xy(arri,xo=lonr,yo=latr,method='nearest',to=tt)

		Res_val = np.ma.masked_where(tempf>=1e36,tempf)
		if len(Res_val.shape)>3:
			LONR[varname]=[]
			LATR[varname]=[]
			f_out[varname]=[]

			for all_l in range(0,Res_val.shape[1]):


				test=Res_val[:,all_l,:,:].reshape(Res_val.shape[0],Res_val.shape[2]*Res_val.shape[3])
				gd_node=test[0,:].nonzero()[0]
				test=np.array(test)
				test = np.ma.masked_where(test==Res_val.fill_value,test)
				test=test[:,gd_node]
				if test.shape[1]>0:
					gd_ts=test[:,int(test.shape[1]/2)].nonzero()[0]
					f_out[varname].append( interp1d(tres[gd_ts], test[gd_ts], axis=0,fill_value=np.nan))

					LONR[varname].append(lonr[gd_node])
					LATR[varname].append(latr[gd_node])
		else:
			test=Res_val[:,:,:].reshape(Res_val.shape[0],Res_val.shape[1]*Res_val.shape[2])
			gd_node=test[0,:].nonzero()[0]
			test=np.array(test)
			test = np.ma.masked_where(test>=1e20,test)
			test=test[:,gd_node]
			gd_ts=test[:,int(test.shape[1]/2)].nonzero()[0]
			LONR[varname]=lonr[gd_node]
			LATR[varname]=latr[gd_node]
			f_out[varname] = interp1d(tres[gd_ts], test[gd_ts], axis=0,fill_value=np.nan)#'extrapolate')



	return f_out,LONR,LATR
