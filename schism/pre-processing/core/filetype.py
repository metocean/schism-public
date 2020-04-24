#!/usr/bin/env python2.7
from netCDF4 import Dataset

def create_ncTH(filename,Nnode,Nlev,Nvar,T):
	if '_nu.nc' not in filename:

		nc = Dataset(filename, 'w', format='NETCDF4',clobber=True)
		# dimensions
		nc.createDimension('nOpenBndNodes', Nnode )
		nc.createDimension('nLevels', Nlev )
		nc.createDimension('nComponents', Nvar)
		nc.createDimension('one', 1)
		nc.createDimension('time', len(T))


		time_step = nc.createVariable('time_step', 'f4', ('one',))
		time_step.long_name="time step in seconds"
		time = nc.createVariable('time', 'f8', ('time',))
		time.long_name="simulation time in seconds"
		time_series = nc.createVariable('time_series', 'f4', ('time','nOpenBndNodes','nLevels','nComponents'))
		time_step[0]=T[1]-T[0]
		time[:]=T
	else:
		nc = Dataset(filename, 'w', format='NETCDF4',clobber=True)
		# dimensions
		nc.createDimension('node', Nnode )
		nc.createDimension('nVert', Nlev )
		nc.createDimension('ntracers', Nvar)
		nc.createDimension('time', len(T))
		time = nc.createVariable('time', 'f8', ('time',))
		time.long_name="simulation time in days"
		time_series = nc.createVariable('tracer_concentration', 'f4', ('time','node','nVert','ntracers'))
		time[:]=T/86400.0

	return time_series,nc

def create_hotstartNC(hot,filename):
	nc = Dataset(filename, 'w', format='NETCDF4',clobber=True)
	# dimensions
	ntracer=hot['tr_nd'].shape[0]
	nc.createDimension('node', hot['np'] )
	nc.createDimension('elem', hot['ne'])
	nc.createDimension('side', hot['ns'])
	nc.createDimension('nVert', hot['nvrt'])
	nc.createDimension('ntracers',ntracer)
	nc.createDimension('one', 1)
	nc.createDimension('three', 3)

	if 'SED_Nbed' in hot:
		nc.createDimension('SED_Nbed', hot['SED_Nbed'])
		nc.createDimension('SED_MBEDP', hot['SED_MBEDP'])
		nc.createDimension('SED_ntr', hot['SED_ntr'])

	var={}

	var['time'] = nc.createVariable('time', 'f8', ('one',))
	var['iths'] = nc.createVariable('iths', 'i4', ('one',))
	var['ifile'] = nc.createVariable('ifile', 'i4', ('one',))
	var['idry_e'] = nc.createVariable('idry_e', 'i4', ('elem',))
	var['idry_s'] = nc.createVariable('idry_s', 'i4', ('side',))
	var['idry'] = nc.createVariable('idry', 'i4', ('node',))
	var['eta2'] = nc.createVariable('eta2', 'f8', ('node',))
	var['we'] = nc.createVariable('we', 'f8', ('elem','nVert'))
	var['tr_el'] = nc.createVariable('tr_el', 'f8', ('elem','nVert','ntracers'))
	var['su2'] = nc.createVariable('su2', 'f8', ('side','nVert'))
	var['sv2'] = nc.createVariable('sv2', 'f8', ('side','nVert'))
	var['tr_nd'] = nc.createVariable('tr_nd', 'f8', ('node','nVert','ntracers'))
	var['tr_nd0'] = nc.createVariable('tr_nd0', 'f8', ('node','nVert','ntracers'))
	var['q2'] = nc.createVariable('q2', 'f8', ('node','nVert')) #diffusivity shit
	var['xl'] = nc.createVariable('xl', 'f8', ('node','nVert')) #diffusivity shit
	var['dfv'] = nc.createVariable('dfv', 'f8', ('node','nVert')) #Compute turbulence diffusivities dfv, dfh
	var['dfh'] = nc.createVariable('dfh', 'f8', ('node','nVert')) #Compute turbulence diffusivities dfv, dfh
	var['dfq1'] = nc.createVariable('dfq1', 'f8', ('node','nVert')) # diffmin
	var['dfq2'] = nc.createVariable('dfq2', 'f8', ('node','nVert')) #diffmax

	if 'SED_Nbed' in hot:
		var['SED3D_dp'] = nc.createVariable('SED3D_dp', 'f8', ('node')) #dp
		var['SED3D_rough'] = nc.createVariable('SED3D_rough', 'f8', ('node')) #rough
		var['SED3D_bed'] = nc.createVariable('SED3D_bed', 'f8', ('elem','SED_Nbed','SED_MBEDP')) #property
		var['SED3D_bedfrac'] = nc.createVariable('SED3D_bedfrac', 'f8', ('elem','SED_Nbed','SED_ntr')) #fraction


	for key in var.keys():
		Z=hot.get(key,'0')
		if type(Z)!=type(int()):
			try:
				Z=Z.T
			except:
				import pdb;pdb.set_trace()
		var[key][:]=Z

	nc.close()
