	#!/usr/bin/env python2.7

import os,sys
from glob import glob
import netCDF4
import numpy as np
import copy
	

def process_method1(ncout,all_files,params,stats,nnodes,part):
	perc=[]
	for S in stats:
		perc.append(float(S.replace('p','')))

	for param in params:

		for nn in range(0,nnodes,part):

			for i,file in enumerate(all_files):
				if nn % 5000==0:print '%s => read %s, nodes: %.f over %.f' % (param,file,nn,nnodes)
				ncin = netCDF4.Dataset(file,  'r')
				D=ncin.variables[param][:,nn:nn+part]
				ncin.close()

				
				if i==0:
					matrix=copy.deepcopy(D)
				else:
					matrix=np.vstack((matrix,D))


			matrix[matrix==D.fill_value]=np.nan
			#import pdb;pdb.set_trace()
			matrix=np.nanpercentile(matrix,perc,axis=0,interpolation='nearest')

			for istat,stat in enumerate(stats):	
				ncout.variables[param+'_'+stat][0,nn:nn+part]=np.squeeze(matrix[istat,:])

	return ncout

					
						
def process_method2(ncout,all_files,params,stats,nnodes):
	

	for param in params:	
		for i,file in enumerate(all_files):
			print '%s => %s' % (param,file)
			ncin = netCDF4.Dataset(file,  'r')	


			D=ncin.variables[param][:]
			ncin.close()
			little_divier=np.sum(D.mask==False,0)
			if type(little_divier)==type(np.int64()):little_divier=D.shape[0]
			if i==0:
				divider=np.ones(shape=D.shape[1:])*0
				matrix=np.ones(shape=(len(stats),)+D.shape[1:])*0
			divider=divider+little_divier
					
			for istat,stat in enumerate(stats):
				if i==0 and stat=='max':
					matrix[istat,:]=matrix[istat,:]-np.inf
				elif i==0 and stat=='min':
					matrix[istat,:]=matrix[istat,:]+np.inf
				elif i==0 and stat=='mean':
					matrix[istat,:]=matrix[istat,:]*0
					stat='sum'
				elif stat=='mean':
					stat='sum'

				str_stat='np.nan'+stat+'(np.stack((matrix[istat,:],np.nan'+stat+'(D,0))),0)'

				matrix[istat,:]=eval(str_stat)


		for istat,stat in enumerate(stats):	
			if stat=='mean':
				divider[divider==0]=np.nan
			else:
				divider=1

		
			
			ncout.variables[param+'_'+stat][:]=np.squeeze(matrix[istat,:]/divider)

	
	return ncout

def create_output_file(filout,first_file,stats,params):

	ncin = netCDF4.Dataset(first_file,  'r')
 	y = ncin.variables['SCHISM_hgrid_node_y'][:]
	x= ncin.variables['SCHISM_hgrid_node_x'][:]
	nnodes=len(x)
	ele=ncin.variables['SCHISM_hgrid_face_nodes'][:]
	sigma=ncin.dimensions['nSCHISM_vgrid_layers'].size

	ncout = netCDF4.Dataset(filout,  'w')
        
	dnode = ncout.createDimension('nSCHISM_hgrid_node', len(x))
	dele = ncout.createDimension('nSCHISM_hgrid_face', len(ele))
	dface = ncout.createDimension('nMaxSCHISM_hgrid_face_nodes', 4)
	dtime =ncout.createDimension('time', 1)
	dsigma =ncout.createDimension('nSCHISM_vgrid_layers', sigma)	
	done =ncout.createDimension('one', 1)
	dtwo =ncout.createDimension('two', 2)

	vtime = ncout.createVariable('time',np.float64,dimensions=('time'))
	vele = ncout.createVariable('SCHISM_hgrid_face_nodes',np.int32,dimensions=('nSCHISM_hgrid_face','nMaxSCHISM_hgrid_face_nodes'))
	vx = ncout.createVariable('SCHISM_hgrid_node_x',np.float64,dimensions=('nSCHISM_hgrid_node'))
	vy = ncout.createVariable('SCHISM_hgrid_node_y',np.float64,dimensions=('nSCHISM_hgrid_node'))
	vdepth = ncout.createVariable('depth',np.float64,dimensions=('nSCHISM_hgrid_node'))
	vsigma = ncout.createVariable('zcor',np.float64,dimensions=('time','nSCHISM_hgrid_node','nSCHISM_vgrid_layers'))



	vtime[:]=1
	vele[:]=ncin.variables['SCHISM_hgrid_face_nodes'][:]
	vx[:]=ncin.variables['SCHISM_hgrid_node_x'][:]
	vy[:]=ncin.variables['SCHISM_hgrid_node_y'][:]
	vdepth[:]=ncin.variables['depth'][:]

	for j,param in enumerate(params):
		i23d=ncin.variables[param].i23d
		ivs=ncin.variables[param].ivs
		for S in stats:
			print param+'_'+S
			if ivs==1 and i23d==1: # elev
				ncout.createVariable(param+'_'+S,np.float64,dimensions=('time','nSCHISM_hgrid_node'))
			if ivs==2 and i23d==1: # dahv
				ncout.createVariable(param+'_'+S,np.float64,dimensions=('time','nSCHISM_hgrid_node','two'))
			if ivs==1 and i23d==2: # zcor
				ncout.createVariable(param+'_'+S,np.float64,dimensions=('time','nSCHISM_hgrid_node','nSCHISM_vgrid_layers'))
			if ivs==2 and i23d==2: # hvel
				ncout.createVariable(param+'_'+S,np.float64,dimensions=('time','nSCHISM_hgrid_node','nSCHISM_vgrid_layers','two'))					

	
	ncin.close()
	return ncout,nnodes
	
	

def process(fileout,dirin,prefix,Istart,Iend,params,stats,part):


	all_files_tmp = [y for x in os.walk(dirin) for y in glob(os.path.join(x[0], prefix+'_*.nc'))]
	all_files=[]
	for file in all_files_tmp:
		[tmp,filenum]=file.replace('.nc','').split('_')
		if int(filenum)>=Istart and int(filenum)<=Iend:
			all_files.append(file)


	print "%.f files found" % len(all_files)
	ncin,nnodes=create_output_file(fileout,all_files[0],stats,params)

	stat_method1=[]
	stat_method2=[]
	for S in stats:
		if 'p' in S:
			stat_method1.append(S) # Need to go throw node by node
		else:
			stat_method2.append(S) # Go throw each file one by one


	if stat_method1!=[]:
		ncout=process_method1(ncin,all_files,params,stat_method1,nnodes,part)
	if stat_method2!=[]:
		ncout=process_method2(ncin,all_files,params,stat_method2,nnodes)

	
	ncout.close()

            
    
    
if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(prog='extract_stats.py', usage='%(prog)s fileout dirout INDstart INDend params stats (i.e python extract_stats.py stats.nc /home/remy/Buisness/Ciberon/ 15 16 ''{"elev":"elev","temp":"temp"}'' mean min max p90')
	## main arguments
	
	parser.add_argument('fileout', type=str,help='name of the output file (without the extension)')
	parser.add_argument('dirout', type=str,help='name of the output where are the SCHISM files')
	parser.add_argument('params', type=str,nargs='+',help='name of the parameter to plot')
	parser.add_argument('-prefix', type=str,help='prefix default:schout_',default='schout')
	parser.add_argument('-stats', type=str,nargs='+',help='stats to extract',default=['mean','max'])
	parser.add_argument('-start', type=int,help='First file to take',default=1)
	parser.add_argument('-end', type=int,help='Last file to take',default=10000)
	parser.add_argument('-part', type=int,help='Last file to take',default=100)
	args = parser.parse_args()
	

	### PRINT IT ALL
	print 'output name : %s' % (args.fileout)
	print 'Direcory : %s' % (args.dirout)
	print 'From file #%i and #%i' % (args.start,args.end)
	print 'Do parameters : %s' % (args.params)
	
	
	process(args.fileout,args.dirout,args.prefix,args.start,args.end,args.params,args.stats,args.part)
