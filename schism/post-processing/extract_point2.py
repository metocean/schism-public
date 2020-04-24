#!/usr/bin/env python2.7

import os,sys
import netCDF4
import numpy
import scipy.interpolate as spint

# import scipy.spatial.qhull as qhull
import matplotlib.tri as mtri
import json



def interp3D(zdata,udata,wanted_depth):
	m,n=zdata.shape
	zz=zdata.flatten(1)
	yinterp=numpy.zeros(shape=(n,1))
	yv=udata.flatten(1)
	for i in range(0,n):
		if((wanted_depth[i] < zz[((i*m)+0)]) or (wanted_depth[i] > zz[(i*m)+(m-1)])):
			yinterp[i] = numpy.NaN
		else:
			for j in range(0,m):
				if(wanted_depth[i]<=zz[(i*m)+j]):
					yinterp[i] = (wanted_depth[i]-zz[(i*m)+(j-1)]) / (zz[(i*m)+j]-zz[(i*m)+(j-1)]) * (yv[(i*m)+j]-yv[(i*m)+(j-1)]) + yv[(i*m)+(j-1)]; 
					           									
					break
				else:
					yinterp[i] = numpy.NaN
	return yinterp

def create_output_file(filout,params,depth,dtype):
	f = file(filout, "w")
	# header
	f.write('Year\tMonth\tDay\tHour\tMin\tSec\tDepth\t')
	col_tot=0
	val=-1

	for vars in params:      
		val=val+1
		if dtype[val][0]==2: # 3D variable 
			if dtype[val][1] ==2: # U and V data		
				for n in range(0,len(depth)):
					f.write('%s\t' % ('u_'+vars+'lev_'+str(depth[n])))
					col_tot=col_tot+1

				for n in range(0,len(depth)):
					f.write('%s\t' % ('v_'+vars+'_lev_'+str(depth[n])))
					col_tot=col_tot+1
			else:
				for n in range(0,len(depth)):
					f.write('%s\t' % (vars+'_lev_'+str(depth[n])))
					col_tot=col_tot+1

		else: # 2D variable
			if dtype[val][1] ==2: # U and V data
				f.write('%s\t%s\t' % ('u_'+vars,'v_'+vars))
				col_tot=col_tot+2
			else:
				f.write('%s\t' % (vars))
				col_tot=col_tot+1

	f.write('\n')
	return f,col_tot
 
def pair(arg):
    # For simplity, assume arg is a pair of integers
    # separated by a comma. If you want to do more
    # validation, raise argparse.ArgumentError if you
    # encounter a problem.
    return [float(x) for x in arg.split(',')]

def init_nc(dirin,pref,Istart):


	## extract static vars
	ncfile=os.path.join(dirin,pref+'_'+str(Istart)+'.nc')
	ncs=netCDF4.Dataset(ncfile)
	X=ncs.variables['SCHISM_hgrid_node_x'][:]
	Y=ncs.variables['SCHISM_hgrid_node_y'][:]

	ele=ncs.variables['SCHISM_hgrid_face_nodes'][...,:3]-1
	triang = mtri.Triangulation(X, Y, ele)
	nt=len(ncs.variables['time'])
	ncs.close()

	return nt,X,Y,triang

def get_depth(x,y,triang,ncs):

    Z=ncs.variables['depth'][:]
    interp_lin = mtri.LinearTriInterpolator(triang, Z)
    return interp_lin(x,y).data

def get_closest_node(X,Y,xi,yi):

    Z=numpy.arange(0,len(X))
    nodes = spint.griddata((X,Y),Z,numpy.vstack((xi,yi)).T, method='nearest')
    return nodes

def process(POS,fileout,dirin,Istart,Iend,params,depth,prefix):

	nt,X,Y,triang=init_nc(dirin,prefix,Istart)
	params_str=','.join(params)

	fout=[]
	nodes=[]
	for k in range(0,len(POS[0]),2):
		x=POS[0][k]
		y=POS[0][k+1]
		# create the output file (one for each point)
		fout.append(fileout+'_'+str(x)+'_'+str(y)+'.txt')

		nodes.append(get_closest_node(X,Y,x,y))
	


#### For each of the file
	for nfile in range(Istart,Iend+1):
		fileIN=os.path.join(dirin,prefix+'_'+str(nfile)+'.nc')
		if not os.path.isfile(fileIN):
			continue

		print ' read file: '+str(nfile)
		for i,node in enumerate(nodes):
			ncs=netCDF4.Dataset(fileIN)
			Total= numpy.zeros(shape=(nt,7))
			dtime = netCDF4.num2date(ncs.variables['time'][:],ncs.variables['time'].units)
			Total[:,0]=[dt.year for dt in dtime.astype(object)]
			Total[:,1]=[dt.month for dt in dtime.astype(object)]
			Total[:,2]=[dt.day for dt in dtime.astype(object)]
			Total[:,3]=[dt.hour for dt in dtime.astype(object)]
			Total[:,4]=[dt.minute for dt in dtime.astype(object)]
			Total[:,5]=[dt.second for dt in dtime.astype(object)]
			Total[:,6]=ncs.variables['depth'][node]
			if nfile==Istart:
				f = file(fout[i], "w")
				f.write('Year\tMonth\tDay\tHour\tMin\tSec\tDepth\t')


			for j,param in enumerate(params):
				data=ncs.variables[param][:,node]
				data = numpy.ma.masked_where(data==-9999,data)
				i23d=ncs.variables[param].i23d
				ivs=ncs.variables[param].ivs

				if ivs==1 and i23d==1: # elev
					Total = numpy.hstack((Total, numpy.zeros((Total.shape[0], 1), dtype=Total.dtype)))
					Total[:,-1]=data[:,0]
					if nfile==Istart:
						f.write('%s\t' % param)

				if ivs==2 and i23d==1: # dahv
					Total = numpy.hstack((Total, numpy.zeros((Total.shape[0], 2), dtype=Total.dtype)))
					Total[:,-2:]=numpy.squeeze(data[:,0,:])
					if nfile==Istart:
						f.write('%s_U\t' % param)
						f.write('%s_V\t' % param)	

				if ivs==1 and i23d==2: # zcor
					data=numpy.squeeze(data[:,0])
					zdata=ncs.variables['zcor'][:,node]
					zdata=numpy.squeeze(zdata[:,0,:])
					zdata = numpy.ma.masked_where(zdata==-9999,zdata)
					if depth is None:
						if nfile==Istart:
							f.write('%s_dahv\t' % (param))
						dz=numpy.diff(zdata)
						dep=numpy.max(zdata,1)-numpy.min(zdata,1)
						g=(data[:,0:-1]+data[:,1:])/2
						Total = numpy.hstack((Total, numpy.zeros((Total.shape[0], 1), dtype=Total.dtype)))
						Total[:,-1]=numpy.sum(g*abs(dz),1)/dep


					else:
						for k in range(0,len(depth)):
							if nfile==Istart:
								f.write('%s_lev_%.1f\t' % (param,depth[k]))
							if depth[k]>0: # above sea bed
								wanted_depth=numpy.min(zdata,1)+depth[k];
							if depth[k]<=0: # below sea surface
								wanted_depth=numpy.max(zdata,1)+depth[k];

							data3d=interp3D(zdata.T,data.T,wanted_depth.T)
							Total = numpy.hstack((Total, numpy.zeros((Total.shape[0], 1), dtype=Total.dtype)))
							Total[:,-1]=data3d[:,0]
					

				if ivs==2 and i23d==2: # hvel
					data=numpy.squeeze(data[:,0])
					zdata=ncs.variables['zcor'][:,node]
					zdata=numpy.squeeze(zdata[:,0,:])
					zdata = numpy.ma.masked_where(zdata==-9999,zdata)
					if depth is None:
						print "enter a depth"
						sys.exit(-1)

					for k in range(0,len(depth)):
						if nfile==Istart:
							f.write('%s_U_lev_%.1f\t' % (param,depth[k]))
							f.write('%s_V_lev_%.1f\t' % (param,depth[k]))
						if depth[k]>0: # above sea bed
							wanted_depth=numpy.min(zdata,1)+depth[k];
						if depth[k]<=0: # below sea surface
							wanted_depth=numpy.max(zdata,1)+depth[k];

						data3dU=interp3D(zdata.T,data[:,:,0].T,wanted_depth.T)
						data3dV=interp3D(zdata.T,data[:,:,1].T,wanted_depth.T)
						Total = numpy.hstack((Total, numpy.zeros((Total.shape[0], 1), dtype=Total.dtype)))
						Total[:,-1]=data3dU[:,0]
						Total = numpy.hstack((Total, numpy.zeros((Total.shape[0], 1), dtype=Total.dtype)))
						Total[:,-1]=data3dV[:,0]

			if nfile==Istart:
				f.write('\n')
			numpy.savetxt(f,Total,fmt='%g',delimiter='\t')
		ncs.close()
 	f.close()


     
  
    
    
if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(prog='extract_point.py', usage='%(prog)s fileout dirout INDstart INDend params (i.e python extract_point.py caca /home/remy/Buisness/Ciberon/ 15 16 ''{"elev":"elev","temp":"temp"}'' -POS 256842,9291527 -depth -5 -4')
	## main arguments
	
	parser.add_argument('fileout', type=str,help='name of the output file (without the extension)')
	parser.add_argument('dirout', type=str,help='name of the output where are the SELFE files')
	parser.add_argument('INDstart', type=int,help='First file to take')
	parser.add_argument('INDend', type=int,help='Last file to take')
	parser.add_argument('params', type=str,nargs='+',help='name of the parameter to plot')
	parser.add_argument('pos', type=pair,nargs='+',help='XY position to extract')
	parser.add_argument('-prefix', type=str,help='prefix default:schout_',default='schout')
	parser.add_argument('-depth', type=float,nargs='+',help='Z position to extract must be negative for bsl or pos for asb',default=None)
	args = parser.parse_args()
	

	### PRINT IT ALL
	print 'output name : %s' % (args.fileout)
	print 'Direcory : %s' % (args.dirout)
	print 'From file #%i and #%i' % (args.INDstart,args.INDend)
	print 'Do parameters : %s' % (args.params)

	
	
	process(args.pos,args.fileout,args.dirout,args.INDstart,args.INDend,args.params,args.depth,args.prefix)

