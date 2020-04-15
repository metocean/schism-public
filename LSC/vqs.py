import numpy as np


def sigma_rutgers_mat(H,KB,sc_w,Cs_w,S_RUTGERS,rutgers):
	VSTRETCHING = rutgers['Vstreching']
	rtheta_s = rutgers['rtheta_s']
	rtheta_b = rutgers['rtheta_b']
	TCLINE  = rutgers['Tcline']

	KBm1 = KB-1
	hc = TCLINE
	z_w=np.ones((S_RUTGERS.shape[1],KB))*0.
	z_w[:,0]=-1.
	Zt_avg1 = 0.

	for k in range(1,KBm1):
		cff_w  = hc*sc_w[k]
		cff1_w = Cs_w[k]
		hinv=1./(hc+H)
		z_w[:,k]=(cff_w+cff1_w*H)*hinv


	S_RUTGERS=np.fliplr(z_w).T

	return S_RUTGERS

def sigma_rutgers_vec(h,KB,sc_w,Cs_w,V_RUTGERS,rutgers):
	VSTRETCHING = rutgers['Vstreching']
	rtheta_s = rutgers['rtheta_s']
	rtheta_b = rutgers['rtheta_b']
	TCLINE  = rutgers['Tcline']

	KBm1 = KB-1
	hc = TCLINE
# ! ------------------------------------------------------------------
# ! Section 2 : Compute Vertical Height as in ROMS RUTGERS/UCLA
# !-----------------------------------------------------------------------
# !  New formulation: Compute vertical depths (meters, negative) at
# !                   RHO- and W-points, and vertical grid thicknesses.
# !  Various stretching functions are possible, as defined above.
# !
# !         z_w(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_w
# !
# !         Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
# !
# !         but with zeta = 0
# !-----------------------------------------------------------------------
	z_w = 0.
	Zt_avg1 = 0.
	
	z_w=np.ones((KB,1))*0.
	z_w[0] = -1.

	for k in range(1,KBm1):
		cff_w  = hc*sc_w[k]
		cff1_w = Cs_w[k]
		hwater=h
		hinv=1./(hc+hwater)
		cff2_w=(cff_w+cff1_w*hwater)*hinv
		z_w[k]=cff2_w


# ! ------------------------------------------------------------------
# ! Section 3 : WRAPPER : ROMS vert. coord. to SCHISM vert. coord.
# !-----------------------------------------------------------------------

	V_RUTGERS = np.flipud(z_w)


	return V_RUTGERS

def sigma_rutgers(KB,sc_w,Cs_w,rutgers):
	VSTRETCHING = rutgers['Vstreching']
	rtheta_s = rutgers['rtheta_s']
	rtheta_b = rutgers['rtheta_b']
	TCLINE  = rutgers['Tcline']

# !
# ! JEROME (IRD Noumea, 5-Oct-2017)
# !
# ! Support to use several Vertical Stretching function for vertical
# ! terrain-following coordinates as implemented in ROMS RUTGERS/UCLA
# ! ------------------------------------------------------------------
# ! Various possible vertical stretching are provided. They are tunable
# ! by using different value for Vstretching (2, 3 or 4), rtheta_s, rtheta_b
# ! and Hc (Tcline). The original vertical stretching function from Song
# ! and
# ! Haidvogel (1994) is not supplied, but can be approached by setting
# ! Vstretching = 2
# !
# ! See: Shchepetkin, A.F. and J.C. McWilliams, 2005: The regional oceanic
# ! modeling system (ROMS): a split-explicit, free-surface,
# ! topography-following-coordinate oceanic model, Ocean
# ! See details: https://www.myroms.org/wiki/Vertical_S-coordinate

	
#     Wrapper :
	KBm1 = KB-1
	hc = TCLINE

	sc_w=np.ones((KB,1))*-1.
	Cs_w=np.ones((KB,1))*-1.

	if VSTRETCHING==2:
# !-----------------------------------------------------------------------
# ! Vstretching = 2 : The New  A. Shchepetkin new vertical stretching
# !                   function.
# ! See :
# !    Shchepetkin, A.F. and J.C. McWilliams, 2005: The regional oceanic !
# !         modeling system (ROMS): a split-explicit, free-surface,      !
# !         topography-following-coordinate oceanic model, Ocean         !
# !         Modelling, 9, 347-404.
# !-----------------------------------------------------------------------
 
		Aweight = 1.
		Bweight = 1.
		ds=1./np.float(KBm1)
		# W-Point layer Interface

		for k in range(KBm1,0,-1):
			cff_w  = ds*np.float(k-KBm1)
			sc_w[k]=cff_w

			if rtheta_s>0:
				Csur=(1.-np.cosh(rtheta_s*cff_w))/(np.cosh(rtheta_s)-1.)

				if rtheta_b>0:
					Cbot=np.sinh(rtheta_b*(cff_w+1.))/np.sinh(rtheta_b)-1.0
					Cweight=(cff_w+1.)**Aweight*(1.+(Aweight/Bweight)* \
						(1.-(cff_w+1.)**Bweight))

					Cs_w[k]=Cweight*Csur+(1.-Cweight)*Cbot
				else:
					Cs_w[k]=Csur

			else:
				Cs_w[k]=cff_w

		
		
		sc_w[0]=-1.
		Cs_w[0]=-1.

# !-----------------------------------------------------------------------
# ! Vstretching = 3 : R. Geyer stretching function for high bottom
# !                   boundary layer resolution
# !-----------------------------------------------------------------------
	if VSTRETCHING==3:

		exp_sur=rtheta_s
		exp_bot=rtheta_b
		Hscale=3.
		ds=1./np.float(KBm1)
# W-Point layer Interface
		sc_w[-1]=0.
		Cs_w[-1]=0.

		for k in range(KBm1-1,0,-1):
			cff_w  = ds*np.float(k-KBm1)
			sc_w[k]= cff_w
			Cbot= np.log(np.cosh(Hscale*(cff_w+1.)**exp_bot))/ \
			   np.log(np.cosh(Hscale))-1.

			Csur=-np.log(np.cosh(Hscale*np.abs(cff_w)**exp_sur))/ \
			   np.log(np.cosh(Hscale))

			Cweight=0.5*(1.-np.tanh(Hscale*(cff_w+0.5)))

			Cs_w[k]=Cweight*Cbot+(1.-Cweight)*Csur

       
		sc_w[0]=-1.
		Cs_w[0]=-1.

# !-----------------------------------------------------------------------
# ! Vstretching = 4 : A. Shchepetkin (UCLA-ROMS, 2010) double vertical
# !                   stretching function
# !-----------------------------------------------------------------------

	if VSTRETCHING==4:
 
		ds=1./np.float(KBm1)
# W-Point layer Interface
		sc_w[-1]=0.
		Cs_w[-1]=0.

		for k in range(KBm1-1,0,-1):
			cff_w = ds*np.float(k-KBm1)
			sc_w[k]  = cff_w

			if rtheta_s> 0.:
				Csur=(1.-np.cosh(rtheta_s*cff_w))/(np.cosh(rtheta_s)-1.)
			else:
				Csur = (cff_w**2)*-1.

			if rtheta_b > 0:
				Cbot=(np.exp(rtheta_b*Csur)-1.)/(1.-np.exp(-rtheta_b))
				Cs_w[k]=Cbot
			else:
				Cs_w[k]=Csur
        
       
		sc_w[0]=-1.
		Cs_w[0]=-1.


	return sc_w,Cs_w 

def quadratic(m,sigma,etal,hsm,a_vqs0):
#         Option 1: quadratic 
	a_vqs=np.max([-1.,a_vqs0-m*0.03])
	tmp=a_vqs*sigma*sigma+(1+a_vqs)*sigma #transformed sigma
	z_mas=tmp*(etal+hsm)+etal
	return z_mas

def sigmaS(m,sigma,hsm,theta_b,theta_f,etal):
#         Option 2: sigma
	#theta_b=0
	#theta_f=1+0.01*(m-1)
	cs=(1-theta_b)*np.sinh(theta_f*sigma)/np.sinh(theta_f)+ \
	theta_b*(np.tanh(theta_f*(sigma+0.5))-np.tanh(theta_f*0.5))/2/np.tanh(theta_f*0.5)
	z_mas=etal*(1+sigma)+hsm[0]*sigma+(hsm[m]-hsm[0])*cs


	return z_mas
class VQS():
	""" A class that manages VQS grid
	"""
	def __init__(self,maxdepth=[],hsm=[],nv_vqs=[]):
		self.m_vqs=10
		self.nv_vqs=4
		if maxdepth == []:
			maxdepth=100
		if hsm==[]:
			self.hsm=np.linspace(2,maxdepth,self.m_vqs,dtype=int)
		else:
			self.hsm=hsm

		
		self.m_vqs=len(self.hsm)
		#self.hsm=[2,6,9,13,17,20,24,28,32,36]
		self.master=[]
		self.sigma_vqs=[]
		self.dp=[]
		if nv_vqs==[]:
			self.nv_vqs=np.ones((self.m_vqs),dtype=int)
			for i,h in enumerate(self.hsm):
				self.nv_vqs[i]=4+1*i # of levels for each master grid (increasing with depth)
		else:
			self.nv_vqs=nv_vqs

	def get_master(self,a_vqs0=-0.3,etal=0,opt=1,theta_b=0,theta_f=1,rutgers=None):

	
		nv_vqs=self.nv_vqs
		nvrt_m=nv_vqs[self.m_vqs-1]

		z_mas=np.ones((int(nvrt_m),int(self.m_vqs)))*-1.e5
		if type(theta_b)!=type(np.array([])):
			theta_b=np.array([theta_b])
		if type(theta_f)!=type(np.array([])):
			theta_f=np.array([theta_f])
		
		if theta_b.shape[0]==1:
			theta_b=np.ones((self.m_vqs,1))*theta_b[0]

		if theta_f.shape[0]==1:
			theta_f=np.ones((self.m_vqs,1))*theta_f[0]

		for m in range(0,self.m_vqs):
			for k in range(0, nv_vqs[m]):
		  		sigma=(k)/(1.0-nv_vqs[m])
		  		if opt==1:
	  				z_mas[k,m]=quadratic(m,sigma,etal,self.hsm[m],a_vqs0)
		  		elif opt>1:
		  			z_mas[k,m]=sigmaS(m,sigma,self.hsm,theta_b[m],theta_f[m],etal)
			
			if opt==3:
				sc_w,Cs_w=sigma_rutgers(nv_vqs[m],0.,0.,rutgers)
				V_RUTGERS=sigma_rutgers_vec(self.hsm[m],nv_vqs[m],sc_w,Cs_w,np.ones((nv_vqs[m],1)),rutgers)
				z_mas[0:np.int(nv_vqs[m]),m] = V_RUTGERS[0:np.int(nv_vqs[m])].T*self.hsm[m]


		z_m=np.ones((int(1+nvrt_m),int(self.m_vqs)))
		z_m[0,:]=nv_vqs[:]
		z_m[1:,:]=z_mas
		self.master=z_m

	def compute_zcor(self,dp,a_vqs0=-0.3,dz_bot_min=.1,opt=1,rutgers=None):


		nodes=len(dp)
		self.dp=dp
		znd=np.ones((self.master.shape[0],nodes))*-dp
		self.sigma_vqs=np.ones((self.master.shape[0],nodes))*99.
		kbp=np.ones((nodes,1),dtype=np.int) 
		nv_vqs=self.master[0,:]

		if opt==3:
			H=np.minimum(dp,self.hsm[0])
			sc_w,Cs_w=sigma_rutgers( np.int(nv_vqs[0]),0.,0. ,rutgers)
			S_RUTGERS=sigma_rutgers_mat(H,np.int(nv_vqs[0]),sc_w,Cs_w,np.ones((np.int(nv_vqs[0]),nodes)),rutgers)

		shallow=dp<=self.hsm[0]
		deep=~shallow

		if np.any(shallow):
			kbp[shallow]=nv_vqs[0]
			for k in range(0,int(nv_vqs[0])):
				if opt==3:
					self.sigma_vqs[k,shallow]=S_RUTGERS[k,shallow]
				else:
					sigma=(k)/(1.0-nv_vqs[0])
					self.sigma_vqs[k,shallow]=a_vqs0*sigma*sigma+(1+a_vqs0)*sigma #transformed sigma
				
				znd[k,shallow]=self.sigma_vqs[k,shallow]*dp[shallow]
		

		zrat=np.ones((nodes,1))*0;
		kbp1=kbp.copy()
		if np.any(deep):
			for m in range(1,self.m_vqs):
				idx=(dp>self.hsm[m-1]) & (dp<=self.hsm[m]) & (deep)
				idxN=idx.nonzero()[0]
				if np.any(idx):
					zrat[idx,0]=((dp[idx])-self.hsm[m-1])/(self.hsm[m]-self.hsm[m-1])
					z1=self.master[1:int(nv_vqs[m])+1,m-1]
					z2=self.master[1:int(nv_vqs[m])+1,m]
					z1=np.tile(z1,(len(idxN),1))
					z2=np.tile(z2,(len(idxN),1))
					Zrat=np.tile(zrat[idx,0],(z1.shape[1],1)).T
					z3=z1+(z2-z1)*Zrat
					minD=np.tile(-dp[idx]+dz_bot_min,(z1.shape[1],1)).T
					D=np.tile(-dp[idx],(z1.shape[1],1)).T
					idx2=z3<minD
					z3[idx2]=D[idx2]
					znd[0:z3.shape[1],idx]=z3.T
					kbp[idxN,0]=np.array([idx2[X,:].nonzero()[0][0] for X in range(idx2.shape[0])])+1
					#kbp1[idxN[idx3[0]],0]=idx3[1]+1
					#import pdb;pdb.set_trace()
					# for ii in idxN:
						# z1=self.master[1:int(nv_vqs[m])+1,m-1]
						# z2=self.master[1:int(nv_vqs[m])+1,m]	
						# z3=z1+(z2-z1)*zrat[ii,0]	
						# bad=np.isnan(z1)
						# z1=z1[~bad]
						# z2=z2[~bad]
						# z3=z3[~bad]
						# idx2=(z3>=-dp[ii]+dz_bot_min)
						# znd[idx2.nonzero()[0],ii]=z3[idx2]
						# k=(~idx2).nonzero()[0][0]
						# kbp[ii,0]=k +1
						# znd[k:nv_vqs[-1],ii]=-dp[ii]

		
		# print znd[:,287]
		# kbp[287]=4
		#kbp[289]=5
		# do sigma_vqs
		nvrt=self.sigma_vqs.shape[0]-1
		self.sigma_vqs[0,deep]=0
		unique_kbp=np.unique(kbp)
		Kbp=kbp[deep]
		deep=deep.nonzero()[0]

		for k in range(1,nvrt):
			idx2=(Kbp>k).nonzero()[0]
			self.sigma_vqs[k-1,deep[idx2]]=znd[k-1,deep[idx2]]/self.dp[deep[idx2]]
			
		#import pdb;pdb.set_trace()
		for k in unique_kbp:
			idx2=(Kbp==k).nonzero()[0]
			self.sigma_vqs[k-1,deep[idx2]]=-1

		
		self.znd=znd
		self.kbp=kbp

	def update_vgrid(self,ind,nlev):
		nvrt=int(self.sigma_vqs.shape[0])
		kbp1=nvrt-nlev-1
		sigma_vqs=np.flipud(self.sigma_vqs)

		for ii in ind:
			sigout=np.ones(nvrt)*99.
			sig0=np.ones(nvrt)*0
			kbp0=int(nvrt-self.kbp[ii])
			sigout[-nlev]=-1.
			sigout[-1]=0.
			#import ipdb;ipdb.set_trace()
			sig0[kbp0:nvrt]=np.linspace(int(-1),int(0),int(nvrt-kbp0))
			sig1=np.linspace(int(-1),int(0),int(nlev))
			for r in range(1,nlev-1):
				for kk in range( int(nvrt-self.kbp[ii]),int(nvrt)):
					if (sig1[r]>=sig0[kk-1]) & (sig1[r]<=sig0[kk]):
						rat=(sig1[r]-sig0[kk-1])/(sig0[kk]-sig0[kk-1])
						sigout[kbp1+r+1]=sigma_vqs[kk,ii]*rat+sigma_vqs[kk-1,ii]*(1-rat)
						break

			self.sigma_vqs[:,ii]=np.flipud(sigout)
			self.kbp[ii]=nvrt-kbp1-1
			znd=np.ones(nvrt)*-1e-6
			znd[0]=0
			znd[0:nlev]=self.dp[ii]*np.flipud(sigout[-nlev:])
			self.znd[:,ii]=znd


		

		#znd[k:nv_vqs[-1]-1,ind]=-self.dp[ind]
		


	def extract_trs(self,x,y,nodes):
		Ltrs=[0]
		Ltrs.extend(np.sqrt(np.diff(x)**2+np.diff(y)**2))
		self.Ltrs=np.cumsum(Ltrs)
		self.kbptrs=self.kbp[nodes]
		self.zndtrs=self.znd[:,nodes]
		
	def export_vgrid(self,fileout):
#     Output in SELFE convention
		a=open(fileout,'w')
		a.write('%12i \n' % 1)
		nvrt=self.sigma_vqs.shape[0]-1
		a.write('%12i \n' %(nvrt))
		sigma_vqs=self.sigma_vqs.copy()

		sigma_vqs=np.flipud(sigma_vqs)
		sigma_vqs[sigma_vqs==0]=0
		tout=np.vstack((np.arange(1,len(self.dp)+1).flatten(),\
			nvrt-self.kbp.flatten()+1,\
			sigma_vqs)).T
		
		a.close()
		with open(fileout,'a') as fh:
			for row in tout:
				line='%11i%11i' % (row[0],row[1])
				row=row[2:]
				lin="".join('%15.6f' % value for value in row[row!=99])
				fh.write(line+lin + '\n')

		fh.close()
		
