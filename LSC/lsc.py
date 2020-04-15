import os
import sys
import yaml

import numpy as np
import numexpr as ne
from scipy.interpolate import griddata
import matplotlib.collections
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from PyQt5 import QtCore, QtGui, uic
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget,QMessageBox,QInputDialog,QFileDialog

from pyschism.mesh import Hgrid
import vqs
import copy

root_path = os.path.dirname(os.path.abspath(__file__))
qtCreatorFile = os.path.join(root_path,"lsc.ui") 

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

global Xtrs,Ytrs,node
Xtrs=[]
Ytrs=[]
node=[]
global view_type
view_type='bathy'
class SelectFromCollection(object):
    """Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool highlights
    selected points by fading them out (i.e., reducing their alpha values).
    If your collection has alpha < 1, this tool will permanently alter them.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        Axes to interact with.

    collection : :class:`matplotlib.collections.Collection` subclass
        Collection you want to select from.

    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to `alpha_other`.
    """

    def __init__(self, ax, collection,parent):
        self.ax=ax
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.xys = collection
        self.Npts = len(self.xys)
        self.parent=parent

        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        if self.parent._active == "ZOOM" or self.parent._active == "PAN":
                return
        path = Path(verts)

        self.ind .extend( np.nonzero([path.contains_point(xy) for xy in self.xys])[0])
        self.ind=list(np.unique(self.ind))
        for line in self.ax.axes.lines:
            if line.get_gid() == 'node':
                    line.remove()
        self.ax.plot(self.xys[self.ind,0],self.xys[self.ind,1],'r.',gid='node')
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.canvas.draw_idle()



class MyToolbar(NavigationToolbar):
    def __init__(self, figure_canvas,axes_canvas, parent= None):
        self.toolitems = (('Home', 'Reset', 'home', 'home'),
            ('Pan', 'Pan', 'move', 'pan'),
            ('Zoom', 'Zoom', 'zoom_to_rect', 'zoom'),
            (None, None, None, None),
            ('Transect', 'make a transect',os.path.join(root_path,'icons','transect'), 'transect_tool'),
            ('View', 'Switch view',os.path.join(root_path,'icons','switch'), 'switch_view'),
            ('Selector', 'Select node', os.path.join(root_path,'icons','select'), 'select_tool'),
            (None, None, None, None),
            ('Save', 'Save', 'filesave', 'save_figure'))


        #('Subplots', 'putamus parum claram', 'subplots', 'configure_subplots'),
        self.figure_canvas=figure_canvas
        self.axes_canvas=axes_canvas
    

        NavigationToolbar.__init__(self, figure_canvas, parent= parent)
        #NavigationToolbar.setWindowIcon(QIcon('a.ico'))



    def remove_series(self,id='transect'):
        plt=self.axes_canvas.get_axes()[0]
        for c in plt.lines:
            if c.get_gid() == id:
                    c.remove()

    def onclick(self,event):
        self.remove_series()
        if event.dblclick:
            self.figure_canvas.mpl_disconnect(self.cid)
            del self.cid
            app=self.parentWidget().parent().parent()
            app.create_transect(x0=np.array(Xtrs),y0=np.array(Ytrs))
            self.remove_series()
            ax=self.axes_canvas.get_axes()[0]
            ax.plot(app.Xtrs,app.Ytrs,'r.-',gid='transect')
            app.vgrid.extract_trs(app.Xtrs,app.Ytrs,app.Ntrs)
            app.draw_map()
            app.draw_vgrid()
    

        else:
            
            Xtrs.append(event.xdata)
            Ytrs.append(event.ydata)
            ax=self.axes_canvas.get_axes()[0]
            ax.plot(Xtrs,Ytrs,'r.-',gid='transect')    


        self.figure_canvas.draw_idle()

        return Xtrs,Ytrs




    def transect_tool(self):

        if hasattr(self,'cid'):
            self.figure_canvas.mpl_disconnect(self.cid)
            del self.cid
        else:
            global Xtrs,Ytrs
            Xtrs=[]
            Ytrs=[]
            self.cid = self.figure_canvas.mpl_connect('button_press_event', self.onclick)

    def select_tool(self):
        
        app=self.parentWidget().parent().parent()
        if hasattr(self,'selector'):
            self.selector.disconnect()
            num,ok = QInputDialog.getInt(self,"Number of layer","enter a number")
            self.selector.disconnect()
            self.remove_series('node')
            ind=self.selector.ind
            del self.selector

            self.figure_canvas.draw_idle()
            app.vgrid.update_vgrid(ind,num)
            app.nlev.values[:]=app.vgrid.kbp[:,0]
            app.vgrid.extract_trs(app.Xtrs,app.Ytrs,app.Ntrs)
            app.draw_map()
            app.draw_vgrid()

        else:
            ax=self.axes_canvas.get_axes()[0]
            Y=app.gr.y
            X=app.gr.x
            self.selector = SelectFromCollection(ax, np.vstack((X.flatten(),Y.flatten())).T,self)    
            self.figure_canvas.draw_idle()


    def switch_view(self):
        global view_type
        if view_type=='bathy':
            view_type='lev'
        else:
            view_type='bathy'
        
        app=self.parent
        app=self.parentWidget().parent().parent()
        app.draw_map()



class MyApp(QMainWindow, Ui_MainWindow):

    def __init__(self,app):
        QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)

        self.setupUi(self)
        self.bot_layer_thickness.valueChanged.connect(self.refresh_it_all)
        self.avqs0.valueChanged.connect(self.refresh_it_all)
        self.sigma_type.currentChanged.connect(self.refresh_it_all)
        self.thetab.valueChanged.connect(self.refresh_it_all)
        self.thetaf.valueChanged.connect(self.refresh_it_all)
        self.rtheta_b.valueChanged.connect(self.refresh_it_all)
        self.rtheta_s.valueChanged.connect(self.refresh_it_all)
        self.tcline.valueChanged.connect(self.refresh_it_all)
        self.vstreching.activated.connect(self.refresh_it_all)
        self.rutgers.toggled.connect(self.refresh_it_all)
        self.hsm.returnPressed.connect(self.refresh_it_all)
        self.nv_vqs.returnPressed.connect(self.refresh_it_all)
        self.thetaf_exp.returnPressed.connect(self.refresh_it_all)
        self.thetab_exp.returnPressed.connect(self.refresh_it_all)
        self.import_vgrid.triggered.connect(self.importV)
        self.import_nlev.triggered.connect(self.importLEV)
        self.import_hgrid.triggered.connect(self.importH)
        self.export_vgrid.triggered.connect(self.exportV)
        self.export_vgrid_params.triggered.connect(self.exportVparams)
        self.export_trs.triggered.connect(self.exportTRS)
        self.export_nlev.triggered.connect(self.exportLEV)
        self.Xtrs=[]
        self.Ytrs=[]
        self.Ztrs=[]
        self.Ntrs=[]
        self.nlev=[]

        self.gr = Hgrid.open(os.path.expanduser(app.arguments()[1]))
        self.gr.values[:]=self.gr.values[:]*-1.
        self.nlev=copy.deepcopy(self.gr)
        self.nlev.values[:]=self.nlev.values[:]*0

        self.create_transect()
        #self.hsm.setText('0:%i:50'% (np.ceil(self.gr.mesh.nodes[:,2].max())+50))
        self.hsm.setText('2 12 22 32 42 52 62 72 82 200 2000')
        self.vgrid=self.create_vgrid(maxdepth=np.ceil(self.gr.values.max()),hsm=self.get_hsm(self.hsm.text()))
        self.vgrid.compute_zcor(self.gr.values,a_vqs0=-0.3,opt=1)
        self.nlev.values[:]=self.vgrid.kbp[:,0]
        self.vgrid.extract_trs(self.Xtrs,self.Ytrs,self.Ntrs)
        self.create_main_frame()
        self.create_vertical_frame()
        self.draw_map()
        self.draw_vgrid()

    def exportTRS(self):
        file_name = QFileDialog.getSaveFileName(self, "Save transect file", "", "Transect (*.bp)")
        Z=self.Ztrs*-1.
        X=self.Xtrs
        Y=self.Ytrs
        with open(file_name,'w') as fh:
            fh.write('%s\n' % ' ')
            fh.write('%i\n' % len(X))
            for n in range(0,len(X)):
                line='%i\t%.2f\t%.2f\t%.2f\n' % (n,X[n],Y[n],Z[n])
                fh.write(line)

        fh.close()
    def exportLEV(self):
        file_name = QFileDialog.getSaveFileName(self, "Number of level", "", "SCHSIM grid file (*.gr3)")
        gr=copy.deepcopy(self.gr)
        gr.values[:]=self.vgrid.kbp[:,0]
        gr.write(str(file_name))

    def importLEV(self):
        file_name = QFileDialog.getOpenFileName(self, "Load nlev.gr3 file", "", "SCHSIM grid file(*.gr3)")
        self.nlev=Hgrid.open(str(file_name))
        self.change_nlev()
        self.draw_map()
        self.draw_vgrid()
    def exportV(self):
        file_name = QFileDialog.getSaveFileName(self, "Save Vgrid file", "", "SCHSIM vgrids (*.in)")
        self.vgrid.export_vgrid(file_name)

    def exportVparams(self):
        file_name = QFileDialog.getSaveFileName(self, "Save Vgrid param file", "", "SCHSIM vgrids params(*.yaml)")
        params=self.get_all_value()
        params['nv_vqs']=params['nv_vqs'].tolist()
        params['hsm']=params['hsm'].tolist()
        if type(params['thetaf'])==np.ndarray:
            params['thetaf']=params['thetaf'].tolist()
        if type(params['thetab'])==np.ndarray:
            params['thetab']=params['thetab'].tolist()
        with open(file_name, 'w') as yaml_file: # this would write the yaml file that my function read probably best so we can track
            yaml.dump(params, yaml_file, default_flow_style=False)


    def importH(self):
        pass
    def importV(self):
        file_name = QFileDialog.getOpenFileName(self, "Load Vgrid param file", "", "SCHSIM vgrids params(*.yaml)")
        with open(file_name ,'r') as f:
            params = yaml.load(f)
        params['nv_vqs']=np.asarray(params['nv_vqs'])
        params['hsm']=np.asarray(params['hsm'])
        self.set_all_value(params)

    def create_transect(self,x0=[],y0=[]):
        x= self.gr.x
        y= self.gr.y 
        if x0==[]: #first transect
            x0=x.min()
            y0=y.min()
            x1=x.max()
            y1=y.max()

            X= np.arange(x0,x1,np.ceil((x1-x0)/100.))
            Y= np.arange(y0,y1,np.ceil((y1-y0)/100.))
        else:
            total_len=0
            for n in range(1,len(x0)):
                total_len=total_len+(np.sqrt((x0[n-1]-x0[n])**2+(y0[n-1]-y0[n])**2))


            X=[x0[0]]
            Y=[y0[0]]
            
            dx=total_len/100.
            for n in range(1,len(x0)):
                sub_len=np.sqrt((x0[n]-x0[n-1])**2+(y0[n]-y0[n-1])**2)
                N=np.floor(sub_len/dx)
                dist_x=(x0[n]-x0[n-1])/N
                dist_y=(y0[n]-y0[n-1])/N


                for I in range(0,int(N)):
                    new_len=np.sqrt((X[-1]+dist_x-x0[n-1])**2+(Y[-1]+dist_y-y0[n-1])**2)
                    if new_len<sub_len:
                        X.append(X[-1]+dist_x)
                        Y.append(Y[-1]+dist_y)
                    


        N=griddata((x,y),np.arange(0,len(self.gr.values),1),(X,Y),method='nearest')
        indexes = np.unique(N, return_index=True)[1]
        N=[N[index] for index in sorted(indexes)]
        gd=(self.gr.values[N]>-1).nonzero()[0]
        
        N= N[gd[0]:gd[-1]]
        self.Ztrs=self.gr.values[N]
        self.Xtrs=self.gr.x[N]
        self.Ytrs=self.gr.y[N]
        self.Ntrs=N

    def create_vgrid(self,maxdepth=[],a_vqs0=-0.3,etal=0,opt=1,theta_b=0,theta_f=1,hsm=[],nv_vqs=[],rutgers=None):
        vgrid=vqs.VQS(maxdepth,hsm=hsm,nv_vqs=nv_vqs)

        vgrid.get_master(a_vqs0,etal,opt,theta_b,theta_f,rutgers)

        return vgrid

    def change_nlev(self):
        old_nlev=self.vgrid.kbp[:,0]
        new_nlev=self.nlev.values

        if sum(new_nlev)>0:
	        di=old_nlev-new_nlev
	        di_lev_unique=set(np.delete(di,np.where(di==0.0),axis=0))
	        for nlev in di_lev_unique:
	            ind=np.where(di==nlev)[0]
	            new=self.vgrid.kbp[ind[0],0]-int(nlev)
	            self.vgrid.update_vgrid(ind,int(new))
	            self.nlev.values[ind]=int(new)


    def create_main_frame(self):
        self.main_frame = QWidget()
        self.fig = Figure((5.0, 4.0), dpi=100)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        self.canvas.setFocusPolicy(Qt.StrongFocus)
        self.canvas.setFocus()

        self.mpl_toolbar = MyToolbar(self.canvas,self.fig, self.main_frame)
        self.mpl_toolbar.children()[15].setCheckable(True)
        self.mpl_toolbar.update()

        self.Vbox.addWidget(self.canvas)  # the matplotlib canvas
        self.Vbox.addWidget(self.mpl_toolbar)


    def create_vertical_frame(self):
        self.vert_frame = QWidget()
        self.vert_fig = Figure((5.0, 4.0), dpi=100)
        self.vert_canvas = FigureCanvas(self.vert_fig)
        self.vert_canvas.setParent(self.vert_frame)
        self.vert_canvas.setFocusPolicy(Qt.StrongFocus)
        self.vert_canvas.setFocus()

        self.mpl_vert_toolbar = NavigationToolbar(self.vert_canvas, self.vert_frame)

        self.Vbox_vert.addWidget(self.vert_canvas)  # the matplotlib canvas
        self.Vbox_vert.addWidget(self.mpl_vert_toolbar)

    def draw_map(self,nval=60,Zmin=[],Zmax=[]):
        gr=copy.deepcopy(self.gr)
        elem=np.array(gr.elements,dtype=int)-1
        

        if hasattr(self,'axes'):
            ylim=self.axes.get_ylim()
            xlim=self.axes.get_xlim()
            self.axes.clear()
            self.axes_cb.clear()
            #pass

        else:
            self.fig.clear()
            self.axes = self.fig.add_axes([0, 0, 1, 1])
            self.axes_cb = self.fig.add_axes([0.85, 0.05, 0.03, 0.9])
            ylim=[]
            self.tri_idx=elem[:,-1]<0
            self.quad=~self.tri_idx
            self.tris=elem[self.tri_idx,0:3]
            self.els=np.ones((self.tris.shape[0]+self.quad.nonzero()[0].shape[0]*2,3),dtype=int)
            self.els[0:self.tris.shape[0],:]=self.tris[:]
            I=self.tris.shape[0]
            for i,q in enumerate(self.quad.nonzero()[0]):
            	self.els[I,:]=elem[q,0:3]
            	self.els[I+1,:]=elem[q,[0,2,3]]
            	I=I+2
            self.triang=Triangulation(gr.x,gr.y,self.els)






        if view_type=='bathy':
            Z=self.gr.values 
            if Zmax==[]:
                Zmax=np.ceil(Z.max())
            if Zmin==[]:
                Zmin=np.ceil(Z.min())
            levels = np.linspace(Zmin, Zmax, nval)
            ticks=np.floor(np.linspace(max(Zmin,0), Zmax, 10))
            tit='Bathymetrie'

        else:

            gr.values[:]= self.vgrid.kbp[:,0]
            Zmax=gr.values.max()+1
            Zmin=gr.values.min()
            nval=Zmax-Zmin+1
            levels = np.linspace(int(Zmin),int(Zmax),int(nval))
            ticks=np.floor(levels)+.5
            tit='Number of Sigma level'

        

        
        #quads=Triangulation(gr.x,gr.y,quads)
        
        tricon=self.axes.tricontourf(self.triang,gr.values,vmin=Zmin,vmax=Zmax,cmap=plt.cm.Spectral_r,levels=levels, origin='lower',antialiased=False)



        tricon.set_clim(Zmin,Zmax)
        self.axes.tick_params(labelbottom='off',labelleft='off')
        self.axes.set_aspect('equal','datalim','C')
        self.axes.plot(self.Xtrs,self.Ytrs,'r-',gid='transect')
    
        self.cb=self.fig.colorbar(tricon,self.axes_cb,ticks=ticks)
        if view_type=='lev':
            self.cb.ax.set_yticklabels(levels)  # horizontal colorbar


        self.titre.setText(tit)
        if ylim!=[]:
            self.axes.set_ylim(ylim)
            self.axes.set_xlim(xlim)

        self.canvas.draw_idle()



    def draw_vgrid(self):
        if hasattr(self,'axes_top'):
            self.axes_top.clear()
            self.axes_bot.clear()

        else:
            self.vert_fig.clear()
            self.axes_top = self.vert_fig.add_subplot(211)
            self.axes_bot = self.vert_fig.add_subplot(212)
            pos1=self.axes_top.get_position()
            pos2=[pos1.x0 , pos1.y0 + 0.05,  pos1.width , pos1.height ]
            self.axes_top.set_position(pos2)

        x=np.arange(1,self.vgrid.master.shape[1]+1);
        zcor_m=self.vgrid.master[1:,]

        zcor_m[zcor_m==-100000.]=np.nan

        self.axes_top.plot(x,zcor_m.T,'b-')
        self.axes_top.plot(x,-self.vgrid.hsm,'ko')

        X=np.tile(x,(self.vgrid.master.shape[0]-1,1))
        self.axes_top.plot(X,zcor_m,'k-')

        self.axes_top.set_title('Master grid')
        self.axes_top.set_xlabel('Grid #')

    
        self.axes_bot.set_title('Transect before adjustment (transect1)')        
        self.axes_bot.set_xlabel('Along transect distance (m)')
        
        self.axes_bot.plot(self.vgrid.Ltrs,self.vgrid.zndtrs.T,'b-')
        self.axes_bot.plot(self.vgrid.Ltrs,-self.Ztrs,'r.')

        L=np.tile(self.vgrid.Ltrs,(self.vgrid.zndtrs.shape[0],1))
        self.axes_bot.plot(L,self.vgrid.zndtrs,'k-')


        
        self.vert_canvas.draw()
    def get_all_value(self):
        params={}
        params['a_vqs0']=self.avqs0.value()
        params['dz_bot_min']=self.bot_layer_thickness.value()
        params['hsm']=self.get_hsm(self.hsm.text())

        if self.thetab_exp.text()=='':
            thetab=self.thetab.value()

        else:
            thetab=self.get_theta(params['hsm'],self.thetab_exp.text())
        if self.thetaf_exp.text()=='':
            thetaf=self.thetaf.value()
        else:
            thetaf=self.get_theta(params['hsm'],self.thetaf_exp.text())

        params['thetab']=thetab
        params['thetaf']=thetaf
        params['nv_vqs']=self.get_nv_vqs(params['hsm'],self.nv_vqs.text())
        opt=self.sigma_type.currentIndex()+1
        if self.rutgers.isChecked() and opt==2:
            opt=3
            params['rtheta_s']=self.rtheta_s.value()
            params['rtheta_b']=self.rtheta_b.value()
            params['Tcline']=self.tcline.value()
            params['Vstreching']=self.vstreching.currentIndex()+2

        opt_label=['quadratic','S-layers','Rutgers']

        params['mode']=opt_label[opt-1]

        return params
    def set_all_value(self,params):
        hsm=[]
        for n in params['hsm']:
            try:
                hsm.append([n[0]])
            except:
                hsm.append([n])

        hsm=str(hsm)
        hsm=hsm[1:-1].replace(',','')
        self.hsm.setText(hsm)

        nv_vqs=[]
        for n in params['nv_vqs']:
            try:
                nv_vqs.append([n[0]])
            except:
                nv_vqs.append([n])


        nv_vqs=str(nv_vqs)
        nv_vqs=nv_vqs[1:-1].replace(',','')
        self.nv_vqs.setText(nv_vqs)


        if type(params['thetaf'])==list:
            thetaf=[]
            for n in params['thetaf']:
                try:
                    thetaf.append([n[0]])
                except:
                    thetaf.append([n])

            thetaf=str(thetaf)
            thetaf=thetaf[1:-1].replace(',','')
            self.thetaf_exp.setText(thetaf)
        else:
            self.thetaf.setValue(params['thetaf'])


        if type(params['thetab'])==list:
            thetab=[]
            for n in params['thetab']:
                try:
                    thetab.append([n[0]])
                except:
                    thetab.append([n])

            thetab=str(thetab)
            thetab=thetab[1:-1].replace(',','')
            self.thetab_exp.setText(thetab)
        else:
            self.thetab.setValue(params['thetab'])



        self.avqs0.setValue(params['a_vqs0'])
        self.bot_layer_thickness.setValue(params['dz_bot_min'])




        if params['mode']=='quadratic':
            self.sigma_type.setCurrentIndex(0)
        else:
            self.sigma_type.setCurrentIndex(1)

        if params['mode']=='Rutgers':
            self.rutgers.setChecked(True)
            self.rtheta_s.setValue(params['rtheta_s'])
            self.rtheta_b.setValue(params['rtheta_b'])
            self.tcline.setValue(params['Tcline'])
            self.vstreching.setCurrentIndex(params['Vstreching']-2)

        self.refresh_it_all()
    def refresh_it_all(self):
        a_vqs0=self.avqs0.value()
        dz_bot_min=self.bot_layer_thickness.value()
        opt=self.sigma_type.currentIndex()+1

        try:
            hsm=self.get_hsm(self.hsm.text())
        except:
            QMessageBox.information(QWidget(), "No", "Syntax not correct for depth" )
            return


        if self.thetab_exp.text()=='':
            thetab=self.thetab.value()
        else:
            thetab=self.get_theta(hsm,self.thetab_exp.text())
        if self.thetaf_exp.text()=='':
            thetaf=self.thetaf.value()
            
            if thetaf<=0:
            	self.thetaf.setValue(.01)
            	thetaf=.1

        else:
            thetaf=self.get_theta(hsm,self.thetaf_exp.text())


        try:
            nv_vqs=self.get_nv_vqs(hsm,self.nv_vqs.text())
        except:
            QMessageBox.information(QWidget(), "No", "Syntax not correct for N lev" )
            return

        maxdepth=np.ceil(self.gr.values.max())
        if hsm.max()<maxdepth:
            QMessageBox.critical(QWidget(), "No", "last depth must be > Max depth of  %.f " % (maxdepth))
            return

        if len(hsm)<2:
            QMessageBox.critical(QWidget(), "No", "You need at least 2 master grid")
            return

        if len(hsm)>100:
            QMessageBox.critical(QWidget(), "No", "Too much")
            return
        
        rutgers={}
        if self.rutgers.isChecked() and opt==2:
            opt=3
            rutgers['rtheta_s']=self.rtheta_s.value()
            rutgers['rtheta_b']=self.rtheta_b.value()
            rutgers['Tcline']=self.tcline.value()
            rutgers['Vstreching']=self.vstreching.currentIndex()+2


        self.vgrid=self.create_vgrid(maxdepth=np.ceil(self.gr.values.max()),\
            a_vqs0=a_vqs0,etal=0,opt=opt,\
            theta_b=thetab,theta_f=thetaf,hsm=hsm,nv_vqs=nv_vqs,rutgers=rutgers)



        self.vgrid.compute_zcor(self.gr.values,a_vqs0=a_vqs0,dz_bot_min=dz_bot_min,opt=opt,rutgers=rutgers)

        self.change_nlev()
            
        self.vgrid.extract_trs(self.Xtrs,self.Ytrs,self.Ntrs)
        self.draw_vgrid()
        self.draw_map()


    def get_theta(self,hsm,pp):
        if type(pp)!=type(str()):
            pp= str(pp)

        if ':' in pp:
            theta=eval('np.r_['+pp+']')
        elif 'N' in pp:
            theta=np.ones((hsm.shape[0],1))
            for N in range(0,len(hsm)):
                theta[N]=eval(pp) # of levels for each master grid (increasing with depth)
        else:
            pp=pp.replace(' ',',')
            theta=np.array(eval(pp))
            if len(theta.shape)==0:
                theta=np.ones((hsm.shape[0],1))*theta

        if len(theta)>len(hsm):
            theta=theta[0:len(hsm)]
        elif len(theta)<len(hsm):
            theta0=np.ones((hsm.shape[0],1))
            theta0[0:len(theta),0]=theta
            theta=theta0


        return theta

    def get_hsm(self,pp):
        if type(pp)!=type(str()):
            pp= str(pp)
        if ':' in pp:
            hsm=eval('np.r_['+pp+']')
        else:
            pp=pp.replace(' ',',')
            hsm=np.array(eval(pp))


        return hsm

    def get_nv_vqs(self,hsm,pp):
        if type(pp)!=type(str()):
            pp= str(pp)

        if ':' in pp:
            nv_vqs=eval('np.r_['+pp+']')
        elif 'N' in pp:
            nv_vqs=np.ones((hsm.shape[0]))
            for N in range(0,len(hsm)):
                nv_vqs[N]=eval(pp) # of levels for each master grid (increasing with depth)
        else:
            pp=pp.replace(' ',',')
            nv_vqs=np.array(eval(pp))
            if len(nv_vqs.shape)==0:
                nv_vqs=np.ones((hsm.shape[0]))*nv_vqs
            

        if len(nv_vqs)>len(hsm):
            nv_vqs=nv_vqs[0:len(hsm)]
        elif len(nv_vqs)<len(hsm):
            nv_vqs0=np.ones((hsm.shape[0]))
            nv_vqs0[0:len(nv_vqs),0]=nv_vqs
            nv_vqs=nv_vqs0


        return nv_vqs.astype(int)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MyApp(app)
    window.show()
    sys.exit(app.exec_())
