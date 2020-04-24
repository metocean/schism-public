import sys,os
from base_io import BaseIO
import numpy as np
import argparse
import re

class VerticalGrid(BaseIO):
    """A class that manages I/O of VerticalGrid files
    """

    def __init__(self, logger=None):
        super(VerticalGrid, self).__init__()

        if logger:
            self.logger = logger



    def load(self,filename):
        '''Read in a vgrid.in file.'''

        self.logger.info("\tReading: %s" %filename)
        try:
            self.vgrid_file = open(filename)
            self.file = filename
        except:
            self.logger.error("Vgrid file %s not found" %filename)
            raise Exception("Vgrid file %s not found" %filename)

        self._read_header()
        self._read_Zlevel()
        self._read_Slevel()
        self._read_LSCcoord()

        return self


    def close(self):
        self.vgrid_file.close()

    def write(self, vfile='vgrid.in', hc=10.,thetab=0.,thetaf=0.0001,ivcor=2,nsigma=2,zcoord=[10000]):
        '''Create vertical grid [vgrid.in]'''

        sigma = np.arange(0, 1+(1./(nsigma-1)), 1./(nsigma-1))*-1
        sigma = sigma[0:nsigma][::-1]
        sigma[-1] = 0.0

        zcoord = sorted([np.abs(x)*-1 for x in zcoord])
        nz = len(zcoord)
        hs = zcoord[nz-1]*-1

        if ivcor == 1:  # LSC2
            self.logger.info("This type of vgrid is not implemented!")
            self.logger.error("This type of vgrid is not implemented!")
            return
            
        elif ivcor == 2:  # S/SZ
            a = open(vfile, 'w')
            a.write('%.f\n' % ivcor)
            a.write('%.f %.f %.f\n' % (nsigma+nz-1, nz, hs))
            a.write('Z levels\n')
            for kz in range(0, nz):
                a.write('%.f %.f\n' % (kz+1, zcoord[kz]))
            a.write('S levels\n')
            a.write('%.1f %.1f %f\n' %
                    (hc, thetab, thetaf))
            for ks in range(0, nsigma):
                a.write('%.f %.4f\n' % (ks+nz, sigma[ks]))
            a.close()

    def get_LSC(self, NP, dep, eta, h0):
        ''' Docstring '''

        coord = self._LSCcoord[NP]
        nvrt = int(self._nvrt)
        kbpl = int(coord[0])

        ztmp = [np.negative(dep)]
        for k in range(2, len(coord) - 1):
            ztmp.append((eta + dep) * coord[k] + eta)
        ztmp.append(eta)

        zz = [np.negative(dep) for k in range(0, kbpl - 1)]

        return zz + ztmp
        # kbpl=kbp(inode)
        # do k=kbpl+1,nvrt-1
        # 	ztmp(k)=(eta2(inode)+dp(inode))*sigma_lcl(k,inode)+eta2(inode)
        # enddo !k
        # ztmp(kbpl)=-dp(inode) !to avoid underflow
        # ztmp(nvrt)=eta2(inode) !to avoid underflow

    def get_pureS(self, dep, eta, h0):
        ''' Docstring '''

        z = list()
        idx1, idx2, idx3 = None, None, None

        if (dep <= self.hc):
            idx1 = True

        if eta <= (np.negative(self.hc) - (dep - self.hc) * self.thetas / np.sinh(self.thetas)):
            idx2 = True
            etam = 0.98 * (np.negative(self.hc) - (dep - self.hc)
                           * self.thetas / np.sinh(self.thetas))

        if (dep + eta) <= h0:
            idx3 = True

        for i in range(0, len(self._Scoord)):

            #---------------------------------------- compute C (sigma) -----------------------------
            sigma = self._Scoord[i]
            C = (1 - self.thetab) * np.sinh(self.thetas * sigma) / np.sinh(self.thetas) + \
                self.thetab * (np.tanh(self.thetas * (sigma+0.5)) - \
                np.tanh(self.thetas / 2)) / (2 * np.tanh(self.thetas / 2))

            if idx1:
                z.append(sigma * (dep + eta) + eta)
            else:
                z.append(eta * (1 + sigma) + self.hc * sigma + (dep - self.hc) * C)

            if idx2:
                sigma_hat = (z[i] - etam) / (dep + etam)
                z[i] = sigma_hat * (dep + eta) + eta
            elif idx3:
                z[i] = dep

        return z

    def _read_header(self):
        '''work to be done'''
        ## Header
        tokens, ok = self._read_and_parse_line(self.vgrid_file, 0)
        res = [float(s) for s in re.split('[\s,]+', tokens[0])]
        if len(res) > 1:
            self._nvrt = res[0]
            self._kz = res[1]
            self._h_s = res[2]
        else:
            self._ivcor = res[0]
            tokens, ok = self._read_and_parse_line(self.vgrid_file, 0)
            res = [float(s) for s in re.split('[\s,]+', tokens[0])]
            self._nvrt = res[0]
            
            if self._ivcor != 1:
                self._kz = res[1]
                self._h_s = res[2]
            else:
                self._kz = 0

        self._Zcoord, self._Scoord, self._LSCcoord = list(), list(), list()

    def _read_LSCcoord(self):
        '''work to be done'''
        self.close()
        vgrid_file = open(self.file)

        # ------------------------------------------------------------------
        
        tokens, ok = self._read_and_parse_line(vgrid_file, 0)
        res = [float(s) for s in re.split('[\s,]+', tokens[0])]
        if len(res) == 1 and self._ivcor == 1:
            tokens, ok = self._read_and_parse_line(vgrid_file, 0)
            ok = True
            while ok:
                tokens, ok = self._read_and_parse_line(vgrid_file, 0)
                if ok:
                    res = [float(s) for s in re.split('[\s,]+', tokens[0])]
                    self._LSCcoord.append(res[1:])

    def _read_Zlevel(self):
        '''work to be done'''
        ## Header
        # First line: ignored
        tokens, ok = self._read_and_parse_line(self.vgrid_file, 0)
        for k in range(0, int(self._kz)):
            tokens, ok = self._read_and_parse_line(self.vgrid_file, 0)
            res = [float(s) for s in re.split('[\s,]+', tokens[0])]
            self._Zcoord.append(res[1])

    def _read_Slevel(self):
        '''work to be done'''
        ## Header
        # First line: ignored
        tokens, ok = self._read_and_parse_line(self.vgrid_file, 0)
        if self._ivcor != 1:
            tokens, ok = self._read_and_parse_line(self.vgrid_file, 0)
            res = [float(kt) for kt in re.split('[\s,]+', tokens[0])]
            self.hc = res[0]
            self.thetab = res[1]
            self.thetas = res[2]
            ok = True
            while ok:
                tokens, ok = self._read_and_parse_line(self.vgrid_file, 0)
                if ok:
                    res = [float(s) for s in re.split('[\s,]+', tokens[0])]
                    self._Scoord.append(res[1])

    def sigma_to_zlayer(self, NP, dep, eta, h0):
        ''' Docstring '''

        if len(self._Zcoord) > 1:
            zs = self.get_pureS(np.min([dep, self._h_s]), eta, h0)
            zz = self._Zcoord[0:-1]
            z = np.array(zz + zs)
            z[z < np.negative(dep)] = np.negative(dep)

        elif len(self._Scoord) > 1:
            z = self.get_pureS(dep, eta, h0)

        elif len(self._LSCcoord) > 1:
            z = self.get_LSC(NP, dep, eta, h0)

        return z

    def link_to_root(self,rootidr):
        pass
    #     if not os.path.isfile(join(nest.rootdir,'vgrid.in')):
    # os.symlink(self.info.file, join(nest.rootdir,'vgrid.in') )

def run_as_script():
    parser = argparse.ArgumentParser(description="SCHISM grid reader")
    parser.add_argument('-g', '--hgrid', type=str,
                        help='grid file')
    parser.add_argument('-f', '--format', type=str,default='gr3',
                        help='grid format')
    parser.add_argument('-e', '--epsg', type=int, default=2193,
                        help="grid coordinate reference")

    args = parser.parse_args()
    
    
    gr=VerticalGrid(args.vgrid,format=args.format,epsg=args.epsg,logger=None)

if __name__ == '__main__':
    run_as_script()
