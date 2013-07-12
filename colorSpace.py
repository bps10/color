#! /usr/bin/env python
# -*- coding: utf-8 *-*
from __future__ import division
import numpy as np
import matplotlib.pylab as plt

from base import plot as pf
from base import spectsens as ss
from base import optics as op

class colorSpace(object):
    '''
    '''
    def __init__(self, stim='wright', fundamental='neitz', 
                 LMSpeaks=[559.0, 530.0, 421.0]):
        
        self.params = {'lights': stim.lower, }
        self.setLights(stim)
        self.genStockmanFilter()
        self.genLMS(fundamental, LMSpeaks)
        self.genConvMatrix()
        
        self.LMStoCMFs()
        self.CMFtoEE_CMF()
        self.EE_CMFtoRGB()
        
    def genLMS(self, fundamental, LMSpeaks=[559.0, 530.0, 421.0]):
        '''
        '''
        if len(LMSpeaks) != 3:
            print 'LMSpeaks must be length 3! Using defaults: 559, 530, 417nm'
            LMSpeaks = [559.0, 530.0, 421.0]
            
        if fundamental.lower() == 'stockman':
            ind = np.max(self.spectrum)

            sens = ss.stockmanfund(minLambda=390, maxLambda=ind)

            LS = sens[:, 0]
            MS = sens[:, 1]
            SS = sens[:, 2]
    
            Lresponse = LS * self.spectrum
            Mresponse = MS * self.spectrum
            Sresponse = SS * self.spectrum
            
        elif fundamental.lower() == 'stockspecsens':
            ind = np.max(self.spectrum)

            sens = ss.stockman(minLambda=390, maxLambda=ind)

            LS = sens[:, 0]
            MS = sens[:, 1]
            SS = sens[:, 2]
            
            Lresponse = LS / self.filters * self.spectrum
            Mresponse = MS / self.filters * self.spectrum
            Sresponse = SS / self.filters * self.spectrum
            
        elif fundamental.lower() == 'neitz':
            minspec = min(self.spectrum)
            maxspec = max(self.spectrum)
            self.Lc = ss.neitz(LMSpeaks[0], 0.5, False, minspec, 
                                             maxspec, 1)
            self.Mc = ss.neitz(LMSpeaks[1], 0.5, False, minspec, 
                                             maxspec, 1)
            self.Sc = ss.neitz(LMSpeaks[2], 0.4, False, minspec, 
                                             maxspec, 1)
                                            
            Lresponse = self.Lc / self.filters * self.spectrum
            Mresponse = self.Mc / self.filters * self.spectrum
            Sresponse = self.Sc / self.filters * self.spectrum
        
        #record param
        self.params['fundamentals'] = fundamental
        self.params['LMSpeaks']= LMSpeaks
        
        self.Lnorm = Lresponse / np.max(Lresponse)
        self.Mnorm = Mresponse / np.max(Mresponse)
        self.Snorm = Sresponse / np.max(Sresponse)

    def genStockmanFilter(self, maxLambda=770):
        '''
        '''
        self.filters, self.spectrum = op.filters.stockman(minLambda=380, 
            maxLambda=maxLambda, RETURN_SPECTRUM=True, 
            resolution=1)

    def genConvMatrix(self, PRINT=False):
        '''
        '''
        self.convMatrix = np.array([
            [np.interp(self.lights['l'], self.spectrum, self.Lnorm),
            np.interp(self.lights['m'], self.spectrum, self.Lnorm),
            np.interp(self.lights['s'], self.spectrum, self.Lnorm)],

            [np.interp(self.lights['l'], self.spectrum, self.Mnorm),
            np.interp(self.lights['m'], self.spectrum, self.Mnorm),
            np.interp(self.lights['s'], self.spectrum, self.Mnorm)],

            [np.interp(self.lights['l'], self.spectrum, self.Snorm),
            np.interp(self.lights['m'], self.spectrum, self.Snorm),
            np.interp(self.lights['s'], self.spectrum, self.Snorm)]])

        if PRINT == True:
            print self.convMatrix

    def genTetraConvMatrix(self, Xpeak):
        '''
        '''
        minspec = min(self.spectrum)
        maxspec = max(self.spectrum)
        Xsens = ss.neitz(Xpeak, 0.5, False, minspec, 
                                         maxspec, 1)
        Xresponse = Xsens / self.filters * self.spectrum
        Xnorm = Xresponse / np.max(Xresponse)
        
        lights = {'l': 600, 'm': 510, 's': 420, 'x': 720}
        convMatrix = np.array([
            [np.interp(lights['l'], self.spectrum, self.Lnorm),
            np.interp(lights['m'], self.spectrum, self.Lnorm),
            np.interp(lights['s'], self.spectrum, self.Lnorm),
            np.interp(lights['x'], self.spectrum, self.Lnorm)],

            [np.interp(lights['l'], self.spectrum, self.Mnorm),
            np.interp(lights['m'], self.spectrum, self.Mnorm),
            np.interp(lights['s'], self.spectrum, self.Mnorm),
            np.interp(lights['x'], self.spectrum, self.Mnorm)],

            [np.interp(lights['l'], self.spectrum, self.Snorm),
            np.interp(lights['m'], self.spectrum, self.Snorm),
            np.interp(lights['s'], self.spectrum, self.Snorm),
            np.interp(lights['x'], self.spectrum, self.Snorm)],

            [np.interp(lights['l'], self.spectrum, Xnorm),
            np.interp(lights['m'], self.spectrum, Xnorm),
            np.interp(lights['s'], self.spectrum, Xnorm),
            np.interp(lights['x'], self.spectrum, Xnorm)]])
        return convMatrix

    def genXYZ(self, plot=True):
        '''
        '''
        rgb = np.array([self.rVal, self.gVal, self.bVal])
        JuddVos = self._genJuddVos()

        convXYZ = np.array([[2.768892, 1.751748, 1.130160],
                            [1.000000, 4.590700, 0.060100],
                            [0,        0.056508, 5.594292]])
         
        neitzXYZ = np.dot(convXYZ, rgb)             
        xyzM = np.linalg.lstsq(neitzXYZ.T, JuddVos)[0]
        xyz = np.dot(xyzM, neitzXYZ)

        self.X = xyz[0, :]
        self.Y = xyz[1, :]
        self.Z = xyz[2, :]
        
        if plot:
            self._plotColorSpace(self.X, self.Y, self.spectrum)
            plt.show()        

    def genKaiser(self, neitz=False):
        '''
        '''
        kaiser = np.genfromtxt('data/kaiser1987.csv', delimiter=",")
        # sort by subject:
        subj1 = np.where(kaiser[:, 3] == 1)
        subj2 = np.where(kaiser[:, 3] == 2)
        s1_xy = kaiser[subj1, :2][0]
        s2_xy = kaiser[subj2, :2][0]
        # solve for z:
        s1_z = 1.0 - (s1_xy[:, 0] + s1_xy[:, 1])
        s2_z = 1.0 - (s2_xy[:, 0] + s2_xy[:, 1])
        # combine:
        s1_xyz = np.zeros((len(s1_xy), 3))
        s1_xyz[:, :2] = s1_xy
        s1_xyz[:, 2] = s1_z
        s2_xyz = np.zeros((len(s2_xy), 3))
        s2_xyz[:, :2] = s2_xy
        s2_xyz[:, 2] = s2_z

        jv = self._genJuddVos()

        if neitz:
            JV_Neitz_transMatrix = self._genJuddVos2Neitz(jv)

            sub1_Neitz = np.dot(JV_Neitz_transMatrix, s1_xyz.T).T 
            sub2_Neitz = np.dot(JV_Neitz_transMatrix, s2_xyz.T).T
        else:
            sub1_Neitz = s1_xyz                     
            sub2_Neitz = s2_xyz 
                        
        return sub1_Neitz, sub2_Neitz, jv        

    def setLights(self, stim):
        '''
        '''
        if (stim.lower() != 'wright' and stim.lower() != 'stiles and burch' 
            and stim.lower() != 'cie 1931'):
            print 'Sorry, stim light not understood, using wright'
            stim = 'wright' 
            
        if stim.lower() == 'wright':
            self.lights = {
                            'l': 650.0,
                            'm': 530.0,
                            's': 460.0,
                           }
        if stim.lower() == 'stiles and burch':
            self.lights = {'l': 645.0, 
                           'm': 526.0, 
                           's': 444.0, }

        if stim.lower() == 'cie 1931':
            self.lights = {'l': 700.0, 
                           'm': 546.1, 
                           's': 435.8, }

    def TrichromaticEquation(self, r, g, b):
        '''
        '''
        rgb = r + g + b
        r_ = r / rgb
        g_ = g / rgb
        b_ = b / rgb
        
        return r_, g_, b_
        
    def LMStoCMFs(self):
        '''
        '''
        LMSsens = np.array([self.Lnorm, self.Mnorm, self.Snorm])
        self.CMFs = np.dot(np.linalg.inv(self.convMatrix), LMSsens)

        #save sums for later normalization:            
        Rnorm = sum(self.CMFs[0, :])
        Gnorm = sum(self.CMFs[1, :])
        Bnorm = sum(self.CMFs[2, :])
        self.EEfactors = {'r': Rnorm, 'g': Gnorm, 'b': Bnorm, }
                           
    def CMFtoEE_CMF(self):
        '''
        '''
        self.CMFs[0, :], self.CMFs[1, :], self.CMFs[2, :] = self._EEcmf(
                                        self.CMFs[0, :], 
                                        self.CMFs[1, :], 
                                        self.CMFs[2, :])

    def lambda2BY(self, lam, verbose=False):
        '''
        '''
        r, g, b = self.find_testLightMatch(lam)
        line = self._lineEq(r, g)
        self.find_copunctuals()
        imagLine = self._lineEq(self.copunctuals['deutan'][0], 
                                self.copunctuals['deutan'][1],
                                self.copunctuals['tritan'][0], 
                                self.copunctuals['tritan'][1])
        xpoints = np.arange(self.copunctuals['tritan'][0], 
                            self.copunctuals['deutan'][0], 0.01)
        ypoints = imagLine(xpoints)
        neutPoint = self._findDataIntercept(xpoints, ypoints, line)

        if verbose is True:
            return neutPoint, [r, g]
        else:
            return neutPoint

    def lambda2RG(self, lam, equal_energy=True, verbose=False):
        '''
        '''
        self.find_copunctuals()
        imagLine = self._lineEq(self.copunctuals['protan'][0], 
                                self.copunctuals['protan'][1],
                                self.copunctuals['deutan'][0], 
                                self.copunctuals['deutan'][1])
        xpoints = np.arange(self.copunctuals['protan'][0], 
                            self.copunctuals['deutan'][0], 0.01)
        ypoints = imagLine(xpoints)

        r, g, b = self.find_testLightMatch(lam)
        if equal_energy:
            line = self._lineEq(r, g)
        else:
            ind = int(len(xpoints) * (2/3))

            line = self._lineEq(r, g, xpoints[ind], ypoints[ind])
            
        #print xpoints, ypoints
        neutPoint = self._findDataIntercept(xpoints, ypoints, line)

        if verbose is True:
            return neutPoint, [r, g]
        else:
            return neutPoint

    def BY2lambda(self, propS, propM, propL=0, verbose=False):
        '''
        '''
        l = propL
        m = -propM
        s = propS
        
        r, g, b = self.find_rgb(np.array([l, m, s]))
        line = self._lineEq(r, g)
        neutPoint = self._findDataIntercept(self.rVal, self.gVal, line)
        
        if verbose is True:
            return neutPoint, [r, g]
        else:
            return neutPoint

    def RG2lambda(self, propS, propM, propL=0, verbose=False):
        '''
        '''
        l = propL
        m = -propM
        s = propS
        
        r, g, b = self.find_rgb(np.array([l, m, s]))
        line = self._lineEq(r, g)
        neutPoint = self._findDataIntercept(self.rVal, self.gVal, line)
        
        if verbose is True:
            return neutPoint, [r, g]
        else:
            return neutPoint

    def EE_CMFtoRGB(self, rgb=None):
        '''
        '''
        if rgb is None:
            self.rVal, self.gVal, self.bVal = self.TrichromaticEquation(
                            self.CMFs[0, :], self.CMFs[1, :], self.CMFs[2, :])
        else:
            return self.TrichromaticEquation(rgb[0], rgb[1], rgb[2])

    def find_copunctuals(self):
        '''
        '''
        protan = self.find_rgb(np.array([1, 0, 0]))
        deutan = self.find_rgb(np.array([0, 1, 0]))
        tritan = self.find_rgb(np.array([0, 0, 1]))
        
        self.copunctuals = {'protan': protan, 
                            'deutan': deutan, 
                            'tritan': tritan, }
        
    def find_testLightMatch(self, testLight=600, R=None, G=None, B=None):
        '''
        '''
        if R == None or G == None or B == None:
            Lnorm = self.Lnorm
            Mnorm = self.Mnorm
            Snorm = self.Snorm
        else:
            Lnorm = R
            Mnorm = G
            Snorm = B
            
        l_ = np.interp(testLight, self.spectrum, Lnorm)
        m_ = np.interp(testLight, self.spectrum, Mnorm)
        s_ = np.interp(testLight, self.spectrum, Snorm)

        if R == None or G == None or B == None:
            rOut, gOut, bOut = self.find_rgb(LMS=np.array([l_, m_, s_]))
        else:
            rOut, gOut, bOut = l_, m_, s_
            
        return [rOut, gOut, bOut]        
    
    def find_testlightFromRG(self, r, g):
        '''
        '''
        err = lambda r, g, lam: ((r - self.rVal[lam])**2 + (g - 
                                self.gVal[lam])**2)
        i = 0
        startE = err(r, g, i)
        error = True
        try:
            while error:
                e = err(r, g, i)
                
                if startE < e and e < 10e-2:
                    error = False
                else:
                    startE = e
                    i += 1
                    
            #linear interpolate between points
            t0 = err(r, g, i)
            t1 = err(r, g, i + 1)
            outLam = self.spectrum[i] + (0 - t0 / (t1 - t0))
            return outLam
        except IndexError:
            raise IndexError('Pure light not found. Are you sure the [r, g] \
                            coords lie on the spectral line?')
        
    def find_rgb(self, LMS):
        '''
        '''
        cmf = np.dot(np.linalg.inv(self.convMatrix), LMS)
        cmf[0], cmf[1], cmf[2] = self._EEcmf(cmf[0], cmf[1], cmf[2])
        out = self.TrichromaticEquation(cmf[0], cmf[1], cmf[2])
        return out

    def find_BYweights(self):
        '''Function not finished
        '''
        neut, RG = self.BY2lambda(s, m, 0, True)
        n, rg  = self.BY2lambda(0.48, 0.52, 0, True)
        print self.find_testlightFromRG(n[0], n[1])

    def trichromaticAnalysis(self, Lmax=560, Smax=417):
        '''
        '''
        M_lamMax = []
        volume = []
        for i in range(420, 561):
            M_lamMax.append(i)
            self.genLMS('Neitz', [Lmax, i, Smax])
            self.genConvMatrix()
            volume.append(np.linalg.det(self.convMatrix))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
        
        ax.plot(M_lamMax, volume, 'k', linewidth=3)
        
        ind = np.where(np.array(M_lamMax) == 530)[0]
        ax.plot(M_lamMax[ind], volume[ind], 'g+', markersize=15, 
                markeredgewidth=3)
        
        ax.set_ylabel('color space volume')        
        ax.set_xlabel('$M_{\lambda_{max}}$') 
        
        ax.text(0.5, 0.5, 
            ('$L_{\lambda_{max}}$: ' + str(Lmax) + 'nm\n$S_{\lambda_{max}}$: '
            + str(Smax) + 'nm'), 
            fontsize=20, 
            horizontalalignment='center',
            verticalalignment='top',
            transform=ax.transAxes)
        
        plt.tight_layout()
        plt.show()

    def tetrachromaticAnalysis(self, Lmax=560, Mmax=530, Smax=417):
        '''
        '''
        X_lamMax = []
        volume = []
        self.genLMS('Neitz', [Lmax, Mmax, Smax])
        
        for i in range(300, 850):
            X_lamMax.append(i)
            
            tetraSystem = self.genTetraConvMatrix(i)
            volume.append(abs(np.linalg.det(tetraSystem)))
            '''
            if i % 10 == 0:
                print i, 'nm, rank: ', np.linalg.matrix_rank(tetraSystem)
            '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
        
        ax.plot(X_lamMax, volume, 'k', linewidth=3)
        
        ind = np.where(np.array(X_lamMax) == 495)[0]
        ax.semilogy(X_lamMax[ind], volume[ind], 'g+', markersize=15, 
                markeredgewidth=3)
        
        ax.set_ylabel('color space volume')        
        ax.set_xlabel('$X_{\lambda_{max}}$') 
        
        ax.text(0.2, 0.95, 
            ('$L_{\lambda_{max}}$: ' + str(Lmax) + 'nm\n$M_{\lambda_{max}}$: '
            + str(Mmax) + 'nm\n$S_{\lambda_{max}}$: ' + str(Smax) + 'nm'), 
            fontsize=20, 
            horizontalalignment='center',
            verticalalignment='top',
            transform=ax.transAxes)
        
        plt.tight_layout()
        plt.show()

    def _EEcmf(self, r_, g_, b_):   
        '''
        '''
        
        r_ *= 100. / self.EEfactors['r'] 
        g_ *= 100. / self.EEfactors['g']
        b_ *= 100. / self.EEfactors['b']
        
        return [r_, g_, b_]
        
    def _lineEq(self, x1, y1, x2=None, y2=None):
        '''Return the equation of a line from a given point that will pass
        through equal energy. Returns a function that takes one variable, x, 
        and returns y.
        '''
        if x2 == None:
            x2 = 1. / 3.
        if y2 == None:
            y2 = 1. / 3.
            
        m_ = (y2 - y1) / (x2 - x1)
        b_ = (y1) - (m_ * (x1))
        return lambda x: (m_ * x) + b_

    def _findDataIntercept(self, x, y, func):
        '''
        '''
        diff = True
        s = np.sign(func(x[0]) - y[0])
        i = 0
        while diff is True:
            err = func(x[i]) - y[i]
            sig = np.sign(err)
            if sig != s:
                diff = False
            else:
                i +=1
        #linear interpolate between points
        t0 = func(x[i - 1]) - y[i - 1]
        t1 = func(x[i]) - y[i]
        outX = x[i - 1] + ((x[i] - x[i - 1]) * (0 - t0 / (t1 - t0)))
        outY = func(outX)
        return [outX, outY]        

    def _genJuddVos2Neitz(self, juddVos):
        '''
        '''
        neitz = np.array([self.rVal, self.gVal, self.bVal]).T
        
        JuddVos_Neitz_lights = np.array([
            [np.interp(self.lights['l'], self.spectrum, juddVos[:, 0]),
            np.interp(self.lights['m'], self.spectrum, juddVos[:, 0]),
            np.interp(self.lights['s'], self.spectrum, juddVos[:, 0])],

            [np.interp(self.lights['l'], self.spectrum, juddVos[:, 1]),
            np.interp(self.lights['m'], self.spectrum, juddVos[:, 1]),
            np.interp(self.lights['s'], self.spectrum, juddVos[:, 1])],

            [np.interp(self.lights['l'], self.spectrum, juddVos[:, 2]),
            np.interp(self.lights['m'], self.spectrum, juddVos[:, 2]),
            np.interp(self.lights['s'], self.spectrum, juddVos[:, 2])]])

        foo = np.dot(np.linalg.inv(JuddVos_Neitz_lights), juddVos.T).T
        
        tempMat = np.linalg.lstsq(foo, neitz)[0]
        JuddVos_Neitz_transMatrix = np.dot(np.linalg.inv(JuddVos_Neitz_lights),
                        (tempMat))
        return JuddVos_Neitz_transMatrix
        
    def _genJuddVos(self):
        '''
        '''
        try:
            from scipy import interpolate as interp
        except ImportError:
            raise ImportError('Sorry cannot import scipy')
            
        #lights = np.array([700, 546.1, 435.8])
        juddVos = np.genfromtxt('data/ciexyzjv.csv', delimiter=',')
        spec = juddVos[:, 0]
        juddVos = juddVos[:, 1:]

        juddVos[:, 0] *= 100. / sum(juddVos[:, 0]) 
        juddVos[:, 1] *= 100. / sum(juddVos[:, 1])
        juddVos[:, 2] *= 100. / sum(juddVos[:, 2])
        r, g, b = self.TrichromaticEquation(juddVos[:, 0], 
                                            juddVos[:, 1],
                                            juddVos[:, 2])
        juddVos[:, 0], juddVos[:, 1], juddVos[:, 2] = r, g, b

        L_spline = interp.splrep(spec, juddVos[:, 0], s=0)
        M_spline = interp.splrep(spec, juddVos[:, 1], s=0)
        S_spline = interp.splrep(spec, juddVos[:, 2], s=0)
        L_interp = interp.splev(self.spectrum, L_spline, der=0)
        M_interp = interp.splev(self.spectrum, M_spline, der=0)
        S_interp = interp.splev(self.spectrum, S_spline, der=0)

        JVinterp = np.array([L_interp, M_interp, S_interp]).T      
        
        return JVinterp
        
    def returnConvMat(self):
        '''
        '''
        return self.convMatrix
    
    def returnCMFs(self):
        '''
        '''
        return {'cmfs': self.CMFs, 'wavelengths': self.spectrum, }
    
    def return_rgb(self):
        '''
        '''
        return {'r': self.rVal, 'g': self.gVal, 'b': self.bVal, }

    def plotConeSpace(self):
        '''
        '''
        self._plotColorSpace(self.Lnorm, self.Mnorm, self.spectrum)
        plt.show()

    def plotLUV(self):
        '''
        '''
        # make sure stim is cie 1931 and fundamentals neitz
        self.__init__(fundamental='neitz',
                                 LMSpeaks=[559.0, 530.0, 421.0],
                                 stim='cie 1931')
        
        self.genXYZ()
        u = 4 * self.X / (-2 * self.X + 12 * self.Y + 3)
        v = 9 * self.Y / (-2 * self.X + 12 * self.Y + 3)
        
        ind1 = np.where(self.spectrum == 420)[0]
        ind2 = np.where(self.spectrum == 700)[0]
        spectrum = self.spectrum[ind1:ind2+1]
        u = u[ind1:ind2]
        v = v[ind1:ind2]
        
        self._plotColorSpace(u, v, spectrum, ee=False, Luv=True, 
                             skipLam=[530, 550])
        self.cs_ax.axes.get_xaxis().set_visible(False)
        self.cs_ax.axes.get_yaxis().set_visible(False)
        plt.axis('off')
        plt.show()

    def plotKaiser(self, neitz=False, showBY=True, clip=True,
                   showSub1=False, showSub2=True, stockman=False,
                   series=True):
        '''
        
        '''
            
        sub1_Neitz, sub2_Neitz, jv = self.genKaiser(neitz)

        if neitz:
            self._plotColorSpace()

        else:
            self._plotColorSpace(rVal=jv[:, 0], gVal=jv[:, 1],
                                 spec=self.spectrum)
            self.cs_ax.plot([jv[-1, 0], jv[0, 0]], 
                            [jv[-1, 1], jv[0, 1]], 'k-', linewidth=3)
            self.cs_ax.set_ylim([0, 0.9])
            self.cs_ax.set_xlim([-0.05, 0.8])
        if showSub1:     
            self.cs_ax.plot(sub1_Neitz[:, 0], sub1_Neitz[:, 1], 'ko', 
                            markersize=8, markeredgewidth=2,
                            markerfacecolor='w',
                            linewidth=2)
        if showSub2:
            self.cs_ax.plot(sub2_Neitz[:, 0], sub2_Neitz[:, 1], 'kx',
                            markersize=8, markeredgewidth=2, linewidth=2)

        if showBY:
            if stockman:
                neut2, RG2 = self.BY2lambda(0, 0, 1., True)
                c2 = (1, 0, 0)
                c3 = (0.5, 0.5, 0)
                neut3, RG3 = self.lambda2RG(522, False, True)
            else:
                neut2, RG2 = self.BY2lambda(1, 0, 0, True)
                c2 = (0, 0, 1)
                c3 = (0, 0.5, 0.5)
                neut3, RG3 = self.lambda2BY(522, True)
            neut1, RG1 = self.BY2lambda(0, 1., 0, True)

            c1 = (0, 1, 0)
            # plot green copunctual line
            self.cs_ax.plot([neut1[0], RG1[0]], [neut1[1], RG1[1]], 
                            '-o', c=c1, markersize=8, linewidth=2)  
            # plot red or blue copunctual depending on neitz or stockman
            self.cs_ax.plot([neut2[0], RG2[0]], [neut2[1], RG2[1]], 
                            '-o', c=c2, markersize=8, linewidth=2)  
            # plot 
            self.cs_ax.plot([neut3[0], RG3[0]], [neut3[1], RG3[1]], 
                            '-o', c=c3, markersize=8, linewidth=2)  

        if stockman and series:

            for lam in [500, 505, 510, 515]:
                neut3, RG3 = self.lambda2RG(lam, False, True)
                self.cs_ax.plot([neut3[0], RG3[0]], [neut3[1], RG3[1]], 
                    '-o', c=c3, markersize=8, linewidth=2)  


        if clip is True:                
            self.cs_ax.set_xlim([-0.4, 1.2])
            self.cs_ax.set_ylim([-0.2, 1.2])
        
        self.cs_ax.set_xlabel('x', fontsize=10)
        self.cs_ax.set_ylabel('y', fontsize=10)
        
        plt.tight_layout()
        plt.show()

    def plotCIE(self):
        '''
        '''        
        try:
            from matplotlib.patches import Wedge
        except ImportError:
            raise ImportError('Sorry no patches module found')
            
        sub1_Neitz, sub2_Neitz, jv = self.genKaiser()
        ## plot confusion lines
        clip_area = Wedge((jv[0, 0], jv[0, 1]), r=10, theta1=0, theta2=360)
        CIEcopunctuals = {'deutan': np.array([1.10, -0.1, 0.1]),
                          'protan': np.array([0.753, 0.247, 0]), 
                          'tritan': np.array([0.17, 0, 0.83]),
                          }
        for deficit in CIEcopunctuals:
            self._plotColorSpace(rVal=jv[:, 0], gVal=jv[:, 1],
                                     spec=self.spectrum)     
            
            print deficit, ': ', CIEcopunctuals[deficit]
            
            if deficit.lower() == 'deutan' or deficit.lower() == 'protan':
                lambdas = [420, 460, 470, 480, 490, 500, 515,]
            elif deficit.lower() == 'tritan':
                lambdas = [420, 460, 480, 500, 520, 535, 545, 555,
                           570, 585, 600, 625, 700]
            
            self.cs_ax.plot(CIEcopunctuals[deficit][0],
                            CIEcopunctuals[deficit][1], 'ko', markersize=8)
            for lam in lambdas:
                R, G, B = jv[:, 0], jv[:, 1], jv[:, 2]
                self.cs_ax.plot([self.find_testLightMatch(lam, 
                                    R, G, B)[0],
                                 CIEcopunctuals[deficit][0]],

                                [self.find_testLightMatch(lam, 
                                    R, G, B)[1],
                                 CIEcopunctuals[deficit][1]],
                                'k-', linewidth=1) 
                self.cs_ax.plot([jv[-1, 0], jv[0, 0]], 
                        [jv[-1, 1], jv[0, 1]], 'k-', linewidth=3)
                        
                self.cs_ax.set_clip_path(clip_area)
                
            self.cs_ax.set_ylim([-0.12, 0.9])
            self.cs_ax.set_xlim([-0.05, 1.15])          
            self.cs_ax.set_xlabel('x', fontsize=10)
            self.cs_ax.set_ylabel('y', fontsize=10)
            self.cs_ax.text(0.8, 1, deficit, fontsize=18,
                            horizontalalignment='right',
                            verticalalignment='top',
                            transform=self.cs_ax.transAxes)
            plt.tight_layout()
            plt.show()        
                
    def plotCompare(self, compare=['stockman', 'stockSpecSens', 'neitz']):
        '''
            '''
        self.genStockmanFilter()
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
        style = ['-', '--', '-.']
        for i, condition in enumerate(compare):
            
            self.genLMS(fundamental=condition)
            
            ax.plot(self.spectrum, self.Lnorm, 'r' + style[i], linewidth=2)
            ax.plot(self.spectrum, self.Mnorm, 'g' + style[i], linewidth=2)
            ax.plot(self.spectrum, self.Snorm, 'b' + style[i], linewidth=2)
        #ax.set_ylim([-0.01, 1.01])
        ax.set_xlim([380, 781])
        ax.set_xlabel('wavelength (nm)')
        ax.set_ylabel('sensitivity')
        plt.tight_layout()
        plt.show()

    def plotFilters(self):
        '''
            '''
        try:
            plt.__version__
        except NameError:
            import matplotlib.pylab as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
        ax.semilogy(self.spectrum, self.filters, 'k', linewidth=2)
        ax.set_ylabel('log density')
        ax.set_xlabel('wavelength (nm)')
        ax.set_xlim([380, 781])
        ax.set_ylim([-10, max(self.filters)])
        plt.tight_layout()
        plt.show()

    def plotSpecSens(self):
        '''
            '''
        try:
            plt.show()
        except NameError:
            import matplotlib.pylab as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
        ax.plot(self.spectrum, self.Lnorm, 'r', linewidth=2)
        ax.plot(self.spectrum, self.Lc, 'r--', linewidth=2)
        ax.plot(self.spectrum, self.Mnorm, 'g', linewidth=2)
        ax.plot(self.spectrum, self.Mc, 'g--', linewidth=2)
        ax.plot(self.spectrum, self.Snorm, 'b', linewidth=2)
        ax.plot(self.spectrum, self.Sc, 'b--', linewidth=2)
        
        ax.set_ylim([-0.01, 1.01])
        ax.set_xlim([380, 781])
        ax.set_xlabel('wavelength (nm)')
        ax.set_ylabel('sensitivity')
        plt.tight_layout()
        plt.show()
    
    def plotCMFs(self):
        '''
            '''
        try:
            plt.__version__
        except NameError:
            import matplotlib.pylab as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
        ax.plot(self.spectrum, self.CMFs[0, :], 'r', linewidth=2)
        ax.plot(self.spectrum, self.CMFs[1, :], 'g', linewidth=2)
        ax.plot(self.spectrum, self.CMFs[2, :], 'b', linewidth=2)
        ax.set_xlim([self.spectrum[0], self.spectrum[-1]])
        ax.set_xlabel('wavelength (nm)')
        ax.set_ylabel('sensitivity')
        plt.tight_layout()
        plt.show()

    def plotcoeff(self):
        '''
            '''
        try:
            plt.__version__
        except NameError:
            import matplotlib.pylab as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
        ax.plot(self.spectrum, self.rVal, 'r', linewidth=2)
        ax.plot(self.spectrum, self.gVal, 'g', linewidth=2)
        ax.plot(self.spectrum, self.bVal, 'b', linewidth=2)
        ax.set_xlim([self.spectrum[0], self.spectrum[-1]])
        ax.set_xlabel('wavelength (nm)')
        ax.set_ylabel('coefficients')
        plt.tight_layout()
        plt.show()

    def plotColorSpace(self):
        '''
        '''
        self._plotColorSpace()
        plt.show()
    
    def plotBYsystem(self, PRINT=False, clip=True):
        '''
        '''
        self._plotColorSpace()
        
        for s in range(0, 11):
            m = (10.0 - s) / 10.0
            s = s / 10.0

            neut, RG = self.BY2lambda(s, m, 0, True)
            if PRINT is True:
                #print RG
                #print neut
                print self.find_testlightFromRG(neut[0], neut[1])
            self.cs_ax.plot([neut[0], RG[0]], [neut[1], RG[1]], 
                            '-o', c=(0, m, s), markersize=8, linewidth=2)
        
        if clip is True:                
            self.cs_ax.set_xlim([-0.4, 1.2])
            self.cs_ax.set_ylim([-0.2, 1.2])
        
        plt.show()
            
    def plotRGsystem(self, PRINT=False, clip=True):
        '''
        '''
        self._plotColorSpace()
        
        for l in range(0, 11):
            m = (10.0 - l) / 10.0
            l = l / 10.0
            
            neut, RG = self.RG2lambda(0, m, l, True)
            
            if PRINT is True:
                #print RG
                #print neut
                print self.find_testlightFromRG(neut[0], neut[1])
            self.cs_ax.plot([neut[0], RG[0]], [neut[1], RG[1]], 
                            '-o', c=(l, m, 0), markersize=8, linewidth=2)
        
        if clip is True:                
            self.cs_ax.set_xlim([-0.4, 1.2])
            self.cs_ax.set_ylim([-0.2, 1.2])
        
        plt.show()

    def plotConfusionLines(self, deficit='tritan', clip=True):
        '''add confusion lines
            '''
        
        self._plotColorSpace()
        self.find_copunctuals()
        print deficit, ': ', self.copunctuals[deficit]
        
        if deficit.lower() == 'deutan' or deficit.lower() == 'protan':
            lambdas = [420, 460, 470, 480, 490, 500, 515,]
        elif deficit.lower() == 'tritan':
            lambdas = [420, 460, 480, 500, 520, 535, 545, 555,
                       570, 585, 600, 625, 700]
        
        self.cs_ax.plot(self.copunctuals[deficit][0],
                        self.copunctuals[deficit][1], 'ko', markersize=8)
        for lam in lambdas:
            self.cs_ax.plot([self.find_testLightMatch(lam)[0],
                             self.copunctuals[deficit][0]],
                            [self.find_testLightMatch(lam)[1],
                             self.copunctuals[deficit][1]],
                            'k-', linewidth=1)   
        
        self.cs_ax.text(0.7, 1, deficit, fontsize=18)
        if clip is True:                
            self.cs_ax.set_xlim([-0.4, 1.2])
            self.cs_ax.set_ylim([-0.2, 1.2])
        plt.show()                 

    def genStockmanAnalysis(self, scale=0.34):
        '''
        '''
        self.genLMS('stockman', [559, 530, 421])

        stage2 = {}
        stage2['L-M'] = self.Lnorm - self.Mnorm
        stage2['M-L'] = self.Mnorm - self.Lnorm
        LM = self.Lnorm + (0.5 * self.Mnorm)
        stage2['S-LM'] = self.Snorm - (0.69 * LM)
        stage2['LM-S'] = (0.69 * LM) - self.Snorm

        stage3 = {}
        stage3['red'] = (2.55 * stage2['L-M']) + stage2['S-LM']
        stage3['green'] = (2.55 * stage2['M-L']) + stage2['LM-S']
        stage3['blue'] = (scale * (2.55 * stage2['M-L']) + 
            (self.Snorm - (scale * 0.69 * LM)))
        stage3['yellow'] = (scale * (2.55 * stage2['L-M']) +
            ((scale * 0.69 * LM) - self.Snorm))

        print 'RG sys: ', np.sum(stage3['red'])
        print 'BY sys: ', np.sum(stage3['blue'])

        return stage2, stage3

    def plotStockmanAnalysis(self):
        '''
        '''

        stage2, stage3 = self.genStockmanAnalysis()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
        for key in stage2:
            ax.plot(self.spectrum, stage2[key])

        ax.set_xlim([self.spectrum[0], self.spectrum[-1]])
        ax.set_xlabel('wavelength (nm)')
        ax.set_ylabel('sensitivity')
        plt.tight_layout()
        plt.show()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
        for key in stage3:
            ax.plot(self.spectrum, stage3[key], c=key)

        ax.set_xlim([self.spectrum[0], self.spectrum[-1]])
        ax.set_xlabel('wavelength (nm)')
        ax.set_ylabel('sensitivity')
        plt.tight_layout()
        plt.show()

        # Unique green series plot
        fig1 = plt.figure()
        fig2 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax2 = fig2.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax1, ['left', 'bottom'], Nticks=[5, 5])
        pf.TufteAxis(ax2, ['left', 'bottom'], Nticks=[5, 5])

        for j in [0.04, 0.34, 0.94, 2, 4, 6]:
            stage2, stage3 = self.genStockmanAnalysis(j)
            ax1.plot(self.spectrum, stage3['blue'], c='b', alpha=0.7)
            ax2.plot(self.spectrum, self.Snorm -
                j / 10 * (self.Lnorm + (0.5 * self.Mnorm)), 
                c='b', alpha=0.7)

        ax1.plot(self.spectrum, np.zeros(len(self.spectrum)), 'k',
            linewidth=1)
        ax2.plot(self.spectrum, np.zeros(len(self.spectrum)), 'k',
            linewidth=1)
        ax1.set_xlim([self.spectrum[0], 650])
        ax1.set_ylim([-0.7, 1.4])
        ax1.set_xlabel('wavelength (nm)')
        ax1.set_ylabel('sensitivity')

        ax2.set_xlim([self.spectrum[0], 700])
        ax2.set_ylim([-0.9, 1.2])
        ax2.set_xlabel('wavelength (nm)')
        ax2.set_ylabel('sensitivity')

        fig1.tight_layout()
        fig2.tight_layout()
        plt.show()

    def _plotColorSpace(self, rVal=None, gVal=None, spec=None, ee=True,
                        invert=False, Luv=False, skipLam=None, color=False):
        '''
        '''           
        downSamp = 10
        minLam = 460
        maxLam = 630
        
        if rVal == None or gVal == None or spec == None:
            rVal = self.rVal
            gVal = self.gVal
            spec = self.spectrum
            
            JuddV = False
            offset = 0.02
            turn = [500, 510]
        
        elif Luv:
            JuddV = False
            offset = 0.015
            turn = [500, 510]
            minLam = 420
            maxLam = 630            
        else:
            JuddV = True
            offset = 0.01
            turn = [510, 520]

        
        fig = plt.figure()
        self.cs_ax = fig.add_subplot(111)
        pf.AxisFormat(FONTSIZE=10, TickSize=6)
        
        if not JuddV:
            pf.AxisFormat(FONTSIZE=10, TickSize=6)
            pf.centerAxes(self.cs_ax)
        if JuddV:
            pf.AxisFormat(FONTSIZE=10, TickSize=8)
            pf.centerAxes(self.cs_ax)
            
        self.cs_ax.plot(rVal, gVal, 'k', linewidth=3.5)
        
        # add equi-energy location to plot
        if ee:
            self.cs_ax.plot(1.0/3.0, 1.0/3.0, 'ko', markersize=5)
            self.cs_ax.annotate(s='{}'.format('E'), xy=(1./3.,1./3.),
                                xytext=(2,8),
                                ha='right', textcoords='offset points',
                                fontsize=14)
        
        
        #rgb = np.reshape([self.Lnorm,self.Mnorm,self.Snorm], 
              #           [len(self.Lnorm) / 2, len(self.Lnorm) / 2, 3])

        if color:
            import clip as clip
            rgb = clip.clip_rgb_color(
                                np.array([self.rVal, self.gVal, self.bVal]))
            clip.rgb_patch_plot(self.cs_ax, rgb[0], color_names=None)
        

        
        # annotate plot
        dat = zip(spec[::downSamp], rVal[::downSamp], gVal[::downSamp])

        for text, X, Y in dat:
            if text > minLam and text < maxLam and not np.any(
                                        text == np.asarray(skipLam)):
                
                if text <= turn[0]:
                    self.cs_ax.scatter(X - offset, Y, marker='_', s=150, c='k')
                    self.cs_ax.annotate(s='{}'.format(int(text)),
                                        xy=(X, Y),
                                        xytext=(-15, -5),
                                        ha='right',
                                        textcoords='offset points', 
                                        fontsize=16)
                elif text > turn[0] and text <= turn[1]:
                    self.cs_ax.scatter(X, Y + offset, marker='|', s=150, c='k')
                    self.cs_ax.annotate(s='{}'.format(int(text)),
                                        xy=(X, Y),
                                        xytext=(5, 20),
                                        ha='right',
                                        textcoords='offset points',
                                        fontsize=16)
                else:
                    self.cs_ax.scatter(X + offset, Y, marker='_', s=150, c='k')
                    self.cs_ax.annotate(s='{}'.format(int(text)),
                                        xy=(X, Y),
                                        xytext=(45, -5),
                                        ha='right',
                                        textcoords='offset points', 
                                        fontsize=16)

        if invert:
            pf.invert(self.cs_ax, fig)
            
        plt.tight_layout()
        #plt.show()


def main(args):
    '''
    '''
    color = colorSpace()

    if args.Compare:
        color.plotCompare()

    if args.Filters:
        color.plotFilters()

    if args.SpecSens:
        color.plotSpecSens()

    if args.CMFs:
        color.plotCMFs()
    
    if args.coeff:
        color.plotcoeff()
    
    if args.ColorSpace:
        color.plotColorSpace()
    
    if args.ConfusionLines:
        color.plotConfusionLines()

    if args.BYsystem:
        color.plotBYsystem(False)

    if args.RGsystem:
        color.plotRGsystem(False)

    if args.Kaiser:
        color.plotKaiser(neitz=True, stockman=True)

    if args.Stockman:
        color.plotStockmanAnalysis()
    
    if args.tri:
        color.trichromaticAnalysis()
    
    if args.tetra:
        color.tetrachromaticAnalysis()
    
    if args.ConeSpace:
        color.plotConeSpace()
    
    if args.LUV:
        color.plotLUV()


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description="Color Space: display Neitz or Stockman\
        derived color spaces")
    
    parser.add_argument("-x", "--Compare", action="store_true", 
                        help="compare stockman and neitz fundamentals")
    parser.add_argument("-a", "--Filters", action="store_true",
                        help="plot lens and macula filters")
    parser.add_argument("-r", "--SpecSens", action="store_true", 
                         help="display spectral sensitivities")
    parser.add_argument("-y", "--CMFs", action="store_true",
                        help="plot color matching functions")    
    parser.add_argument("-f", "--coeff", action="store_true",
                        help="plot x,y,z coefficients")
    parser.add_argument("-i", "--ColorSpace", action="store_true",
                        help="plot color space")
    parser.add_argument("-c", "--ConfusionLines", action="store_true",
                        help="plot color space with confusion lines")
    parser.add_argument("-q", "--BYsystem", action="store_true",
                        help="plot blue-yellow system on color space")   
    parser.add_argument("-p", "--RGsystem", action="store_true",
                        help="plot red-green system on color space") 

    parser.add_argument("-k", "--Kaiser", action="store_true",
                        help="plot Kaiser data in Neitz or CIE space")
    parser.add_argument("-m", "--Stockman", action="store_true",
                        help="plot Stockman model.")
    parser.add_argument("-t", "--tri", action="store_true",
                        help="trichromatic analysis plot")
    parser.add_argument("-e", "--tetra", action="store_true",
                        help="tetrachromatic analysis plot")   
    parser.add_argument("-o", "--ConeSpace", action="store_true",
                        help="displace cone space plot")
    parser.add_argument("-l", "--LUV", action="store_true",
                        help="display best fit LUV space") 
    
    args = parser.parse_args()
    main(args)
 