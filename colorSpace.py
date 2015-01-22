#! /usr/bin/env python
# -*- coding: utf-8 *-*
from __future__ import division
import numpy as np
import matplotlib.pylab as plt

from base import plot as pf
from base import spectsens as ss
from base import optics as op

from genLMS import genLMS

class colorSpace(object):
    '''
    '''
    def __init__(self, stim='wright', fundamental='neitz', 
                 LMSpeaks=[559.0, 530.0, 421.0]):
        
        self.param = {'lights': stim.lower, 'stim': stim, 
                       'fund': fundamental, 'LMSpeaks': LMSpeaks}
        self.gen_space()

    def gen_space(self):
        '''
        '''
        # get fundamentals
        self.genfund(self.param['fund'])
        # gen fund to CMF matrix
        self.genConvMatrix()
        # convert fund to CMF
        self.LMStoCMFs()
        # convert CMF to Equal Energy
        self.CMFtoEE_CMF()
        # convert EE CMF to chromaticity coords (i.e. x + y + z = 1)
        self.EE_CMFtoRGB()

    
    def genfund(self, fundamental):
        '''
        ''' 
        self.fund = fundamental.lower()
        if self.fund in ['neitz', 'stockman']:
            self.genStockmanFilter()
            self.Lnorm, self.Mnorm, self.Snorm = genLMS(self.spectrum, 
                                            self.filters, 
                                            fundamental=self.fund, 
                                            LMSpeaks=self.param['LMSpeaks'])

        elif self.fund == 'smithpokorny' or self.fund == 'sp':
            sens, spectrum = ss.smithpokorny(minLambda=390, 
                                             maxLambda=720, 
                                             return_spect=True, 
                                             quanta=False)
            self.spectrum = spectrum
            self.Lnorm = sens[:, 0]
            self.Mnorm = sens[:, 1]
            self.Snorm = sens[:, 2]

        else:
            raise InputError('fundamentals not supported: must be neitz, \
            stockman or smithpokorny')

    def genStockmanFilter(self, maxLambda=770):
        '''
        '''
        self.filters, self.spectrum = op.filters.stockman(minLambda=390, 
            maxLambda=maxLambda, RETURN_SPECTRUM=True, 
            resolution=1)

    def genConvMatrix(self, PRINT=False):
        '''
        '''
        if self.fund in ['neitz', 'stockman']:
            self.setLights(stim)        
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
        else:
            self.convMatrix = ss.JuddVosCIE_to_sp()

        if PRINT == True:
            print self.convMatrix

    def genTetraConvMatrix(self, Xpeak):
        '''
        '''
        self.setLights(stim)        
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
            self.plotColorSpace(self.X, self.Y, self.spectrum)
            plt.show()           

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

    def plotColorSpace(self, rVal=None, gVal=None, spec=None, ee=True,
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
            if self.fund in ['neitz', 'stockman']:
                JuddV = False
                offset = 0.02
                turn = [500, 510]
            else:
                JuddV = True
                offset = 0.01
                turn = [510, 520]
        
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
        fig.set_tight_layout(True)
        self.cs_ax = fig.add_subplot(111)
        pf.AxisFormat(fontsize=10, ticksize=6)
        
        if JuddV:
            pf.AxisFormat(fontsize=10, ticksize=8)
            pf.TufteAxis(self.cs_ax, ['left', 'bottom'], [4, 4])
        else:
            pf.AxisFormat(fontsize=10, ticksize=6)
            pf.centerAxes(self.cs_ax)

        if color:
            import matplotlib.nxutils as nx
            verts = []
            for i, val in enumerate(rVal[:-10]):
                verts.append([rVal[i], gVal[i]])

            verts = np.asarray(verts)
            white = np.linalg.norm([1 / 3, 1 / 3, 1 / 3])
            for x in np.arange(-0.3, 1.1, 0.005):
                for y in np.arange(-0.15, 1.1, 0.01):
                    if x + y <= 1:
                        if nx.points_inside_poly(np.array([[x, y]]), verts):

                            _x = _boundval(x)
                            _y = _boundval(y)
                            _z = 1 - (_x + _y)

                            norm = np.linalg.norm([x, y, 1 - (x + y)])
                            dist = abs(norm - white)
                            
                            if dist <= (1 / 3):
                                _x += ((1 / 3) - dist)
                                _y += ((1 / 3) - dist)
                                _z += ((1 / 3) - dist)

                            self.cs_ax.plot((x), (y), 
                                'o', c=[_x, _y, _z], 
                                ms=6, mec='none', alpha=0.7)

            self.cs_ax.plot(rVal[:-10], gVal[:-10], 'k', linewidth=5)
            self.cs_ax.plot([rVal[0], rVal[-10]], [gVal[0], gVal[-10]], 'k', linewidth=5)

        self.cs_ax.plot(rVal[:-10], gVal[:-10], 'k', linewidth=3.5)
        self.cs_ax.plot([rVal[0], rVal[-10]], [gVal[0], gVal[-10]], 'k', linewidth=3.5)

        # add equi-energy location to plot
        if ee:
            self.cs_ax.plot(1.0/3.0, 1.0/3.0, 'ko', markersize=5)
            self.cs_ax.annotate(s='{}'.format('E'), xy=(1./3.,1./3.),
                                xytext=(2,8),
                                ha='right', textcoords='offset points',
                                fontsize=14)
        
        
        #rgb = np.reshape([self.Lnorm,self.Mnorm,self.Snorm], 
              #           [len(self.Lnorm) / 2, len(self.Lnorm) / 2, 3])

        
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
            
            
def _boundval(v):

    if v > 1:
        v = 1
    if v < 0:
        v = 0
    return round(v, 5)
 
