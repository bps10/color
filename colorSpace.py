 # -*- coding: utf-8 *-*
from __future__ import division
import numpy as np

from spectsens import spectsens

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
        
    def genLMS(self, fundamental, LMSpeaks):
        '''
        '''
        if len(LMSpeaks) != 3:
            print 'LMSpeaks must be length 3! Using defaults: 559, 530, 417nm'
            LMSpeaks = [559.0, 530.0, 421.0]
            
        if fundamental.lower() == 'stockman':
            ind = len(self.spectrum)
            foo = np.genfromtxt(STATIC_ROOT + 
                                    '/stockman/fundamentals2deg.csv', 
                                 delimiter=',')[::10, :]
            self.Lc = 10.0 ** foo[:ind, 1]
            self.Mc = 10.0 ** foo[:ind, 2]
            self.Sc = 10.0 ** foo[:ind, 3]
    
            Lresponse = self.Lc * self.spectrum
            Mresponse = self.Mc * self.spectrum
            Sresponse = self.Sc * self.spectrum
            
        elif fundamental.lower() == 'stockspecsens':
            ind = len(self.spectrum)
            foo = np.genfromtxt(STATIC_ROOT + 
                                    '/stockman/specSens.csv', 
                                delimiter=',')[::10, :]

            LS = np.log10((1.0 - 10.0 ** -((10.0 ** foo[:, 1]) *
                    0.5)) / (1.0 - 10 ** -0.5))
            MS = np.log10((1.0 - 10.0 ** -((10.0 ** foo[:, 2]) *
                    0.5)) / (1.0 - 10 ** -0.5))
            SS = np.log10((1.0 - 10.0 ** -((10.0 ** foo[:, 3]) *
                    0.4)) / (1.0 - 10 ** -0.4))
          
            self.Lc = 10.0 ** LS[:ind]
            self.Mc = 10.0 ** MS[:ind]
            self.Sc = 10.0 ** SS[:ind]
            
            Lresponse = self.Lc / self.filters * self.spectrum
            Mresponse = self.Mc / self.filters * self.spectrum
            Sresponse = self.Sc / self.filters * self.spectrum
            
        elif fundamental.lower() == 'neitz':
            minspec = min(self.spectrum)
            maxspec = max(self.spectrum)
            self.Lc = spectsens(LMSpeaks[0], 0.5, 'anti-log', minspec, 
                                             maxspec, 1)[0]
            self.Mc = spectsens(LMSpeaks[1], 0.5, 'anti-log', minspec, 
                                             maxspec, 1)[0]
            self.Sc = spectsens(LMSpeaks[2], 0.4, 'anti-log', minspec, 
                                             maxspec, 1)[0]
                                                         
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
        lens = np.genfromtxt(STATIC_ROOT + '/stockman/lens.csv', 
                             delimiter=',')[::10, :]
        macula = np.genfromtxt(STATIC_ROOT + 
                                '/stockman/macular.csv', 
                                delimiter=',')[::10, :]

        spectrum = lens[:, 0]
        ind = np.where(spectrum == maxLambda)[0]
        self.spectrum = spectrum[:ind+1]
        
        self.filters = 10.0 ** (lens[:ind + 1, 1] +  macula[:ind + 1, 1])
        
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
        Xsens = spectsens(Xpeak, 0.5, 'anti-log', minspec, 
                                         maxspec, 1)[0]
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
        kaiser = np.genfromtxt('static/data/kaiser1987.csv', delimiter=",")
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
            
    def BY2lambda(self, propS, propM, verbose=False):
        '''
        '''
        l = 0
        m = -propM
        s = propS
        
        r, g, b = self.find_rgb(np.array([l, m, s]))
        line = self._lineEq(r, g)
        neutPoint = self._findDataIntercept(self.rVal, self.gVal, line)
        
        if verbose is True:
            return neutPoint, [r, g]
        else:
            return neutPoint
            
    def EE_CMFtoRGB(self):
        '''
        '''
        self.rVal, self.gVal, self.bVal = self.TrichromaticEquation(
                            self.CMFs[0, :], self.CMFs[1, :], self.CMFs[2, :])

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
        juddVos = np.genfromtxt('static/data/ciexyzjv.csv', delimiter=',')
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
        
        self._plotColorSpace(u, v, spectrum, ee=False)
        self.cs_ax.axes.get_xaxis().set_visible(False)
        self.cs_ax.axes.get_yaxis().set_visible(False)
        plt.axis('off')
        plt.show()

    def plotKaiser(self, neitz=False, confusion=None):
        '''
        
        '''
        try:
            from matplotlib.patches import Wedge
        except ImportError:
            raise ImportError('Sorry no patches module found')
            
        sub1_Neitz, sub2_Neitz, jv = self.genKaiser(neitz)

        if neitz:
            self._plotColorSpace()
            confusion = False
        else:
            self._plotColorSpace(rVal=jv[:, 0], gVal=jv[:, 1],
                                 spec=self.spectrum)
            self.cs_ax.plot([jv[-1, 0], jv[0, 0]], 
                            [jv[-1, 1], jv[0, 1]], 'k-', linewidth=3)
            self.cs_ax.set_ylim([0, 0.9])
            self.cs_ax.set_xlim([-0.05, 0.8])
            if not confusion:
                confusion = True
            
        self.cs_ax.plot(sub1_Neitz[:, 0], sub1_Neitz[:, 1], 'rx', 
                        markersize=8, linewidth=2)
        self.cs_ax.plot(sub2_Neitz[:, 0], sub2_Neitz[:, 1], 'bx',
                        markersize=8, linewidth=2)

        self.cs_ax.set_xlabel('x', fontsize=10)
        self.cs_ax.set_ylabel('y', fontsize=10)
        
        plt.tight_layout()
        plt.show()
        
        if confusion:
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
            print condition
            self.genLMS(fund=condition)
            
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
            
            neut, RG = self.BY2lambda(s, m, True)
            
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
            
    def plotConfusionLines(self, deficit='protan'):
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
        plt.show()                 

    def _plotColorSpace(self, rVal=None, gVal=None, spec=None, ee=True,
                        invert=True):
        '''
        '''      
        try:
            plt.__version__
        except NameError:
            import matplotlib.pylab as plt
        
        downSamp = 10
        
        if rVal == None or gVal == None or spec == None:
            rVal = self.rVal
            gVal = self.gVal
            spec = self.spectrum
            
            JuddV = False
            offset = 0.02
            turn = [500, 510]
            
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
        '''
        import matplotlib.patches as patch
        rgb = np.random.random((100,100))
        print rgb.shape
        im = plt.imshow(rgb, origin='lower', 
                interpolation='spline36',
                extent=([-2, 2, -2, 2]))
        path = patch.Path(list(zip(self.rVal, self.gVal)))
        p = patch.PathPatch(path, facecolor='none')
        im.set_clip_path(p)
        #self.cs_ax.add_patch(p)
        scaled_z = (z - z.min()) / z.ptp()
        colors = plt.cm.coolwarm(scaled_z)
        
        plt.scatter(x, y, marker='+', edgecolors=colors, s=150, linewidths=4)
        '''
        # annotate plot
        dat = zip(spec[::downSamp], rVal[::downSamp], gVal[::downSamp])

        for text, X, Y in dat:
            if text > 460 and text < 630:
                
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


if __name__ != '__main__':
    from NeitzModel import settings
    STATIC_ROOT = settings.STATIC_ROOT
## todo:: Create a logging function.

if __name__ == '__main__':
    
    STATIC_ROOT = './static'
    
    color = colorSpace()
    import matplotlib.pylab as plt
    import PlottingFun as pf
    #color.plotCompare()
    #color.plotFilters()
    #color.plotSpecSens()
    #color.plotCMFs()
    #color.plotcoeff()
    #color.plotColorSpace()
    #color.plotConfusionLines()
    #color.plotBYsystem(PRINT=True)
    color.plotKaiser(neitz=True)
    #color.trichromaticAnalysis()
    #color.tetrachromaticAnalysis()
    #color.plotConeSpace()
    color.plotLUV()