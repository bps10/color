# -*- coding: utf-8 *-*
from __future__ import division
import numpy as np

from spectsens import spectsens
import PlottingFun as pf

## todo:: Create a logging function.

class colorSpace(object):

    def __init__(self, stim='wright'):

        if stim.lower() == 'wright':
            self.lights = {
                            'l': 650.0,
                            'm': 530.0,
                            's': 460.0,
                           }
        if stim.lower() == 'stockman':
            self.lights = {'l': 645.0, 
                           'm': 526.0, 
                           's': 444.0, }
        
    def genLMS(self, LMSpeaks=[559.0, 530.0, 421.0], fund='neitz'):
        '''
        '''
        if len(LMSpeaks) != 3:
            print 'LMSpeaks must be length 3! Using defaults: 559, 530, 417nm'
            LMSpeaks = [559.0, 530.0, 421.0]
        
        if fund.lower() == 'stockman':
            ind = len(self.spectrum)
            foo = np.genfromtxt('stockman/fundamentals2deg.csv', 
                                 delimiter=',')[::10, :]
            self.Lc = 10.0 ** foo[:ind, 1]
            self.Mc = 10.0 ** foo[:ind, 2]
            self.Sc = 10.0 ** foo[:ind, 3]
    
            Lresponse = self.Lc * self.spectrum
            Mresponse = self.Mc * self.spectrum
            Sresponse = self.Sc * self.spectrum
            
        elif fund.lower() == 'stockspecsens':
            ind = len(self.spectrum)
            foo = np.genfromtxt('stockman/specSens.csv', 
                                delimiter=',')[::10, :]

            LS = np.log10((1.0 - 10.0 ** -((10.0 ** foo[:, 1]) *
                    0.5)) / (1.0 - 10 ** -0.5))
            MS = np.log10((1.0 - 10.0 ** -((10.0 ** foo[:, 2]) *
                    0.5)) / (1.0 - 10 ** -0.5))
            SS = np.log10((1.0 - 10.0 ** -((10.0 ** foo[:, 3]) *
                    0.5)) / (1.0 - 10 ** -0.4))
          
            self.Lc = 10.0 ** LS[:ind]
            self.Mc = 10.0 ** MS[:ind]
            self.Sc = 10.0 ** SS[:ind]
            
            Lresponse = self.Lc / self.filters * self.spectrum
            Mresponse = self.Mc / self.filters * self.spectrum
            Sresponse = self.Sc / self.filters * self.spectrum
            
        elif fund.lower() == 'neitz':
            
            self.Lc = spectsens(LMSpeaks[0], 0.5, 'anti-log',
                                             min(self.spectrum), 
                                             max(self.spectrum), 1)[1]
            self.Mc = spectsens(LMSpeaks[1], 0.5, 'anti-log',
                                             min(self.spectrum), 
                                             max(self.spectrum), 1)[1]
            self.Sc = spectsens(LMSpeaks[2], 0.4, 'anti-log',
                                             min(self.spectrum), 
                                             max(self.spectrum), 1)[1]
                                                         
            Lresponse = self.Lc / self.filters * self.spectrum
            Mresponse = self.Mc / self.filters * self.spectrum
            Sresponse = self.Sc / self.filters * self.spectrum
            
        for i, light in enumerate(self.spectrum):

            if light == self.lights['l']:
                self.lights['lInd'] = i
            if light == self.lights['m']:
                self.lights['mInd'] = i
            if light == self.lights['s']:
                self.lights['sInd'] = i
        
        self.Lnorm = Lresponse / np.max(Lresponse)
        self.Mnorm = Mresponse / np.max(Mresponse)
        self.Snorm = Sresponse / np.max(Sresponse)

    def TrichromaticEquation(self, r, g, b):
        '''
        '''
        rgb = r + g + b
        r_ = r / rgb
        g_ = g / rgb
        b_ = b / rgb
        
        return r_, g_, b_

    def genStockmanFilter(self, maxLambda=770):
        '''
        '''
        lens = np.genfromtxt('stockman/lens.csv', delimiter=',')[::10, :]
        macula = np.genfromtxt('stockman/macular.csv', delimiter=',')[::10, :]

        spectrum = lens[:, 0]
        ind = np.where(spectrum == maxLambda)[0]
        self.spectrum = spectrum[:ind+1]
        
        self.filters = 10.0 ** (lens[:ind+1, 1] +  macula[:ind+1, 1])
        
    def LMStoCMFs(self):
        '''
        '''
        self.genStockmanFilter()
        self.genLMS()
        self.genConvMatrix()
        
        LMSsens = np.array([self.Lnorm.T, self.Mnorm.T, self.Snorm.T])
        self.CMFs = np.dot(np.linalg.inv(self.convMatrix), LMSsens)
        
    def genConvMatrix(self):
        '''
        '''
        self.convMatrix = np.array([
            [self.Lnorm[self.lights['lInd']],
            self.Lnorm[self.lights['mInd']],
            self.Lnorm[self.lights['sInd']]],

            [self.Mnorm[self.lights['lInd']],
            self.Mnorm[self.lights['mInd']],
            self.Mnorm[self.lights['sInd']]],

            [self.Snorm[self.lights['lInd']],
            self.Snorm[self.lights['mInd']],
            self.Snorm[self.lights['sInd']]]])

        print self.convMatrix
        

    def CMFtoEE_CMF(self):   
        '''
        '''
        self.CMFs[0, :] *= 100. / sum(self.CMFs[0, :]) 
        self.CMFs[1, :] *= 100. / sum(self.CMFs[1, :])
        self.CMFs[2, :] *= 100. / sum(self.CMFs[2, :])

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
    
    def find_testLight(self, testLight=600):
        '''
        '''
        rOut = np.interp(testLight, self.spectrum, self.rVal)
        gOut = np.interp(testLight, self.spectrum, self.gVal)
        bOut = np.interp(testLight, self.spectrum, self.bVal)

        return [rOut, gOut, bOut]        
        
    def find_rgb(self, LMS=np.array([1/3, 1/3, 1/3])):
        '''
        '''
        out = np.dot(np.linalg.inv(self.convMatrix), LMS)
        out = self.TrichromaticEquation(out[0], out[1], out[2])
        return out
        
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
        
    def _plotColorSpace(self):
        '''
        '''
        
        try:
            plt.__version__
        except NameError:
            import matplotlib.pylab as plt
            
        fig = plt.figure()
        self.cs_ax = fig.add_subplot(111)
        pf.AxisFormat(FONTSIZE=10, TickSize=6)
        pf.centerAxes(self.cs_ax)

        self.cs_ax.plot(self.rVal, self.gVal, 'k', linewidth=3.5)

        # add equi-energy location to plot
        self.cs_ax.plot(1.0/3.0, 1.0/3.0, 'ko', markersize=5)
        self.cs_ax.annotate(s='{}'.format('E'), xy=(1./3.,1./3.), xytext=(2,8),
                    ha='right', textcoords='offset points', fontsize=14)

        # annotate plot
        dat = zip(self.spectrum[::10], self.rVal[::10], self.gVal[::10])
        for text, X, Y in dat:
            if text > 460 and text < 630:

                if text <= 500: 
                    self.cs_ax.scatter(X - 0.02, Y, marker='_', s=150, c='k')
                    self.cs_ax.annotate(s='{}'.format(int(text)),
                                xy=(X, Y), 
                                xytext=(-15, -5), 
                                ha='right', 
                                textcoords='offset points', fontsize=16)
                elif text > 500 and text <= 510:
                    self.cs_ax.scatter(X, Y + 0.02, marker='|', s=150, c='k')
                    self.cs_ax.annotate(s='{}'.format(int(text)),
                                xy=(X, Y), 
                                xytext=(5, 20), 
                                ha='right', 
                                textcoords='offset points', fontsize=16) 
                else:
                    self.cs_ax.scatter(X + 0.02, Y, marker='_', s=150, c='k')
                    self.cs_ax.annotate(s='{}'.format(int(text)),
                                xy=(X, Y), 
                                xytext=(45, -5), 
                                ha='right', 
                                textcoords='offset points', fontsize=16)
        
        #self.cs_ax.set_xlim([-0.4, 1.2])
        #self.cs_ax.set_ylim([-0.2, 1.2])
        
        plt.tight_layout()
        #plt.show()

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
            self.cs_ax.plot([self.find_testLight(lam)[0],
                     self.copunctuals[deficit][0]],
                     [self.find_testLight(lam)[1], 
                      self.copunctuals[deficit][1]],
                     'k-', linewidth=1)   
                     
        self.cs_ax.text(0.7, 1, deficit, fontsize=18)
        plt.show()                 


if __name__ == '__main__':

    import matplotlib.pylab as plt
    color = colorSpace()
    color.LMStoCMFs()
    color.CMFtoEE_CMF()
    color.EE_CMFtoRGB()
    
    #color.plotCompare()
    #color.plotSpecSens()
    #color.plotCMFs()
    #color.plotcoeff()
    #color.plotColorSpace()
    color.plotConfusionLines()
    