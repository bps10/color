# -*- coding: utf-8 *-*
from __future__ import division
import numpy as np

from spectsens import spectsens
import PlottingFun as pf


## todo:: Create a logging function.

class colorSpace(object):

    def __init__(self, stim='wright'):

        if stim.lower() == 'wright':
            self.lights = {'l': 650.0, 'm': 530.0, 's': 460.0, }
        if stim.lower() == 'stockman':
            self.lights = {'l': 645.0, 'm': 526.0, 's': 444.0, }

    def genLMS(self, LMSpeaks=[561.0, 534.0, 422.0], fund='compare'):

        if len(LMSpeaks) != 3:
            print 'LMSpeaks must be length 3! Using defaults: 559, 530, 417nm'
            LMSpeaks = [559.0, 530.0, 421.0]
        
        if fund.lower() == 'stockman':
            #foo = np.genfromtxt('stockman/fundamentals2deg.csv', delimiter=',')[::10, :]
            foo = np.genfromtxt('stockman/specSens.csv', delimiter=',')[::10, :]
            self.Lc = 10.0 ** foo[:, 1]
            self.Mc = 10.0 ** foo[:, 2]
            self.Sc = 10.0 ** foo[:, 3]
    
            Lresponse = self.Lc / self.filters * self.spectrum
            Mresponse = self.Mc / self.filters * self.spectrum
            Sresponse = self.Sc / self.filters * self.spectrum
            
        if fund.lower() == 'neitz' or fund.lower() == 'compare':
            
            self.Lc = spectsens(LMSpeaks[0], 0.6, 'anti-log',
                                             min(self.spectrum), 
                                             max(self.spectrum), 1)[0]
            self.Mc = spectsens(LMSpeaks[1], 0.6, 'anti-log',
                                             min(self.spectrum), 
                                             max(self.spectrum), 1)[0]
            self.Sc = spectsens(LMSpeaks[2], 0.4, 'anti-log',
                                             min(self.spectrum), 
                                             max(self.spectrum), 1)[0]
                                             
            Lresponse = self.Lc / self.filters * self.spectrum
            Mresponse = self.Mc / self.filters * self.spectrum
            Sresponse = self.Sc / self.filters * self.spectrum
            
        if fund.lower() == 'compare':
            comp = [False, True, True]
            if comp[0] == True:
                foo = np.genfromtxt('stockman/fundamentals2deg.csv', delimiter=',')[::10, :]
                
                LN = 10.0 ** foo[:, 1]
                MN = 10.0 ** foo[:, 2]
                SN = 10.0 ** foo[:, 3]  
                LN *= self.spectrum
                MN *= self.spectrum
                SN *= self.spectrum
                LN /= max(LN)
                MN /= max(MN)
                SN /= max(SN)
                
            if comp[1] == True:
                foo = np.genfromtxt('stockman/specSens.csv', delimiter=',')[::10, :]
                LS = 10.0 ** foo[:, 1]
                MS = 10.0 ** foo[:, 2]
                SS = 10.0 ** foo[:, 3]  
                LS /= self.filters * self.spectrum
                MS /= self.filters * self.spectrum
                SS /= self.filters * self.spectrum
                LS /= sum(LS)
                MS /= sum(MS)
                SS /= sum(SS)

                
            if comp[2] == True:
                LN = spectsens(LMSpeaks[0], 0.6, 'anti-log',
                                                 min(self.spectrum), 
                                                 max(self.spectrum), 1)[0]
                MN = spectsens(LMSpeaks[1], 0.6, 'anti-log',
                                                 min(self.spectrum), 
                                                 max(self.spectrum), 1)[0]
                SN = spectsens(LMSpeaks[2], 0.4, 'anti-log',
                                                 min(self.spectrum), 
                                                 max(self.spectrum), 1)[0]
                                                 
                LN /= self.filters * self.spectrum
                MN /= self.filters * self.spectrum
                SN /= self.filters * self.spectrum
                LN /= sum(LN)
                MN /= sum(MN)
                SN /= sum(SN)
    
            fig = plt.figure()
            ax = fig.add_subplot(111)
            pf.AxisFormat()
            pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
            ax.plot(self.spectrum, LS, 'r--', linewidth=2)
            ax.plot(self.spectrum, LN, 'r-', linewidth=2)
            ax.plot(self.spectrum, MS, 'g--', linewidth=2)
            ax.plot(self.spectrum, MN, 'g-', linewidth=2)
            ax.plot(self.spectrum, SS, 'b--', linewidth=2)
            ax.plot(self.spectrum, SN, 'b-', linewidth=2)
            #ax.set_ylim([-0.01, 1.01])
            ax.set_xlim([380, 781])
            ax.set_xlabel('wavelength (nm)')
            ax.set_ylabel('sensitivity')        
            plt.tight_layout()
            plt.show()
            
        for i, light in enumerate(self.spectrum):

            if light == self.lights['l']:
                self.lights['lInd'] = i
                print self.Lc[i]
            if light == self.lights['m']:
                self.lights['mInd'] = i
                print self.Mc[i]
            if light == self.lights['s']:
                self.lights['sInd'] = i
                print self.Sc[i]


        
        print np.sum(Lresponse), np.sum(Mresponse), np.sum(Sresponse)
        
        self.Lnorm = Lresponse / np.max(Lresponse)
        self.Mnorm = Mresponse / np.max(Mresponse)
        self.Snorm = Sresponse / np.max(Sresponse)

    def genTrichromaticEquation(self):

        self.f_r = lambda r, g, b: r / (r + g + b)
        self.f_g = lambda r, g, b: g / (r + g + b)
        self.f_b = lambda r, g, b: b / (r + g + b)

    def genStockmanFilter(self):

        lens = np.genfromtxt('stockman/lens.csv', delimiter=',')[::10, :]
        macula = np.genfromtxt('stockman/macular.csv', delimiter=',')[::10, :]
        self.filters = 10.0 ** lens[:, 1] + 10.0 ** macula[:, 1]
        self.spectrum = lens[:, 0]

    def quanta2energy(self, quanta, lambdas):
        h = 6.628e-34
        c = 2.998e+08
        energy = (quanta * h * c) / 1e-9 * lambdas
        return energy
        
    def LMStoCMFs(self):
        
        self.genStockmanFilter()
        self.genLMS()
        self.genTrichromaticEquation()
        
        LMSsens = np.array([self.Lnorm.T, self.Mnorm.T, self.Snorm.T])

        convMatrix = np.array([
            [self.Lnorm[self.lights['lInd']],
            self.Lnorm[self.lights['mInd']],
            self.Lnorm[self.lights['sInd']]],

            [self.Mnorm[self.lights['lInd']],
            self.Mnorm[self.lights['mInd']],
            self.Mnorm[self.lights['sInd']]],

            [self.Snorm[self.lights['lInd']],
            self.Snorm[self.lights['mInd']],
            self.Snorm[self.lights['sInd']]]])

        print convMatrix
        
        self.CMFs = np.dot(np.linalg.inv(convMatrix), LMSsens)

    def CMFtoEE_CMF(self):   
        
        self.CMFs[0, :] *= 100. / sum(self.CMFs[0, :]) 
        self.CMFs[1, :] *= 100. / sum(self.CMFs[1, :])
        self.CMFs[2, :] *= 100. / sum(self.CMFs[2, :])

    def EE_CMFtoRGB(self):
        
        self.rVal = self.f_r(self.CMFs[0, :], self.CMFs[1, :], self.CMFs[2, :])
        self.gVal = self.f_g(self.CMFs[0, :], self.CMFs[1, :], self.CMFs[2, :])
        self.bVal = self.f_b(self.CMFs[0, :], self.CMFs[1, :], self.CMFs[2, :])
    
    def find_rgb(self, testLight=600, LMS=[1, 1, 1]):
        
        rOut = np.interp(testLight, self.spectrum, self.rVal)
        gOut = np.interp(testLight, self.spectrum, self.gVal)
        bOut = np.interp(testLight, self.spectrum, self.bVal)

        total = sum(LMS)
        rOut *= LMS[0] * total
        gOut *= LMS[1] * total        
        bOut *= LMS[2] * total        
        
        return [rOut, gOut, bOut]        
        
    def returnCMFs(self):
        return {'cmfs': self.CMFs, 'wavelengths': self.spectrum, }
    
    def return_xyz(self):
        return {'x': self.xVal, 'y': self.yVal, 'z': self.zVal, }

    def plotFilters(self):
        
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

    def plotColorSpace(self):
        
        try:
            plt.__version__
        except NameError:
            import matplotlib.pylab as plt
            
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat(FONTSIZE=10, TickSize=6)
        pf.centerAxes(ax)

        # filter for positive y vals
        # ind = self.yVal >= 0
        ax.plot(self.rVal, self.gVal, 'k', linewidth=3)

        # add equi-energy location to plot
        ax.plot(1.0/3.0, 1.0/3.0, 'ko', markersize=5)
        ax.annotate(s='{}'.format('E'), xy=(1./3.,1./3.), xytext=(2,8),
                    ha='right', textcoords='offset points', fontsize=14)
                    
        # annotate plot
        dat = zip(self.spectrum[::10], self.rVal[::10], self.gVal[::10])
        for text, X, Y in dat:
            if text > 460 and text < 630:

                if text <= 500: 
                    ax.scatter(X-0.02, Y, marker='_', s=150, c='k')
                    ax.annotate(s='{}'.format(int(text)),
                                xy=(X, Y), 
                                xytext=(-15, -5), 
                                ha='right', 
                                textcoords='offset points', fontsize=16)
                elif text > 500 and text <= 510:
                    ax.scatter(X, Y+0.02, marker='|', s=150, c='k')
                    ax.annotate(s='{}'.format(int(text)),
                                xy=(X, Y), 
                                xytext=(5, 20), 
                                ha='right', 
                                textcoords='offset points', fontsize=16) 
                else:
                    ax.scatter(X+0.02, Y, marker='_', s=150, c='k')
                    ax.annotate(s='{}'.format(int(text)),
                                xy=(X, Y), 
                                xytext=(45, -5), 
                                ha='right', 
                                textcoords='offset points', fontsize=16)
                    
        #ax.set_ylim([-0.01, 2.25])
        plt.tight_layout()
        plt.show()


if __name__ == '__main__':

    import matplotlib.pylab as plt
    color = colorSpace()
    color.LMStoCMFs()
    #color.plotSpecSens()
    color.CMFtoEE_CMF()
    color.EE_CMFtoRGB()
    color.plotCMFs()
    color.plotColorSpace()
