# -*- coding: utf-8 *-*
from __future__ import division
import numpy as np

from spectsens import spectsens
import PlottingFun as pf


class colorSpace():

    def __init__(self, stim='wright'):

        if stim.lower() == 'wright':
            self.lights = {'l': 650.0, 'm': 530.0, 's': 460.0, }
        if stim.lower() == 'stockman':
            self.lights = {'l': 645.0, 'm': 526.0, 's': 444.0, }

    def genLMS(self, LMSpeaks=[559.0, 530.0, 417.0]):

        if len(LMSpeaks) != 3:
            print 'LMSpeaks must be length 3! Using defaults: 559, 530, 417nm'
            LMSpeaks = [559.0, 530.0, 417.0]
        
        self.Lcone = lambda x: spectsens(LMSpeaks[0], 0.35, 'anti-log',
                                         x, x, 1)[0]
        self.Mcone = lambda x: spectsens(LMSpeaks[1], 0.35, 'anti-log',
                                         x, x, 1)[0]
        self.Scone = lambda x: spectsens(LMSpeaks[2], 0.35, 'anti-log',
                                         x, x, 1)[0]

    def genTrichromaticEquation(self):

        self.X = lambda R, G, B: R / (R + G + B)
        self.Y = lambda R, G, B: G / (R + G + B)
        self.Z = lambda R, G, B: B / (R + G + B)

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
        
    def genCMFs(self):

        self.genLMS()
        self.genTrichromaticEquation()
        self.genStockmanFilter()

        Lresponse = np.zeros((len(self.spectrum), 1))
        Mresponse = np.zeros((len(self.spectrum), 1))
        Sresponse = np.zeros((len(self.spectrum), 1))

        self.Lc = np.zeros((len(self.spectrum), 1))
        self.Mc = np.zeros((len(self.spectrum), 1))
        self.Sc = np.zeros((len(self.spectrum), 1))

        for i, light in enumerate(self.spectrum):

            self.Lc[i] = self.Lcone(light)
            self.Mc[i] = self.Mcone(light)
            self.Sc[i] = self.Scone(light)

            Lresponse[i] = self.Lc[i] / self.filters[i] / (light ** -1.0)
            Mresponse[i] = self.Mc[i] / self.filters[i] / (light ** -1.0)
            Sresponse[i] = self.Sc[i] / self.filters[i] / (light ** -1.0)

            if light == self.lights['l']:
                self.lights['lInd'] = i
            if light == self.lights['m']:
                self.lights['mInd'] = i
            if light == self.lights['s']:
                self.lights['sInd'] = i

        self.Lnorm = Lresponse / np.max(Lresponse)
        self.Mnorm = Mresponse / np.max(Mresponse)
        self.Snorm = Sresponse / np.max(Sresponse)
        LMSsens = np.array([self.Lnorm.T[0], self.Mnorm.T[0], self.Snorm.T[0]])

        convMatrix = np.array([
            [self.Lnorm[self.lights['lInd']][0],
            self.Lnorm[self.lights['mInd']][0],
            self.Lnorm[self.lights['sInd']][0]],

            [self.Mnorm[self.lights['lInd']][0],
            self.Mnorm[self.lights['mInd']][0],
            self.Mnorm[self.lights['sInd']][0]],

            [self.Snorm[self.lights['lInd']][0],
            self.Snorm[self.lights['mInd']][0],
            self.Snorm[self.lights['sInd']][0]]])

        print convMatrix
        
        self.CMFs = np.dot(np.linalg.inv(convMatrix), LMSsens)
        
        self.CMFs[0, :] /= sum(self.CMFs[0, :]) / 100.
        self.CMFs[1, :] /= sum(self.CMFs[1, :]) / 100.
        self.CMFs[2, :] /= sum(self.CMFs[2, :]) / 100.

        self.xVal = self.X(self.CMFs[0, :], self.CMFs[1, :], self.CMFs[2, :])
        self.yVal = self.Y(self.CMFs[0, :], self.CMFs[1, :], self.CMFs[2, :])
        self.zVal = self.Z(self.CMFs[0, :], self.CMFs[1, :], self.CMFs[2, :])

    def findRGB(self, RGB=[1,1,1]):
        tot = RGB / sum(RGB)
        print sum(tot)
    
    def returnCMFs(self):
        return {'cmfs': self.CMFs, 'wavelengths': self.spectrum, }
    
    def return_xyz(self):
        return self.xVal, self.yVal, self.zVal
                
    def plotColorSpace(self, plotFilters=False, plotSpecSens=False, 
                       plotCMFs=False):
        import matplotlib.pylab as plt

        if plotFilters:
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

        if plotSpecSens:
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

        if plotCMFs:
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

        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat(FONTSIZE=10, TickSize=6)
        pf.centerAxes(ax)

        # filter for positive y vals
        # ind = self.yVal >= 0
        ax.plot(self.xVal, self.yVal, 'k', linewidth=3)

        ax.plot(1.0/3.0, 1.0/3.0, 'ko', markersize=5)
        ax.annotate(s='{}'.format('S'), xy=(1./3.,1./3.), xytext=(2,8),
                    ha='right', textcoords='offset points', fontsize=14)
                    
        # annotate plot
        dat = zip(self.spectrum[::10], self.xVal[::10], self.yVal[::10])
        for text, X, Y in dat:
            if text > 460 and text < 630:

                if text <= 500: 
                    ax.scatter(X-0.02, Y, marker='_', s=150, c='k')
                elif text > 500 and text <= 510:
                    ax.scatter(X, Y+0.02, marker='|', s=150, c='k')
                else:
                    ax.scatter(X+0.02, Y, marker='_', s=150, c='k')
                   
                if text <510:
                    ax.annotate(s='{}'.format(int(text)),
                                xy=(X, Y), 
                                xytext=(-15, -5), 
                                ha='right', 
                                textcoords='offset points', fontsize=16)
                elif text >= 510 and text < 520:
                    ax.annotate(s='{}'.format(int(text)),
                                xy=(X, Y), 
                                xytext=(5, 20), 
                                ha='right', 
                                textcoords='offset points', fontsize=16)                    
                else:
                    ax.annotate(s='{}'.format(int(text)),
                                xy=(X, Y), 
                                xytext=(45, -5), 
                                ha='right', 
                                textcoords='offset points', fontsize=16)
                    
        #ax.set_ylim([-0.01, 2.25])
        plt.tight_layout()
        plt.show()


if __name__ == '__main__':

    color = colorSpace()
    color.genCMFs()
    color.plotColorSpace()
