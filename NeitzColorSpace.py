# -*- coding: utf-8 *-*
from __future__ import division
import numpy as np

from spectsens import spectsens
import PlottingFun as pf


class colorSpace():

    def __init__(self):
        self.lights = {'l': 650.0, 'm': 530.0, 's': 460.0, }

    def genLMS(self):

        self.Lcone = lambda x: spectsens(559.0, 0.35, 'anti-log', x, x, 1)[0]
        self.Mcone = lambda x: spectsens(530.0, 0.35, 'anti-log', x, x, 1)[0]
        self.Scone = lambda x: spectsens(417.0, 0.35, 'anti-log', x, x, 1)[0]

    def genTrichromaticEquation(self):

        self.X = lambda R, G, B: R / (R + G + B)
        self.Y = lambda R, G, B: G / (R + G + B)
        self.Z = lambda R, G, B: B / (R + G + B)

    def getStockmanFilter(self):
        
        lens = np.genfromtxt('stockman/lens.csv', delimiter=',')[::10, :]
        macula = np.genfromtxt('stockman/macular.csv', delimiter=',')[::10, :]
        self.filters = 10.0 ** lens[:, 1] + 10.0 ** macula[:, 1]
        self.spectrum = lens[:,0]
    
    def genCMFs(self):

        self.genLMS()
        self.genTrichromaticEquation()
        self.getStockmanFilter()

        Lresponse = np.zeros((len(self.spectrum), 1))
        Lval = np.zeros((len(self.spectrum), 1))
        
        Mresponse = np.zeros((len(self.spectrum), 1))
        Mval = np.zeros((len(self.spectrum), 1))
        
        Sresponse = np.zeros((len(self.spectrum), 1))
        Sval = np.zeros((len(self.spectrum), 1))
        
        self.Lc = np.zeros((len(self.spectrum), 1))
        self.Mc = np.zeros((len(self.spectrum), 1))
        self.Sc = np.zeros((len(self.spectrum), 1))
        
        i = 0
        for light in self.spectrum:
            
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
        
            i += 1

        self.Lnorm = Lresponse / np.max(Lresponse)
        self.Mnorm = Mresponse / np.max(Mresponse)
        self.Snorm = Sresponse / np.max(Sresponse)

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
        self.CMFs = np.dot(np.linalg.inv(convMatrix),
                       np.array([self.Lnorm.T[0],
                                self.Mnorm.T[0],
                                self.Snorm.T[0]]))

        self.xVal = self.X(self.CMFs[0,:], self.CMFs[1,:], self.CMFs[2,:])
        self.yVal = self.Y(self.CMFs[0,:], self.CMFs[1,:], self.CMFs[2,:])
        self.zVal = self.Z(self.CMFs[0,:], self.CMFs[1,:], self.CMFs[2,:])

    def plotColorSpace(self, plotFilters=False, plotSpecSens=False):
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
            ax.set_ylim([-10, max(filters)])
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
        pf.AxisFormat()
        pf.centerAxes(ax)
        ax.plot(self.xVal, self.yVal, 'k', linewidth=3)

        #ax.set_ylim([-0.01, 2.25])
        plt.tight_layout()
        plt.show()




if __name__ == '__main__':

    color = colorSpace()
    color.genCMFs()
    color.plotColorSpace()
