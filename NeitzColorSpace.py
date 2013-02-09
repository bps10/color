# -*- coding: utf-8 *-*
from __future__ import division
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
from spectsens import spectsens
import PlottingFun as pf


class colorSpace():

    def __init__(self):
        self.lights = {'l': 650, 'm': 530, 's': 460, }

    def genLMS(self):

        self.Lcone = lambda x: spectsens(559, 0.35, 'anti-log', x, x, 1)[0]
        self.Mcone = lambda x: spectsens(530, 0.35, 'anti-log', x, x, 1)[0]
        self.Scone = lambda x: spectsens(417, 0.35, 'anti-log', x, x, 1)[0]

    def genTrichromaticEquation(self):

        self.Xtri = lambda R, G, B: R / (R + G + B)
        self.Ytri = lambda R, G, B: G / (R + G + B)
        self.Ztri = lambda R, G, B: B / (R + G + B)

    def genCMFs(self):

        self.genLMS()
        self.genTrichromaticEquation()
        filters, spectrum = self.getStockmanFilter()

        Lresponse = np.zeros((len(spectrum), 1))
        Lval = np.zeros((len(spectrum), 1))
        Mresponse = np.zeros((len(spectrum), 1))
        Mval = np.zeros((len(spectrum), 1))
        Sresponse = np.zeros((len(spectrum), 1))
        Sval = np.zeros((len(spectrum), 1))
        Lc, Mc, Sc, filts = [], [], [], []
        
        i = 0
        for light in spectrum:
            Lc.append(self.Lcone(light))
            Mc.append(self.Mcone(light))
            Sc.append(self.Scone(light))
            
            Lresponse[i] = Lc[i] / filters[i]
            Mresponse[i] = Mc[i] / filters[i]
            Sresponse[i] = Sc[i] / filters[i]
            
            if light == self.lights['l']:
                self.lights['lFilt'] = i
            if light == self.lights['m']:
                self.lights['mFilt'] = i
            if light == self.lights['s']:
                self.lights['sFilt'] = filters[i]
            i += 1

        Lnorm = Lresponse / np.max(Lresponse) #* (spectrum ** -1.0)
        Mnorm = Mresponse / np.max(Mresponse) #* (spectrum ** -1.0)
        Snorm = Sresponse / np.max(Sresponse) #* (spectrum ** -1.0)

        convMatrix = np.array([
            [Lnorm[self.lights['lFilt']][0] / Lnorm[self.lights['sFilt']][0],
            Lnorm[self.lights['mFilt']][0] / Lnorm[self.lights['sFilt']][0],
            Lnorm[self.lights['sFilt']][0] / Lnorm[self.lights['sFilt']][0]],

            [Mnorm[self.lights['lFilt']][0] / Mnorm[self.lights['sFilt']][0],
            Mnorm[self.lights['mFilt']][0] / Mnorm[self.lights['sFilt']][0],
            Mnorm[self.lights['sFilt']][0] / Mnorm[self.lights['sFilt']][0]],

            [Snorm[self.lights['lFilt']][0] / Snorm[self.lights['sFilt']][0],
            Snorm[self.lights['mFilt']][0] / Snorm[self.lights['sFilt']][0],
            Snorm[self.lights['sFilt']][0] / Snorm[self.lights['sFilt']][0]]
                               ])
  
        print convMatrix
        CMFs = np.dot(np.linalg.inv(convMatrix),
                       np.array([Lnorm.T[0] / Lnorm[self.lights['sFilt']][0],
                                Mnorm.T[0] / Mnorm[self.lights['sFilt']][0],
                                Snorm.T[0] / Snorm[self.lights['sFilt']][0]]))
        #for i in range(3):
        #        CMFs[i,:] /= spectrum ** -1.0

        #for i in range(len(CMFs[0,:])):
        Ltri = self.Xtri(CMFs[0,:], CMFs[1,:], CMFs[2,:])
        Mtri = self.Ytri(CMFs[0,:], CMFs[1,:], CMFs[2,:])
        Stri = self.Ztri(CMFs[0,:], CMFs[1,:], CMFs[2,:])

        plotFilters, plotSpecSens = False, False
        if plotFilters:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            pf.AxisFormat()
            pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
            ax.semilogy(spectrum, filters, 'k', linewidth=2)
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
            ax.plot(spectrum, Lnorm, 'r', linewidth=2)
            ax.plot(spectrum, Lc, 'r--', linewidth=2)
            ax.plot(spectrum, Mnorm, 'g', linewidth=2)
            ax.plot(spectrum, Mc, 'g--', linewidth=2)
            ax.plot(spectrum, Snorm, 'b', linewidth=2)
            ax.plot(spectrum, Sc, 'b--', linewidth=2)

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
        ax.plot(spectrum, CMFs[0, :], 'r', linewidth=2)
        ax.plot(spectrum, CMFs[1, :], 'g', linewidth=2)
        ax.plot(spectrum, CMFs[2, :], 'b', linewidth=2)
        #ax.set_ylim([-0.01, 2.25])
        ax.set_xlabel('wavelength (nm)')
        ax.set_ylabel('sensitivity')
        plt.tight_layout()
        plt.show()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.centerAxes(ax)
        ax.plot(Ltri, Mtri, 'k', linewidth=3)
        #ax.set_ylim([-0.01, 2.25])
        plt.tight_layout()
        plt.show()

    def getStockmanFilter(self):

        lens = np.genfromtxt('stockman/lens.csv', delimiter=',')
        macula = np.genfromtxt('stockman/macular.csv', delimiter=',')
        filters = 10 ** lens[::10, 1] + 10 ** macula[::10, 1]
        wavelengths = lens[::10,0]
        return filters, wavelengths

if __name__ == '__main__':

    color = colorSpace()
    color.genCMFs()
