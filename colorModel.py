# -*- coding: utf-8 *-*
from __future__ import division
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pylab as plt
from math import factorial

from spectsens import spectsens
import PlottingFun as pf
from sompy import SOM


class colorModel():

    def __init__(self,
                ConeRatio={'fracLvM': 0.50, 's': 0.05, },
                maxSens={'l': 559.0, 'm': 530.0, 's': 417.0, }):

        if ConeRatio['fracLvM'] > 1 or ConeRatio['fracLvM'] < 0:
            raise IOError('Fraction of LvM must be between 0 and 1!')

        self.sRatio = ConeRatio['s']
        self.lRatio = (1 - self.sRatio) * (ConeRatio['fracLvM'])
        self.mRatio = (1 - self.sRatio) * (1 - self.lRatio)

        self.maxSens = maxSens

    def genModel(self):

        self.genFirstStage()
        self.genSecondStage()
        self.genThirdStage()

    def genFirstStage(self, startLambda=390, endLambda=750, step=1,
                        Out='anti-log'):
        """Compute the first stage in the model
        """

        lambdas = np.arange(startLambda, endLambda + step, step)

        L_cones = spectsens(LambdaMax=self.maxSens['l'], Output=Out,
                            StartWavelength=startLambda,
                            OpticalDensity=0.35,
                            EndWavelength=endLambda, Res=step)[1]
        M_cones = spectsens(LambdaMax=self.maxSens['m'], Output=Out,
                            StartWavelength=startLambda,
                            OpticalDensity=0.35,
                            EndWavelength=endLambda, Res=step)[1]
        S_cones = spectsens(LambdaMax=self.maxSens['s'], Output=Out,
                            StartWavelength=startLambda,
                            OpticalDensity=0.35,
                            EndWavelength=endLambda, Res=step)[1]

        self.FirstStage = {
            'lambdas': lambdas,
            'wavelen': {'startWave': startLambda, 'endWave': endLambda,
                        'step': step, },
            'L_cones': L_cones,
            'M_cones': M_cones,
            'S_cones': S_cones,
            }

    def genSecondStage(self):
        """Compute the second stage in the model
        """
        L_cones = self.FirstStage['L_cones']
        M_cones = self.FirstStage['M_cones']
        cones = {
            's': self.FirstStage['S_cones'],
            'm': self.FirstStage['M_cones'],
            'l': self.FirstStage['L_cones'], }

        self.SecondStage = {'lmsV_L': {}, 'lmsV_M': {}, 'ratio': {}, }
        i = 0
        for s in range(0, 101, 10):
            for m in range(0, 101, 10):
                for l in range(0, 101, 10):
                    if (l + m + s) == 100:
                        ratio = {'s': s, 'm': m, 'l': l, }
                        lmsV_L = self.optimizeChannel(cones, ratio,
                                                        Vcone=L_cones)
                        lmsV_M = self.optimizeChannel(cones, ratio,
                                                        Vcone=M_cones)
                        self.SecondStage['lmsV_L'][i] = lmsV_L
                        self.SecondStage['lmsV_M'][i] = lmsV_M
                        self.SecondStage['ratio'][i] = ratio
                        i += 1

    def genThirdStage(self):
        """Compute the third stage in the model
        """
        self.ThirdStage = {'redGreen': {}, 'blueYellow': {}, }

        gauss = lambda mu, k: (np.exp(-mu) * mu ** k) / factorial(k)
        gaussS, gaussM, gaussL = [], [], []
        for i in range(0, 101, 10):
            gaussS.append(gauss(self.sRatio, i))
            gaussM.append(gauss(self.mRatio, i))
            gaussL.append(gauss(self.lRatio, i))
        gaussS = gaussS / sum(gaussS)
        gaussM = gaussM / sum(gaussM)
        gaussL = gaussL / sum(gaussL)

        lCenterProb = (self.mRatio + self.lRatio) / self.lRatio
        self.ThirdStage = {
            'redGreen': np.zeros(len(self.SecondStage['lmsV_L'][0])),
            'blueYellow': np.zeros(len(self.SecondStage['lmsV_L'][0])),
            }
        for i in self.SecondStage['lmsV_L']:
            lRat = self.SecondStage['ratio'][i]['l'] / 100.
            mRat = self.SecondStage['ratio'][i]['m'] / 100.
            sRat = self.SecondStage['ratio'][i]['s'] / 100.

            prob = (sRat * gaussS[sRat] *
                    mRat * gaussM[mRat] *
                    lRat * gaussL[lRat])

            BY = self.SecondStage['lmsV_L'][i]
            RG = self.SecondStage['lmsV_M'][i]

            self.ThirdStage['redGreen'] += RG * prob * lCenterProb
            self.ThirdStage['blueYellow'] += BY * prob * (1 - lCenterProb)

    def optimizeChannel(self, cones, ratio, Vcone):

        lensMacula = self.getStockmanFilter()

        fun = lambda w, Vcone: (w * (ratio['s'] * cones['s'] +
                                            ratio['m'] * cones['m'] +
                                            ratio['l'] + cones['l']) -
                                         Vcone) / lensMacula
        # error function to minimize
        err = lambda w, Vcone: (fun(w, Vcone)).sum()

        con = fsolve(err, 1, args=(Vcone))
        out = fun(con, Vcone)

        return out

    def getStockmanFilter(self):

        lens = np.genfromtxt('stockman/lens.csv', delimiter=',')[::10]
        macula = np.genfromtxt('stockman/macular.csv', delimiter=',')[::10]

        return 10 ** lens[:361, 1] + 10 ** macula[:361, 1]

    def rectify(self, plot=True):

        self.red = self.final['redGreen']
        self.green = -1. * self.final['redGreen']
        self.blue = self.final['blueYellow']
        self.yellow = -1. * self.final['blueYellow']

        self.red[self.red < 0] = 0
        self.green[self.green < 0] = 0
        self.blue[self.blue < 0] = 0
        self.yellow[self.yellow < 0] = 0

        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            pf.AxisFormat()
            pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])

            ax.plot(self.ThirdStage['lambdas'], self.red, 'r', linewidth=3)
            ax.plot(self.ThirdStage['lambdas'], self.green, 'g', linewidth=3)
            ax.plot(self.ThirdStage['lambdas'], self.blue, 'b', linewidth=3)
            ax.plot(self.ThirdStage['lambdas'], self.yellow, 'y', linewidth=3)
            ax.set_xlabel('wavelength (nm)')
            ax.set_xlim([self.FirstStage['wavelen']['startWave'],
                         self.FirstStage['wavelen']['endWave']])
            plt.tight_layout()
            plt.show()

    def trainNeurons(self):

        ganglion = []

        for i in range(0, len(self.ThirdStage['lambdas'])):
            ganglion.append([self.red[i], self.green[i],
                             self.blue[i], self.yellow[i]])

        width = 10
        height = 10
        color_som = SOM(width, height, 4, 0.05)
        color_som.train(500, ganglion)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        pf.TufteAxis(ax, ['', ], [0, 0])
        plt.imshow(color_som.nodes, interpolation='nearest')
        plt.show()

        return color_som.nodes

    def returnFirstStage(self):
        return self.FirstStage

    def returnSecondStage(self):
        return self.SecondStage

    def returnThirdStage(self):
        return self.ThirdStage

    def returnFinal(self):
        return self.final


def plotModel(FirstStage, SecondStage, ThirdStage):
    """Plot cone spectral sensitivies and first stage predictions.
    """

    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_subplot(211)
    pf.AxisFormat()
    pf.TufteAxis(ax1, ['left', ], Nticks=[5, 5])

    ax1.plot(FirstStage['lambdas'], FirstStage['L_cones'],
            'r', linewidth=3)
    ax1.plot(FirstStage['lambdas'], FirstStage['M_cones'],
            'g', linewidth=3)
    ax1.plot(FirstStage['lambdas'], FirstStage['S_cones'],
            'b', linewidth=3)
    ax1.set_ylim([-0.05, 1.05])
    ax1.set_xlim([FirstStage['wavelen']['startWave'],
                  FirstStage['wavelen']['endWave']])
    #ax1.set_ylabel('sensitivity')

    if 'redGreen' in ThirdStage:

        ax3 = fig.add_subplot(212)
        pf.TufteAxis(ax3, ['left', 'bottom'], Nticks=[5, 5])
        ax3.plot(FirstStage['lambdas'], ThirdStage['redGreen'],
                'r', linewidth=3)
        ax3.plot(FirstStage['lambdas'], ThirdStage['blueYellow'],
                'b', linewidth=3)

    ax3.set_xlim([FirstStage['wavelen']['startWave'],
                     FirstStage['wavelen']['endWave']])
    #ax3.set_ylabel('activity')
    ax3.set_xlabel('wavelength (nm)')

    plt.tight_layout()
    plt.show()

    if 'lmsV_L' in SecondStage:
        fig = plt.figure(figsize=(8, 8))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        pf.AxisFormat()
        pf.TufteAxis(ax1, ['left', ], Nticks=[5, 5])
        pf.TufteAxis(ax2, ['left', 'bottom'], Nticks=[5, 5])

        for i in SecondStage['lmsV_L']:
            ax1.plot(FirstStage['lambdas'], SecondStage['lmsV_M'][i],
                    'b', linewidth=1)
            ax2.plot(FirstStage['lambdas'], SecondStage['lmsV_L'][i],
                    'r', linewidth=1)

        ax1.set_xlim([FirstStage['wavelen']['startWave'],
                      FirstStage['wavelen']['endWave']])
        ax2.set_xlim([FirstStage['wavelen']['startWave'],
                      FirstStage['wavelen']['endWave']])
        plt.tight_layout()
        plt.show()

if __name__ == '__main__':

    model = colorModel()
    model.genModel()
    FirstStage = model.returnFirstStage()
    SecondStage = model.returnSecondStage()
    ThirdStage = model.returnThirdStage()
    plotModel(FirstStage, SecondStage, ThirdStage)
    #model.rectify()
