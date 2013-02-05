# -*- coding: utf-8 *-*
from __future__ import division
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pylab as plt
#from math import factorial

from spectsens import spectsens
import PlottingFun as pf
from sompy import SOM


class colorModel():

    def __init__(self,
                ConeRatio={'fracLvM': 0.50, 's': 0.05, },
                maxSens={'l': 559.0, 'm': 530.0, 's': 417.0, }):

        self.test = True
        self.step = 1
        self.findConeRatios(ConeRatio['fracLvM'], ConeRatio['s'])
        self.maxSens = maxSens
        self.getStockmanFilter()

    def findConeRatios(self, fracLvM, fracS=None):
        if fracS > 1 or fracLvM < 0:
            raise IOError('Fraction of LvM must be between 0 and 1!')

        if fracS is not None:
            self.sRatio = fracS
        self.lRatio = (1 - self.sRatio) * (fracLvM)
        self.mRatio = (1 - self.sRatio) * (1 - fracLvM)

        if self.test:
            if round(self.sRatio + self.mRatio + self.lRatio, 7) != 1.0:
                print 'lms ratios: ', self.sRatio, self.mRatio, self.lRatio
                raise IOError('cone ratios must sum to 1.0!')

    def genModel(self):

        self.genFirstStage()
        self.genSecondStage()
        self.genThirdStage()

    def findUniqueHues(self):
        lambdas = self.FirstStage['lambdas']
        uniqueRed, uniqueGreen, uniqueBlue, uniqueYellow = [], [], [], []
        LMratio = []
        if not self.SecondStage:
            self.genSecondStage()
        else:

            for i in range(0, 101, self.step):

                self.findConeRatios(fracLvM=(i / 100.))
                self.genThirdStage()
                temp = self.returnThirdStage()
                RG = temp['lCenter']
                BY = temp['mCenter']

                if i == 0:
                    uniqueGreen.append(555)
                    uniqueRed.append(595)
                else:
                    zero_cross = np.where(np.diff(np.sign(RG)))[0]
                    uniqueGreen.append(lambdas[zero_cross[0]])
                    uniqueRed.append(lambdas[np.argmin(RG)])

                if i == 100:
                    uniqueBlue.append(420)
                    uniqueYellow.append(555)
                else:

                    zero_cross = np.where(np.diff(np.sign(BY)))[0]
                    uniqueBlue.append(lambdas[zero_cross[0]])
                    uniqueYellow.append(lambdas[zero_cross[1]])

                LMratio.append(i)

        self.uniqueHues = {
            'red': uniqueRed,
            'blue': uniqueBlue,
            'green': uniqueGreen,
            'yellow': uniqueYellow,
            'LMratio': LMratio,
            }

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
        for s in range(0, 101, self.step):
            for m in range(0, 101, self.step):
                for l in range(0, 101, self.step):
                    if (l + m + s) == 100:
                        ratio = {'s': s, 'm': m, 'l': l, }
                        lmsV_L = self.optimizeChannel(cones, ratio,
                                                        Center=L_cones)
                        lmsV_M = self.optimizeChannel(cones, ratio,
                                                        Center=M_cones)
                        self.SecondStage['lmsV_L'][i] = lmsV_L
                        self.SecondStage['lmsV_M'][i] = lmsV_M
                        self.SecondStage['ratio'][i] = ratio
                        i += 1

    def genThirdStage(self):
        """Compute the third stage in the model
        """

        #poisson = lambda mu, k: (np.exp(-mu) * mu ** k) / factorial(k)
        gauss = lambda mu, k: (np.exp(-1 * ((k - mu) ** 2) / 4 ** 2))
        gaussS, gaussM, gaussL = [], [], []
        for i in range(0, 101, self.step):
            gaussS.append(gauss(self.sRatio * 100., i))
            gaussM.append(gauss(self.mRatio * 100., i))
            gaussL.append(gauss(self.lRatio * 100., i))
        gaussS = gaussS / sum(gaussS)
        gaussM = gaussM / sum(gaussM)
        gaussL = gaussL / sum(gaussL)
        '''
        plt.figure()
        plt.plot(range(0, 101, self.step), gaussS, 'b', linewidth=3)
        plt.plot(range(0, 101, self.step), gaussM, 'g', linewidth=3)
        plt.plot(range(0, 101, self.step), gaussL, 'r', linewidth=3)
        plt.show()
        '''
        if self.test:
            if self.mRatio != 0 and self.lRatio != 0:
                if round(sum(gaussS) + sum(gaussM) + sum(gaussL), 5) != 3.0:
                    print 'sum s: ', sum(gaussS)
                    print 'sum m: ', sum(gaussM)
                    print 'sum l: ', sum(gaussL)
                    raise ValueError('gaussian prob distributions \
                        must sum to 1')

        lCenterProb = self.lRatio / (self.mRatio + self.lRatio)
        print lCenterProb
        self.ThirdStage = {
            'mCenter': np.zeros(len(self.SecondStage['lmsV_L'][0])),
            'lCenter': np.zeros(len(self.SecondStage['lmsV_M'][0])),
            'redGreen': np.zeros(len(self.SecondStage['lmsV_L'][0])),
            'blueYellow': np.zeros(len(self.SecondStage['lmsV_M'][0])),
            }
        p = 0
        for i in self.SecondStage['lmsV_L']:
            lRat = self.SecondStage['ratio'][i]['l'] / self.step
            mRat = self.SecondStage['ratio'][i]['m'] / self.step
            sRat = self.SecondStage['ratio'][i]['s'] / self.step

            prob = (gaussS[sRat] * gaussM[mRat] * gaussL[lRat])
            p += prob
            lCenter = self.SecondStage['lmsV_L'][i]
            mCenter = self.SecondStage['lmsV_M'][i]

            self.ThirdStage['mCenter'] += mCenter * prob * (1 - lCenterProb)
            self.ThirdStage['lCenter'] += lCenter * prob * (lCenterProb)

        self.ThirdStage['blueYellow'] = (self.ThirdStage['lCenter'] -
                                        self.ThirdStage['mCenter'])
        self.ThirdStage['redGreen'] = -1 * (self.ThirdStage['lCenter'] -
                                        self.ThirdStage['mCenter'])
        #print p

    def optimizeChannel(self, cones, ratio, Center):

        fun = lambda w, Center: (w * (ratio['s'] * cones['s'] +
                                    ratio['m'] * cones['m'] +
                                    ratio['l'] * cones['l']) -
                                 Center) / self.lensMacula
        # error function to minimize
        err = lambda w, Center: (fun(w, Center)).sum()
        w = fsolve(err, 1, args=(Center))
        out = fun(w, Center)

        if self.test:
            temp = err(w, Center)
            if temp > 1e-8:
                print ratio
                raise ValueError('error function not minimized properly')

        return out

    def getStockmanFilter(self):

        lens = np.genfromtxt('stockman/lens.csv', delimiter=',')[::10]
        macula = np.genfromtxt('stockman/macular.csv', delimiter=',')[::10]

        self.lensMacula = 10 ** lens[:361, 1] + 10 ** macula[:361, 1]

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
        color_som.train(1000, ganglion)

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

    def returnUniqueHues(self):
        return self.uniqueHues


def plotModel(FirstStage, SecondStage, ThirdStage, UniqueHues):
    """Plot cone spectral sensitivies and first stage predictions.
    """

    fig = plt.figure(figsize=(8, 9))
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

        ax2 = fig.add_subplot(212)
        pf.TufteAxis(ax2, ['left', ], Nticks=[5, 5])
        ax2.plot(FirstStage['lambdas'], ThirdStage['lCenter'],
                'r', linewidth=3)
        ax2.plot(FirstStage['lambdas'], ThirdStage['mCenter'],
                'b', linewidth=3)
        """
        ax2.plot(FirstStage['lambdas'], -1 * ThirdStage['lCenter'],
                'g', linewidth=3)
        ax2.plot(FirstStage['lambdas'], -1 * ThirdStage['mCenter'],
                'y', linewidth=3)

        ax3 = fig.add_subplot(313)
        pf.TufteAxis(ax3, ['left', 'bottom'], Nticks=[5, 5])
        ax3.plot(FirstStage['lambdas'], ThirdStage['redGreen'],
                'r', linewidth=3)
        ax3.plot(FirstStage['lambdas'], ThirdStage['blueYellow'],
                'b', linewidth=3)
        """
    ax2.set_xlim([FirstStage['wavelen']['startWave'],
                     FirstStage['wavelen']['endWave']])
    """
    ax3.set_xlim([FirstStage['wavelen']['startWave'],
                     FirstStage['wavelen']['endWave']])
    """
    #ax3.set_ylabel('activity')
    ax2.set_xlabel('wavelength (nm)')

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
            if i % 20 == 0:
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

    if 'green' in UniqueHues:
        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax1, ['left', 'bottom'], Nticks=[5, 5])

        ax1.plot(UniqueHues['LMratio'], UniqueHues['red'],
                'r', linewidth='3')
        ax1.plot(UniqueHues['LMratio'], UniqueHues['green'],
                'g', linewidth='3')
        ax1.plot(UniqueHues['LMratio'], UniqueHues['blue'],
                'b', linewidth='3')
        ax1.plot(UniqueHues['LMratio'], UniqueHues['yellow'],
                'y', linewidth='3')

        ax1.set_ylabel('wavelength (um)')
        ax1.set_xlabel('percent L vs M')

        plt.tight_layout()
        plt.show()

if __name__ == '__main__':

    model = colorModel()
    model.genModel()
    FirstStage = model.returnFirstStage()
    SecondStage = model.returnSecondStage()
    ThirdStage = model.returnThirdStage()

    model.findUniqueHues()
    UniqueHues = model.returnUniqueHues()
    plotModel(FirstStage, SecondStage, ThirdStage, UniqueHues)
    #model.rectify()
