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
                ConeRatio={'fracLvM': 0.25, 's': 0.05, },
                maxSens={'l': 559.0, 'm': 530.0, 's': 417.0, }):

        self.test = False
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
                    print zero_cross
                    uniqueBlue.append(lambdas[zero_cross[0]])
                    try:
                        uniqueYellow.append(lambdas[zero_cross[1]])
                    except:
                        uniqueYellow.append(599)
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

        self.SecondStage = {'lmsV_L': {}, 'lmsV_M': {}, 'percent': {}, }
        i = 0
        for s in range(0, 101, self.step):
            for m in range(0, 101, self.step):
                for l in range(0, 101, self.step):
                    if (l + m + s) == 100:
                        percent = {'s': s, 'm': m, 'l': l, }
                        lmsV_L = self.optimizeChannel(cones, percent,
                                                        Center=L_cones)
                        lmsV_M = self.optimizeChannel(cones, percent,
                                                        Center=M_cones)
                        self.SecondStage['lmsV_L'][i] = lmsV_L
                        self.SecondStage['lmsV_M'][i] = lmsV_M
                        self.SecondStage['percent'][i] = percent
                        i += 1

    def genThirdStage(self, plotGaussians=False):
        """Compute the third stage in the model
        """

        #poisson = lambda mu, k: (np.exp(-mu) * mu ** k) / factorial(k)
        gauss = lambda mu, x, SD: 1. / (SD * (2. * np.pi) ** 0.5) * np.exp(-
                                        (x - mu) ** 2. / (2. * SD ** 2.))
        binom = lambda k, n, p: ((factorial(n) /
                                (factorial(k) * factorial(n - k))
                                    * (p ** k)) * (1 - p) ** (n - k))

        gaussS = []
        for i in range(0, 101, self.step):
            gaussS.append(gauss(self.sRatio * 100., i, 2))  # + gaussS1[i])

        if plotGaussians:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            pf.AxisFormat()
            pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
            ax.plot(range(0, 101, self.step), gaussS, 'b', linewidth=3)
            plt.tight_layout()
            plt.show()

        lCenterProb = self.lRatio / (self.mRatio + self.lRatio)
        
        self.ThirdStage = {
            'mCenter': np.zeros(len(self.SecondStage['lmsV_L'][0])),
            'lCenter': np.zeros(len(self.SecondStage['lmsV_M'][0])),
            }
        p = 0
        for i in self.SecondStage['lmsV_L']:
            lNum = self.SecondStage['percent'][i]['l']
            mNum = self.SecondStage['percent'][i]['m']
            sNum = self.SecondStage['percent'][i]['s']

            probSur = (gauss(self.sRatio * 100, sNum / self.step, 2) *
                        binom(mNum, lNum + mNum, self.mRatio)                         
                            )
 
            p += probSur
            lCenter = self.SecondStage['lmsV_L'][i]
            mCenter = self.SecondStage['lmsV_M'][i]

            self.ThirdStage['mCenter'] += mCenter * (1 - lCenterProb) * probSur 
            self.ThirdStage['lCenter'] += lCenter * (lCenterProb) * probSur 
        '''
        self.ThirdStage['blueYellow'] = (-1 * self.ThirdStage['mCenter'] -
                                        -1 * self.ThirdStage['lCenter'])
        self.ThirdStage['redGreen'] = (self.ThirdStage['lCenter'] -
                                        self.ThirdStage['mCenter'])
        '''
        print self.sRatio, lCenterProb, 'prob :', p

        if self.test:
            if round(p, 2) != 1.0:
                print 'sum p: ', p
                raise ValueError('prob distribution must sum to 1')

    def optimizeChannel(self, cones, percent, Center):

        fun = lambda w, Center: (w * (percent['s'] / 100. * cones['s'] +
                                    percent['m'] / 100. * cones['m'] +
                                    percent['l'] / 100. * cones['l']) -
                                 Center) / self.lensMacula

        # error function to minimize
        err = lambda w, Center: (fun(w, Center)).sum()
        w = fsolve(err, 1, args=(Center))
        out = fun(w, Center)

        if self.test:
            temp = err(w, Center)
            if temp > 1e-8:
                print percent
                raise ValueError('error function not minimized properly')

        return out

    def getStockmanFilter(self):

        lens = np.genfromtxt('stockman/lens.csv', delimiter=',')[::10]
        macula = np.genfromtxt('stockman/macular.csv', delimiter=',')[::10]

        self.lensMacula = 10 ** (lens[:361, 1] + macula[:361, 1])

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

    if 'redGreen' in ThirdStage:
        fig = plt.figure(figsize=(8.5, 11))
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)

    else:
        fig = plt.figure(figsize=(8.5, 8))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)        

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
    ax1.set_ylabel('sensitivity')
    ax1.yaxis.set_label_coords(-0.2, 0.5)

    if 'lCenter' in ThirdStage:
        #pf.centerAxes(ax)        
        pf.TufteAxis(ax2, ['left', 'bottom'], Nticks=[5, 5])
        ax2.plot(FirstStage['lambdas'], 
                 np.zeros((len(FirstStage['lambdas']))), 'k', linewidth=1.0)
        ax2.plot(FirstStage['lambdas'], ThirdStage['lCenter'],
                'r', linewidth=3)
        ax2.plot(FirstStage['lambdas'], ThirdStage['mCenter'],
                'b', linewidth=3)
        ax2.set_xlim([FirstStage['wavelen']['startWave'],
                         FirstStage['wavelen']['endWave']])
        ax2.set_ylabel('activity')
        ax2.yaxis.set_label_coords(-0.2, 0.5)
        
    if 'redGreen' in ThirdStage:
        
        pf.TufteAxis(ax3, ['left', 'bottom'], Nticks=[5, 5])
        ax3.plot(FirstStage['lambdas'], -1 * ThirdStage['redGreen'],
                'r', linewidth=3)
        ax3.plot(FirstStage['lambdas'], -1 * ThirdStage['blueYellow'],
                'b', linewidth=3)
        ax3.set_xlim([FirstStage['wavelen']['startWave'],
                         FirstStage['wavelen']['endWave']])
        ax3.set_xlabel('wavelength (nm)')

    plt.tight_layout()
    plt.show()

    if 'lmsV_L' in SecondStage:
        fig = plt.figure(figsize=(8.5, 8))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        pf.AxisFormat()

        pf.TufteAxis(ax1, ['left', ], Nticks=[5, 5])
        pf.TufteAxis(ax2, ['left', 'bottom'], Nticks=[5, 5])

        ax1.plot(FirstStage['lambdas'], 
                 np.zeros((len(FirstStage['lambdas']))), 'k', linewidth=1.0)
        ax2.plot(FirstStage['lambdas'], 
                 np.zeros((len(FirstStage['lambdas']))), 'k', linewidth=1.0)
        for i in SecondStage['lmsV_L']:
            if i % 50 == 0:
                ax1.plot(FirstStage['lambdas'], SecondStage['lmsV_M'][i],
                        'b', linewidth=1, alpha=0.25)
                ax2.plot(FirstStage['lambdas'], SecondStage['lmsV_L'][i],
                        'r', linewidth=1, alpha=0.25)

        ax1.set_xlim([FirstStage['wavelen']['startWave'],
                      FirstStage['wavelen']['endWave']])
        ax2.set_xlim([FirstStage['wavelen']['startWave'],
                      FirstStage['wavelen']['endWave']])

        ax1.set_ylabel('sensitivity')
        ax2.set_ylabel('sensitivity')

        plt.tight_layout()
        plt.show()

    if 'green' in UniqueHues:
        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax1, ['left', 'bottom'], Nticks=[5, 5])

        ax1.plot(UniqueHues['LMratio'], UniqueHues['red'],
                'r', linewidth=3)
        ax1.plot(UniqueHues['LMratio'], UniqueHues['green'],
                'g', linewidth=3)
        ax1.plot(UniqueHues['LMratio'], UniqueHues['blue'],
                'b', linewidth=3)
        ax1.plot(UniqueHues['LMratio'], UniqueHues['yellow'],
                'y', linewidth=3)

        ax1.set_ylabel('wavelength ($\\mu$m)')
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
