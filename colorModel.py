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
                ConeRatio={'fracLvM': 0.5, 's': 0.05, },
                maxSens={'l': 559.0, 'm': 530.0, 's': 417.0, }):

        if ConeRatio['fracLvM'] > 1 or ConeRatio['fracLvM'] < 0:
            raise IOError('Fraction of LvM must be between 0 and 1!')

        self.maxSens = maxSens
        self.ConeRatio = ConeRatio

        self.genModel()

    def genModel(self):

        self.genFirstStage(self.maxSens)

        collection = {'smVl': {}, 'slVm': {}, 'lVm': {}, 'mVl': {}}
        for i in range(0, 101, 1):
            temp = self.genSecondStage(ConeRatio={'fracLvM': i / 100.0,
                                                  's': 0.05, })

            collection['smVl'][i] = temp['smVl']
            collection['slVm'][i] = temp['slVm']
            collection['lVm'][i] = temp['lVm']
            collection['mVl'][i] = temp['mVl']
        self.SecondStage = collection

        self.genThirdStage()

        poisson = lambda mu, k: (np.exp(-mu) * mu ** k) / factorial(k)
        mu = self.ConeRatio['fracLvM'] * 100
        self.final = {'redGreen':
                        np.zeros(len(self.ThirdStage['redGreen'][i])),
                    'blueYellow':
                        np.zeros(len(self.ThirdStage['blueYellow'][i]))}

        for i in range(0, 101, 1):
            self.final['redGreen'] += (self.ThirdStage['redGreen'][i] *
                                        (poisson(mu, i) / 100.))
            self.final['blueYellow'] += (self.ThirdStage['blueYellow'][i] *
                                        (poisson(mu, i) / 100.))

    def genFirstStage(self, maxSens, startLambda=390,
                        endLambda=750, step=1, Out='anti-log'):
        """Compute the first stage in the model
        """

        lambdas = np.arange(startLambda, endLambda + step, step)

        L_cones = spectsens(LambdaMax=maxSens['l'], Output=Out,
                            StartWavelength=startLambda,
                            OpticalDensity=0.35,
                            EndWavelength=endLambda, Res=step)[1]
        M_cones = spectsens(LambdaMax=maxSens['m'], Output=Out,
                            StartWavelength=startLambda,
                            OpticalDensity=0.35,
                            EndWavelength=endLambda, Res=step)[1]
        S_cones = spectsens(LambdaMax=maxSens['s'], Output=Out,
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

    def genSecondStage(self, ConeRatio):
        """Compute the second stage in the model
        """
        L_cones = self.FirstStage['L_cones']
        M_cones = self.FirstStage['M_cones']
        S_cones = self.FirstStage['S_cones']

        smVl = self.optimizeChannel(cones={'s': S_cones, 'm': M_cones,
                             'l': L_cones, }, xCone=L_cones,
                                 ConeRatio=ConeRatio)

        slVm = self.optimizeChannel(cones={'s': S_cones, 'm': M_cones,
                             'l': L_cones, }, xCone=M_cones,
                                 ConeRatio=ConeRatio)

        mVl = self.optimizeChannel(cones={'s': 0, 'm': M_cones,
                             'l': L_cones, }, xCone=L_cones,
                                 ConeRatio=ConeRatio)

        lVm = self.optimizeChannel(cones={'s': 0, 'm': M_cones,
                             'l': L_cones, }, xCone=M_cones,
                                 ConeRatio=ConeRatio)

        SecondStage = {
            'lambdas': self.FirstStage['lambdas'],
            'smVl': smVl,
            'slVm': slVm,
            'lVm': lVm,
            'mVl': mVl,
            }
        return SecondStage

    def genThirdStage(self):
        """Compute the third stage in the model
        """
        self.ThirdStage = {'redGreen': {}, 'blueYellow': {},
                            'lambdas': self.FirstStage['lambdas']}

        for i in range(0, len(self.SecondStage['smVl'])):
            redGreen = self.SecondStage['smVl'][i] - self.SecondStage['mVl'][i]
            blueYellow = (self.SecondStage['slVm'][i] -
                        self.SecondStage['lVm'][i])

            self.ThirdStage['redGreen'][i] = redGreen
            self.ThirdStage['blueYellow'][i] = blueYellow

    def optimizeChannel(self, cones, xCone, ConeRatio):

        sRatio = ConeRatio['s']
        lRatio = (1 - sRatio) * ConeRatio['fracLvM']
        mRatio = (1 - sRatio) * (1 - ConeRatio['fracLvM'])

        lensMacula = self.getStockmanFilter()

        #if 'S' in cones:
        fun = lambda w, cones, xCone: (w * (sRatio * cones['s'] +
                                            mRatio * cones['m'] +
                                            lRatio + cones['l']) -
                                        (1 - w) * xCone) / lensMacula
        # error function to minimize
        err = lambda w, cones, xCone: (fun(w, cones, xCone)).sum()

        con = fsolve(err, 1, args=(cones, xCone))
        out = fun(con, cones, xCone)

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

    if 'smVl' in SecondStage:

        ax2 = fig.add_subplot(312)
        pf.TufteAxis(ax2, ['left', ], Nticks=[5, 5])
        ax2.plot(FirstStage['lambdas'], SecondStage['smVl'],
                'b', linewidth=3)
        ax2.plot(FirstStage['lambdas'], SecondStage['slVm'],
                'r', linewidth=3)
        ax2.plot(FirstStage['lambdas'], SecondStage['mVl'],
                'orange', linewidth=3)
        ax2.plot(FirstStage['lambdas'], SecondStage['lVm'],
                'g', linewidth=3)
        ax2.set_xlim([FirstStage['wavelen']['startWave'],
                      FirstStage['wavelen']['endWave']])

    if 'redGreen' in ThirdStage:

        ax3 = fig.add_subplot(212)
        pf.TufteAxis(ax3, ['left', 'bottom'], Nticks=[5, 5])
        ax3.plot(FirstStage['lambdas'], ThirdStage['redGreen'],
                'b', linewidth=3)
        ax3.plot(FirstStage['lambdas'], ThirdStage['blueYellow'],
                'r', linewidth=3)

    ax3.set_xlim([FirstStage['wavelen']['startWave'],
                     FirstStage['wavelen']['endWave']])
    #ax3.set_ylabel('activity')
    ax3.set_xlabel('wavelength (nm)')

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':

    model = colorModel()
    FirstStage = model.returnFirstStage()
    SecondStage = model.returnSecondStage()
    ThirdStage = model.returnThirdStage()
    final = model.returnFinal()
    plotModel(FirstStage, [], final)
    model.rectify()
