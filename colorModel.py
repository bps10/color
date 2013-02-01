# -*- coding: utf-8 *-*
from __future__ import division
import numpy as np
import matplotlib.pylab as plt

from spectsens import spectsens
import PlottingFun as pf
from sompy import SOM


class colorModel():

    def __init__(self, LMSpercent={'l': 61, 'm': 32, 's': 7, },
                maxSens={'l': 559.0, 'm': 530.0, 's': 417.0, },
                startLambda=390, endLambda=750, step=1, physio=False):

        if startLambda < 390 or endLambda > 750:
            raise IOError('lambdas must between 390 and 830')
        if LMSpercent['l'] + LMSpercent['m'] + LMSpercent['s'] != 100:
            raise IOError('LMSpercent must sum to 100!')

        self.FirstStage = {}
        self.SecondStage = {}
        self.ThirdStage = {}

        self.genFirstStage(maxSens, startLambda, endLambda, step)
        self.genSecondStage(LMSpercent, physio)
        self.genThirdStage()

    def genFirstStage(self, maxSens, startLambda, endLambda, step,
                Out='anti-log'):
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

    def genSecondStage(self, LMSpercent, physio=False):

        L_cones = self.FirstStage['L_cones']
        M_cones = self.FirstStage['M_cones']
        S_cones = self.FirstStage['S_cones']

        if physio:
            LMS = (LMSpercent['l'] * L_cones +
            LMSpercent['m'] * M_cones +
            LMSpercent['s'] * S_cones)

            SBG_ON = (100 * S_cones) - LMS
            SBG_OFF = LMS - (100 * S_cones)

            mid_ON_Y = (100 * L_cones) - LMS
            mid_OFF_B = LMS - (100 * L_cones)

            mid_ON_G = (100 * M_cones) - LMS
            mid_OFF_R = LMS - (100 * M_cones)

            self.SecondStage = {
                'lambdas': self.FirstStage['lambdas'],
                'SBG_ON': SBG_ON,
                'SBG_OFF': SBG_OFF,
                'mid_ON_Y': mid_ON_Y,
                'mid_ON_G': mid_ON_G,
                'mid_OFF_R': mid_OFF_R,
                'mid_OFF_B': mid_OFF_B,
                }
        else:

            lensMacula = self.getStockmanFilter()

            const = {
                'S': 0.05,
                'Q1': 0.570641738067942,
                'Q2': 0.451415754843724,
                'Q3': 0.55983450938795,
                'Q4': 0.44016549061205,
                }

            smVl = ((const['Q1'] * (const['S'] * S_cones + (1. - const['S'])
                     * M_cones) - (1. - const['Q1']) * L_cones) / lensMacula)
            smVl = (smVl) - np.mean(smVl)

            slVm = ((const['Q2'] * (const['S'] * S_cones + (1. - const['S'])
             * L_cones) - (1. - const['Q2']) * M_cones) / lensMacula)
            slVm = slVm - np.mean(slVm)

            mVl = ((const['Q3'] * M_cones) - ((1. - const['Q3']) *
                     L_cones)) / lensMacula
            mVl = mVl - np.mean(mVl)

            lVm = ((const['Q4'] * L_cones) - ((1. - const['Q4'])
                    * M_cones)) / lensMacula
            lVm = lVm - np.mean(lVm)

            self.SecondStage = {
                'lambdas': self.FirstStage['lambdas'],
                'const': const,
                'smVl': smVl,
                'slVm': slVm,
                'lVm': lVm,
                'mVl': mVl,
                }

    def genThirdStage(self):

        if 'SBG_ON' in self.SecondStage:
            dark = (self.SecondStage['mid_OFF_R']
                    + self.SecondStage['mid_OFF_B']
                    + self.SecondStage['SBG_OFF'])

            light = (self.SecondStage['mid_ON_Y']
                    + self.SecondStage['mid_ON_G']
                    + self.SecondStage['SBG_ON'])

            red = (self.SecondStage['mid_ON_Y']
                    + self.SecondStage['mid_OFF_R']
                    + self.SecondStage['SBG_ON'])

            yellow = (self.SecondStage['mid_ON_Y']
                    + self.SecondStage['mid_OFF_R']
                    + self.SecondStage['SBG_OFF'])

            green = (self.SecondStage['mid_OFF_B']
                    + self.SecondStage['mid_ON_G']
                    + self.SecondStage['SBG_OFF'])

            blue = (self.SecondStage['mid_OFF_B']
                    + self.SecondStage['mid_ON_G']
                    + self.SecondStage['SBG_ON'])

            self.ThirdStage = {
                'lambdas': self.FirstStage['lambdas'],
                'dark': dark,
                'light': light,
                'red': red,
                'green': green,
                'blue': blue,
                'yellow': yellow,
                }

        if 'smVl' in self.SecondStage:

            redGreen = self.SecondStage['smVl'] - self.SecondStage['mVl']
            blueYellow = self.SecondStage['slVm'] - self.SecondStage['lVm']

            self.ThirdStage = {
                'lambdas': self.FirstStage['lambdas'],
                'redGreen': redGreen,
                'blueYellow': blueYellow,
                }

    def getStockmanFilter(self):

        lens = np.genfromtxt('stockman/lens.csv', delimiter=',')[::10]
        macula = np.genfromtxt('stockman/macular.csv', delimiter=',')[::10]

        return 10 ** lens[:361, 1] + 10 ** macula[:361, 1]

    def trainNeurons(self):

        ganglion = []
        for i in range(0, 1000):
            ganglion.append([self.SecondStage['SBG_ON'][i],
                             #self.SecondStage['SBG_OFF'][i],
                             self.SecondStage['mid_ON_Y'][i],
                             #self.SecondStage['mid_ON_G'][i],
                             self.SecondStage['mid_OFF_R'][i],
                             #self.SecondStage['mid_OFF_B'][i]
                             ])
        width = 1
        height = 1
        color_som = SOM(width, height, 3, 0.05)
        color_som.train(500, ganglion)
        plt.imshow(color_som.nodes, interpolation='nearest')
        plt.show()
        return color_som.nodes

    def returnFirstStage(self):
        return self.FirstStage

    def returnSecondStage(self):
        return self.SecondStage

    def returnThirdStage(self):
        return self.ThirdStage


def plotModel(FirstStage, SecondStage, ThirdStage):
    """Plot cone spectral sensitivies and first stage predictions.
    """

    fig = plt.figure(figsize=(8, 12))
    ax1 = fig.add_subplot(311)
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

    if 'SBG_ON' in SecondStage:

        ax2 = fig.add_subplot(312)
        pf.TufteAxis(ax2, ['left', ], Nticks=[5, 5])
        ax2.plot(SecondStage['lambdas'], SecondStage['SBG_ON'],
                'b', linewidth=3)
        ax2.plot(SecondStage['lambdas'], SecondStage['mid_ON_Y'],
                'r', linewidth=3)
        ax2.plot(SecondStage['lambdas'], SecondStage['mid_ON_G'],
                'g', linewidth=3)
        ax2.set_xlim([FirstStage['wavelen']['startWave'],
                      FirstStage['wavelen']['endWave']])

    if 'smVl' in SecondStage:

        ax2 = fig.add_subplot(312)
        pf.TufteAxis(ax2, ['left', ], Nticks=[5, 5])
        ax2.plot(SecondStage['lambdas'], SecondStage['smVl'],
                'b', linewidth=3)
        ax2.plot(SecondStage['lambdas'], SecondStage['slVm'],
                'r', linewidth=3)
        ax2.plot(SecondStage['lambdas'], SecondStage['mVl'],
                'orange', linewidth=3)
        ax2.plot(SecondStage['lambdas'], SecondStage['lVm'],
                'g', linewidth=3)
        ax2.set_xlim([FirstStage['wavelen']['startWave'],
                      FirstStage['wavelen']['endWave']])

    if 'red' in ThirdStage:

        ax3 = fig.add_subplot(313)
        pf.TufteAxis(ax3, ['left', 'bottom'], Nticks=[5, 5])
        ax3.plot(ThirdStage['lambdas'], ThirdStage['red'],
                'r', linewidth=3)
        ax3.plot(ThirdStage['lambdas'], ThirdStage['blue'],
                'b', linewidth=3)

    if 'redGreen' in ThirdStage:

        ax3 = fig.add_subplot(313)
        pf.TufteAxis(ax3, ['left', 'bottom'], Nticks=[5, 5])
        ax3.plot(ThirdStage['lambdas'], ThirdStage['redGreen'],
                'b', linewidth=3)
        ax3.plot(ThirdStage['lambdas'], ThirdStage['blueYellow'],
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
    plotModel(FirstStage, SecondStage, ThirdStage)
