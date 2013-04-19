# -*- coding: utf-8 *-*
from __future__ import division
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pylab as plt
from math import factorial

from spectsens import spectsens
import PlottingFun as pf
#from sompy import SOM


class colorModel():
    '''
    '''
    def __init__(self, q=1.300):

        self.test = False
        self.step = 1
        self._q = q
        self.getStockmanFilter()

    def findConeRatios(self, fracLvM, fracS=None):
        '''
        '''
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

    def genModel(self, ConeRatio={'fracLvM': 0.75, 's': 0.05, },
                 maxSens={'l': 559.0, 'm': 530.0, 's': 417.0, }):
        '''
        '''
        self.findConeRatios(ConeRatio['fracLvM'], ConeRatio['s'])
        self.maxSens = maxSens
        
        self.genFirstStage()
        self.genSecondStage()
        self.genThirdStage()

    def findUniqueHues(self):
        '''
        '''
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
                RG = temp['mCenter']
                BY = temp['lCenter']

                if i == 0:
                    uniqueGreen.append(555)
                    uniqueRed.append(592)
                else:
                    zero_cross = np.where(np.diff(np.sign(BY)))[0]
                    uniqueGreen.append(lambdas[zero_cross[0]])
                    uniqueRed.append(lambdas[np.argmin(BY)])

                if i == 100:
                    uniqueBlue.append(474)
                    uniqueYellow.append(575)
                else:
                    zero_cross = np.where(np.diff(np.sign(RG)))[0]
                    uniqueBlue.append(lambdas[zero_cross[0]])
                    try:
                        uniqueYellow.append(lambdas[zero_cross[1]])
                    except:
                        uniqueYellow.append(600)
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
                            OpticalDensity=0.5,
                            EndWavelength=endLambda, Res=step)[0]
        L_cones /= self.lensMacula
        
        M_cones = spectsens(LambdaMax=self.maxSens['m'], Output=Out,
                            StartWavelength=startLambda,
                            OpticalDensity=0.5,
                            EndWavelength=endLambda, Res=step)[0]
        M_cones /= self.lensMacula
        
        S_cones = spectsens(LambdaMax=self.maxSens['s'], Output=Out,
                            StartWavelength=startLambda,
                            OpticalDensity=0.4,
                            EndWavelength=endLambda, Res=step)[0]
        S_cones /= self.lensMacula

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
                    if (l + m + s) == 100 and s == 5:
                        percent = {'s': s, 'm': m, 'l': l, }
                        lmsV_L = self.optimizeChannel(cones, percent,
                                                        Center=L_cones)
                        lmsV_M = self.optimizeChannel(cones, percent,
                                                        Center=M_cones)
                        self.SecondStage['lmsV_L'][i] = lmsV_L
                        self.SecondStage['lmsV_M'][i] = lmsV_M
                        self.SecondStage['percent'][i] = percent
                        i += 1

    def genThirdStage(self):
        """Compute the third stage in the model
        """
        #gauss = lambda mu, x, SD: 1. / (SD * (2. * np.pi) ** 0.5) * np.exp(-
        #                                (x - mu) ** 2. / (2. * SD ** 2.))
        binom = lambda k, n, p: ((factorial(n) /
                                (factorial(k) * factorial(n - k))
                                    * (p ** k)) * (1 - p) ** (n - k))
        #trinom = lambda l, m, s, L, M, S: (
        #                    factorial(L + M + S) /
        #                    (factorial(L) * factorial(M) * factorial(S)) *
        #                    (l ** L * m ** M * s ** S))
        lCenterProb = self.lRatio / (self.mRatio + self.lRatio)
        
        self.ThirdStage = {
            'mCenter': np.zeros(len(self.SecondStage['lmsV_L'][0])),
            'lCenter': np.zeros(len(self.SecondStage['lmsV_M'][0])),
            }
        p = 0
        for i in self.SecondStage['lmsV_L']:
            
            lNum = self.SecondStage['percent'][i]['l']
            mNum = self.SecondStage['percent'][i]['m']
            #sNum = self.SecondStage['percent'][i]['s']

            probSur = (#gauss(self.sRatio * 100, sNum / self.step, 0.5) *
                        #trinom(self.lRatio, self.mRatio, self.sRatio,
                        #   lNum, mNum, sNum)
                        binom(lNum, lNum + mNum, lCenterProb)                         
                            )
                            
            self.SecondStage['percent'][i]['probSurround'] = probSur
            
            p += probSur
            lCenter = self.SecondStage['lmsV_L'][i]
            mCenter = self.SecondStage['lmsV_M'][i]

            self.ThirdStage['mCenter'] += mCenter * (1 - lCenterProb) * probSur 
            self.ThirdStage['lCenter'] += lCenter * (lCenterProb) * probSur 

        #print self.sRatio, lCenterProb, 'prob :', p

        if self.test:
            if round(p, 2) != 1.0:
                print 'sum p: ', p
                raise ValueError('prob distribution must sum to 1')

    def optimizeChannel(self, cones, percent, Center):
        '''
        '''
        m_ = percent['m'] / (percent['m'] + percent['l'])
        l_ = percent['l'] / (percent['m'] + percent['l'])
        fun = lambda w, Center: (w * (self._q * cones['s'] +
                                    m_ * cones['m'] +
                                    l_ * cones['l']) -
                                 Center)

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
        
    def getStockmanFilter(self, maxLambda=750):
        '''
        '''
        lens = np.genfromtxt('static/stockman/lens.csv', delimiter=',')[::10]
        macula = np.genfromtxt('static/stockman/macular.csv', 
                               delimiter=',')[::10]
        spectrum = lens[:, 0]
        ind = np.where(spectrum == maxLambda)[0] + 1
                                       #just take upto a given index (750nm)
        self.lensMacula = 10 ** (lens[:ind, 1] + macula[:ind, 1])

    def returnFirstStage(self):
        '''
        '''
        return self.FirstStage

    def returnSecondStage(self):
        '''
        '''
        return self.SecondStage

    def returnThirdStage(self):
        '''
        '''
        return self.ThirdStage

    def returnUniqueHues(self):
        '''
        '''
        return self.uniqueHues


def optimizeUniqueHues():
    '''
    '''
    error= []
    parameter = np.arange(1.0, 1.5, 0.01)
    simY = np.array([0, 0, 0, 0.1, 0.4, 0,4, 0.1, 0, 0, 0, 0])
    for q in parameter:
        color = colorModel(q=q)
        freq, freqGreen, freqYellow, freqV = color.getCarroll_LMratios(
                            Volbrecht1997=True, returnVals=True, plot=False)
        
        error.append(np.sum(np.sqrt((freqGreen - freqV) ** 2)) 
                    +  np.sum(np.sqrt((freqYellow - simY) ** 2))
                    )
                    
        if q % 0.1 < 10e-6:
            print q
            
    print min(error), np.argmin(error), parameter[np.argmin(error)]
    np.savetxt('errors.txt', np.array([parameter, error]).T, delimiter='\t')


def LMratiosAnalysis(Volbrecht1997=True, returnVals=False, 
                        plot=True, savefigs=True):
    '''Creates a dictionary like object.
    '''
    carroll = np.genfromtxt('static/data/Carroll2002_lmRatios.txt', 
                        delimiter='\t', dtype=None, skip_header=0, 
                        names=True)
    model = colorModel(q=1.300)
    
    model.genModel()
    model.findUniqueHues()
    uniqueHues = model.returnUniqueHues()
    
    green, yellow, blue = [], [], []
    for subject in carroll['L']:
        green.append(np.interp(subject, uniqueHues['LMratio'], 
                          uniqueHues['green']))
        yellow.append(np.interp(subject, uniqueHues['LMratio'], 
                          uniqueHues['yellow']))
        blue.append(np.interp(subject, uniqueHues['LMratio'],
                              uniqueHues['blue']))

    print 'green: ', np.mean(green), np.std(green)
    print 'yellow: ', np.mean(yellow), np.std(yellow)
    print 'blue: ', np.mean(blue), np.std(blue)   
                       
    BINS = np.arange(0, 101, 5)
    if Volbrecht1997:
        BINS_G = np.arange(488, 564, 3)
        volb = np.genfromtxt('static/data/Volbrecht1997.txt', delimiter='\t',
                      dtype=None, skip_header=0, names=True)
    else:
        BINS_G = np.arange(490, 560, 5)
        
    BINS_Y = np.arange(514, 590, 3)
    BINS_B = np.arange(460, 536, 3)
    
    freq, bins = np.histogram(carroll['L'], bins=BINS)
    freqGreen, bins =np.histogram(green, bins=BINS_G)
    freqYellow, bins = np.histogram(yellow, bins=BINS_Y)
    freqBlue, bins = np.histogram(blue, bins=BINS_B)
                                        
    if plot:
        fig = plt.figure(figsize=(8.5, 11.5))
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        ax4 = ax3.twiny()
        
        pf.AxisFormat()
        plt.rc('legend',**{'fontsize': 12})
        pf.TufteAxis(ax1, ['left', 'bottom'], Nticks=[5, 5])
        pf.TufteAxis(ax2, ['left', 'bottom'], Nticks=[5, 5])
        pf.TufteAxis(ax3, ['left', 'bottom'], Nticks=[5, 5])
        pf.TufteAxis(ax4, ['left', 'bottom'], Nticks=[5, 5])

        ax3.spines['bottom'].set_position(('outward', 55))
        
        
        BINS, freq = pf.histOutline(freq / sum(freq), BINS)
        BINS_Gout, freqGreen = pf.histOutline(freqGreen / sum(freqGreen), 
                                           BINS_G)
        BINS_Yout, freqYellow = pf.histOutline(
                                freqYellow / sum(freqYellow), BINS_Y)
        BINS_Bout, freqBlue = pf.histOutline(freqBlue / sum(freqBlue), 
                                            BINS_B)

        ax1.plot(BINS, freq, 'k', linewidth=3)
        
        if Volbrecht1997:
            binV, freqV = pf.histOutline(volb['count'] / sum(
                                        volb['count']), BINS_G)
            ax2.plot(binV, freqV, c='0.8', linewidth=3,
                             label='Volbrecht 1997')
            ax2.fill_between(binV, freqV, 0, color='0.8')
            
        ax2.plot(BINS_Gout, freqGreen, 'g', linewidth=3,
                 label='predicted')
        ax3.plot(BINS_Yout, freqYellow, 'y', linewidth=3) 
        ax4.plot(BINS_Bout, freqBlue, 'b', linewidth=3) 
          
        ax1.set_xlim([0, 100])
        ax1.set_ylim([-0.002, max(freq) + 0.01])
        ax1.set_ylabel('proportion')
        ax1.set_xlabel('% L v M')
        ax1.yaxis.set_label_coords(-0.2, 0.5)
        
        ax2.set_ylabel('proportion')
        ax2.set_xlabel('unique green (nm)')
        ax2.set_xlim([490, 560])
        ax3.set_xlim([514, 590])
        ax4.set_xlim([460, 536])
        if Volbrecht1997:
            ax2.set_ylim([-0.002, max(max(freqV), max(freqGreen)) + 0.01])
        else:
            ax2.set_ylim([-0.002, max(freqGreen) + 0.01])
        ax2.yaxis.set_label_coords(-0.2, 0.5)
        
        ax3.set_ylabel('proportion')
        #ax3.set_xlabel('unique yellow (nm)')        
        #ax3.set_xlim([460, 590])
        ax3.tick_params(axis='x', colors='y')
        ax3.set_ylim([-0.005, max(max(freqBlue), max(freqYellow)) + 0.02])
        ax3.yaxis.set_label_coords(-0.2, 0.5)
        
        ax4.tick_params(axis='x', colors = 'b')
        ax3.set_xlabel('unique blue, yellow (nm)')

        #ax4.spines['bottom'].set_visible(True)
        ax3.spines['bottom'].set_visible(False)
        #ax4.set_visible(True)
        ax3.edgecolor  = 'y'
        plt.tight_layout()
        
        firsthalf = '../bps10.github.com/presentations/static/figures/'
        secondhalf = 'colorModel/uniqueHues_LMcomparison.png'
        if Volbrecht1997:
            secondhalf = 'colorModel/uniqueHues_LMcomparison_Volbrecht.png'
            ax2.legend()
            
        if savefigs:
            plt.savefig(firsthalf + secondhalf)
        plt.show()
    
    if returnVals:
        return freq, freqGreen, freqYellow, (volb['count'] / 
                                                sum(volb['count']))
        

def plotModel(plotSpecSens=False, plotCurveFamily=False,
              plotUniqueHues=False, savefigs=False, fracLvM=0.25):
    """Plot cone spectral sensitivies and first stage predictions.
    """
    
    model = colorModel()
    model.genModel(ConeRatio={'fracLvM': fracLvM, 's': 0.05, })
    FirstStage = model.returnFirstStage()   
    
    if plotSpecSens:
        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)
    
        pf.AxisFormat()
        pf.TufteAxis(ax1, ['left', 'bottom'], Nticks=[5, 5])
    
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
        ax1.set_xlabel('wavelength (nm)')
        
        plt.tight_layout()
        if savefigs:
            firsthalf = '../bps10.github.com/presentations/static/figures/'
            secondhalf = 'colorModel/specSens.png'
            plt.savefig(firsthalf + secondhalf)
        plt.show()

    if plotCurveFamily:
        SecondStage = model.returnSecondStage()
        
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


        #Not quite there. Need to be able to reference lms_Vl. Also consider
        #the center weight.
        
        from operator import itemgetter
        
        sortedlist = []
        for key in SecondStage['percent']:
            sortedlist.append(SecondStage['percent'][key])
        sortedlist = sorted(sortedlist, key=itemgetter('probSurround'), 
                            reverse=True)
        thresh = sortedlist[15]['probSurround']
        print thresh
        
        for i in SecondStage['lmsV_L']:
            if i % 2 == 0 or SecondStage['percent'][i][
                    'probSurround'] >= thresh:
                if SecondStage['percent'][i]['probSurround'] >= thresh:
                    print SecondStage['percent'][i]
                    ax1.plot(FirstStage['lambdas'], SecondStage['lmsV_M'][i],
                            c=(1,0,0), linewidth=1, alpha=0.25)
                    ax2.plot(FirstStage['lambdas'], SecondStage['lmsV_L'][i],
                            c=(0,0,1), linewidth=1, alpha=0.25)
                else:
                    ax1.plot(FirstStage['lambdas'], SecondStage['lmsV_M'][i],
                            c=(0,0,0), linewidth=1, alpha=0.10)
                    ax2.plot(FirstStage['lambdas'], SecondStage['lmsV_L'][i],
                            c=(0,0,0), linewidth=1, alpha=0.10)
                

        ax1.set_xlim([FirstStage['wavelen']['startWave'],
                      FirstStage['wavelen']['endWave']])
        ax2.set_xlim([FirstStage['wavelen']['startWave'],
                      FirstStage['wavelen']['endWave']])

        ax1.set_ylabel('sensitivity')
        ax2.set_ylabel('sensitivity')
        ax2.set_xlabel('wavelength (nm)')
        
        plt.tight_layout()
        if savefigs:
            firsthalf = '../bps10.github.com/presentations/static/figures/'
            secondhalf = 'colorModel/familyLMS_' + str(int(
                                                fracLvM * 100)) + 'L.png'
            plt.savefig(firsthalf + secondhalf)
        plt.show()

    fig = plt.figure(figsize=(8.5, 11))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    
    model.genModel(ConeRatio={'fracLvM': 0.25, 's': 0.05, })
    ThirdStage = model.returnThirdStage()
    
    pf.AxisFormat()     
    pf.TufteAxis(ax1, ['left', ], Nticks=[5, 3])
    ax1.plot(FirstStage['lambdas'], 
             np.zeros((len(FirstStage['lambdas']))), 'k', linewidth=1.0)
    ax1.plot(FirstStage['lambdas'], ThirdStage['lCenter'],
            'b', linewidth=3)
    ax1.plot(FirstStage['lambdas'], ThirdStage['mCenter'],
            'r', linewidth=3)
    ax1.set_xlim([FirstStage['wavelen']['startWave'],
                     FirstStage['wavelen']['endWave']])
    ax1.set_ylabel('activity')
    ax1.yaxis.set_label_coords(-0.2, 0.5)
    ax1.set_ylim([-0.20, 0.21])
    ax1.text(0.95, 0.95, '25% L', fontsize=16,
        horizontalalignment='right',
        verticalalignment='top',
        transform=ax1.transAxes)

    model.genModel(ConeRatio={'fracLvM': 0.5, 's': 0.05, })
    ThirdStage = model.returnThirdStage()
    
    pf.AxisFormat()     
    pf.TufteAxis(ax2, ['left', ], Nticks=[5, 3])
    ax2.plot(FirstStage['lambdas'], 
             np.zeros((len(FirstStage['lambdas']))), 'k', linewidth=1.0)
    ax2.plot(FirstStage['lambdas'], ThirdStage['lCenter'],
            'b', linewidth=3)
    ax2.plot(FirstStage['lambdas'], ThirdStage['mCenter'],
            'r', linewidth=3)
    ax2.set_xlim([FirstStage['wavelen']['startWave'],
                     FirstStage['wavelen']['endWave']])
    ax2.set_ylabel('activity')
    ax2.yaxis.set_label_coords(-0.2, 0.5)
    ax2.set_ylim([-0.20, 0.21])
    ax2.text(0.95, 0.95, '50% L', fontsize=16, 
        horizontalalignment='right',
        verticalalignment='top',
        transform=ax2.transAxes)


    model.genModel(ConeRatio={'fracLvM': 0.75, 's': 0.05, })
    ThirdStage = model.returnThirdStage()
    
    pf.AxisFormat()     
    pf.TufteAxis(ax3, ['left', 'bottom'], Nticks=[5, 3])
    ax3.plot(FirstStage['lambdas'], 
             np.zeros((len(FirstStage['lambdas']))), 'k', linewidth=1.0)
    ax3.plot(FirstStage['lambdas'], ThirdStage['lCenter'],
            'b', linewidth=3)
    ax3.plot(FirstStage['lambdas'], ThirdStage['mCenter'],
            'r', linewidth=3)
    ax3.set_xlim([FirstStage['wavelen']['startWave'],
                     FirstStage['wavelen']['endWave']])
    ax3.set_ylabel('activity')
    ax3.yaxis.set_label_coords(-0.2, 0.5)
    ax3.set_ylim([-0.20, 0.21])
    ax3.text(0.95, 0.95, '75% L', fontsize=16, 
        horizontalalignment='right',
        verticalalignment='top',
        transform=ax3.transAxes)
    ax3.set_xlabel('wavelength (nm)')
    
    plt.tight_layout()
    if savefigs:
        firsthalf = '../bps10.github.com/presentations/static/figures/'
        secondhalf = 'colorModel/PercentL.png'
        plt.savefig(firsthalf + secondhalf)
        
    plt.show()      
    
    if plotUniqueHues:
        model.findUniqueHues()
        UniqueHues = model.returnUniqueHues()
        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax1, ['left', 'bottom'], Nticks=[5, 5])

        #ax1.plot(UniqueHues['LMratio'], UniqueHues['red'],
        #        'r', linewidth=3)
        ax1.plot(UniqueHues['LMratio'], UniqueHues['green'],
                'g', linewidth=3)
        ax1.plot(UniqueHues['LMratio'], UniqueHues['blue'],
                'b', linewidth=3)
        ax1.plot(UniqueHues['LMratio'], UniqueHues['yellow'],
                'y', linewidth=3)

        ax1.set_ylabel('wavelength (nm)')
        ax1.set_xlabel('% L vs M')

        plt.tight_layout()
        if savefigs:
            firsthalf = '../bps10.github.com/presentations/static/figures/'
            secondhalf = 'colorModel/uniqueHues.png'
            plt.savefig(firsthalf + secondhalf)
        plt.show()

if __name__ == '__main__':
    #optimizeUniqueHues()
    #color = colorModel()
    LMratiosAnalysis(Volbrecht1997=True)
    #plotModel(plotSpecSens=False, plotCurveFamily=True,
    #          plotUniqueHues=False, savefigs=False)
