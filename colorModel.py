# -*- coding: utf-8 *-*
from __future__ import division
import numpy as np
from scipy.optimize import fsolve
from math import factorial

from spectsens import spectsens


class colorModel():
    '''
    '''
    def __init__(self, center_cones=1, q=1.300):

        self.test = False
        self.step = 1
        self._q = q
        self.center_cones = center_cones
        self.lensMacula = getStockmanFilter()

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
        self.fracLvM = ConeRatio['fracLvM']
        
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

    def get_current_uniqueHues(self):
        '''
        '''
        lambdas = self.FirstStage['lambdas']
        RG = self.ThirdStage['mCenter']
        BY = self.ThirdStage['lCenter']

        if self.lRatio == 0:
            uniqueGreen = 555
            uniqueRed = 592
        else:
            zero_cross = np.where(np.diff(np.sign(BY)))[0]
            uniqueGreen = lambdas[zero_cross[0]]
            uniqueRed = lambdas[np.argmin(BY)]

        if self.mRatio == 100:
            uniqueBlue = 474
            uniqueYellow = 575
        else:
            zero_cross = np.where(np.diff(np.sign(RG)))[0]
            uniqueBlue = lambdas[zero_cross[0]]
            try:
                uniqueYellow = lambdas[zero_cross[1]]
            except:
                uniqueYellow = 600 

        hues = {
            'red': uniqueRed,
            'blue': uniqueBlue,
            'green': uniqueGreen,
            'yellow': uniqueYellow,
            'LMratio': self.fracLvM,
            }
        return hues

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
        # get cones into dict for optimization below
        L_cones = self.FirstStage['L_cones']
        M_cones = self.FirstStage['M_cones']
        cones = {
            's': self.FirstStage['S_cones'],
            'm': self.FirstStage['M_cones'],
            'l': self.FirstStage['L_cones'], }

        self.SecondStage = {'lmsV_L': {}, 
                            'lmsV_M': {}, 
                            'percent': {}, 
                            'center_cones': self.center_cones, }

        s, _m, _l, i = 5, 0, 0, 0
        for m in range(0, 101, self.step):
            for l in range(0, 101, self.step): 
                if (l + m + s) == 100:
                    percent = {'s': s, 'm': m, 'l': l, }
                    
                    self.SecondStage['lmsV_L'][_l] = {}
                    self.SecondStage['lmsV_M'][_m] = {}
                    for Lcenter in range(0, self.center_cones + 1):

                        Mcenter = self.center_cones - Lcenter
                        propL = Lcenter / self.center_cones
                        propM = Mcenter / self.center_cones
                        center_cones = (L_cones * propL) + (M_cones * propM)

                        if propL > 0.5:

                            lmsV_L = optimizeChannel(self._q,
                                                    cones, percent,
                                                    Center=center_cones)

                            self.SecondStage['lmsV_L'][_l][Lcenter] = lmsV_L

                        if propM > 0.5:

                            lmsV_M = optimizeChannel(self._q, 
                                                    cones, percent,
                                                    Center=center_cones)
                            
                            self.SecondStage['lmsV_M'][_m][Mcenter] = lmsV_M
                            

                    self.SecondStage['percent'][i] = percent
                    _m += 1
                    _l += 1
                    i += 1
                        

    def genThirdStage(self):
        """Compute the third stage in the model
        """

        percentL = self.lRatio / (self.mRatio + self.lRatio)
        percentM = self.mRatio / (self.mRatio + self.lRatio)
        
        t = self.SecondStage['lmsV_L'][0].keys()
        self.ThirdStage = {
            'mCenter': np.zeros(len(self.SecondStage['lmsV_M'][0][t[0]])),
            'lCenter': np.zeros(len(self.SecondStage['lmsV_L'][0][t[0]])),
            }

        for i in self.SecondStage['percent']:

            # surround probabilities:
            lNum = self.SecondStage['percent'][i]['l']
            mNum = self.SecondStage['percent'][i]['m']
            probSur = binom(lNum, lNum + mNum, percentL)                            
            self.SecondStage['percent'][i]['probSurround'] = probSur

            
            for num_L, lCenter in self.SecondStage['lmsV_L'][i].iteritems():

                centerProb = binom(num_L, self.center_cones, percentL)
                self.ThirdStage['lCenter'] += (lCenter * 
                                                centerProb * probSur) 

            for num_M, mCenter in self.SecondStage['lmsV_M'][i].iteritems():

                centerProb = binom(num_M, self.center_cones, percentM)
                self.ThirdStage['mCenter'] += (mCenter * 
                                                centerProb * probSur) 

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


def getStockmanFilter(maxLambda=750):
    '''
    '''
    lens = np.genfromtxt('static/stockman/lens.csv', delimiter=',')[::10]
    macula = np.genfromtxt('static/stockman/macular.csv', 
                           delimiter=',')[::10]
    spectrum = lens[:, 0]
    ind = np.where(spectrum == maxLambda)[0] + 1
                #just take upto a given index (750nm)
    return 10 ** (lens[:ind, 1] + macula[:ind, 1])


def getCarroll_LMratios():
    '''
    '''
    return np.genfromtxt('static/data/Carroll2002_lmRatios.txt', 
                delimiter='\t', dtype=None, skip_header=0, 
                names=True)


def optimizeChannel(q, cones, percent, Center, test=False):
    '''
    '''
    m_ = percent['m'] / (percent['m'] + percent['l'])
    l_ = percent['l'] / (percent['m'] + percent['l'])
    fun = lambda w, Center: (w * (q * cones['s'] +
                                m_ * cones['m'] +
                                l_ * cones['l']) -
                             Center)

    # error function to minimize
    err = lambda w, Center: (fun(w, Center)).sum()
    w = fsolve(err, 1, args=(Center))
    out = fun(w, Center)

    if test:
        temp = err(w, Center)
        if temp > 1e-8:
            print percent
            raise ValueError('error function not minimized properly')

    return out


def binom(k, n, p):
    '''
    '''
    return (factorial(n) / (factorial(k) * factorial(n - k))
                        * (p ** k) * (1 - p) ** (n - k))


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



if __name__ == '__main__':

    color = colorModel()

