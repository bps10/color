# -*- coding: utf-8 *-*
from __future__ import division
import os
import numpy as np
from scipy.optimize import fsolve
from math import factorial

from base import spectsens as ss
from base import optics as op


class colorModel():
    '''
    '''
    def __init__(self, center_cones=1, q=1.300, age=None, mac_const=1.0):

        self.test = False
        self._q = q
        self.center_cones = center_cones
        self.lensMacula = getLensMaculaFilter(maxLambda=750, age=age, k=mac_const)

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
                 maxSens={'l': 559.0, 'm': 530.0, 's': 417.0, }, OD=None):
        '''
        '''
        self.findConeRatios(ConeRatio['fracLvM'], ConeRatio['s'])
        self.maxSens = maxSens
        self.fracLvM = ConeRatio['fracLvM']
        
        self.genFirstStage(OD=OD)
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

            for i in range(0, 101):

                self.findConeRatios(fracLvM=(i / 100.))
                self.genThirdStage()
                temp = self.returnThirdStage()
                RG = temp['mCenter']
                BY = temp['lCenter']

                if self.lRatio == 0:
                    uniqueGreen.append(555)
                    uniqueRed.append(592)
                else:
                    zero_cross = np.where(np.diff(np.sign(BY)))[0]
                    uniqueGreen.append(lambdas[zero_cross[0]])
                    uniqueRed.append(lambdas[np.argmin(BY)])

                if self.mRatio == 0:
                    uniqueBlue.append(467)
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
            uniqueGreen = 553
            uniqueRed = 592
        else:
            zero_cross = np.where(np.diff(np.sign(BY)))[0]
            uniqueGreen = lambdas[zero_cross[0]]
            uniqueRed = lambdas[np.argmin(BY)]

        if self.mRatio == 0:
            uniqueBlue = 474
            uniqueYellow = 577
        else:
            zero_cross = np.where(np.diff(np.sign(RG)))[0]
            uniqueBlue = lambdas[zero_cross[0]]
            try:
                uniqueYellow = lambdas[zero_cross[1]]
            except:
                uniqueYellow = np.nan 

        hues = {
            'red': uniqueRed,
            'blue': uniqueBlue,
            'green': uniqueGreen,
            'yellow': uniqueYellow,
            'LMratio': self.fracLvM,
            }
        return hues

    def genFirstStage(self, startLambda=390, endLambda=750, step=1,
                        Out='anti-log', OD=None):
        """Compute the first stage in the model
        """
        if OD is None:
            OD = [0.4, 0.38, 0.33]

        lambdas = np.arange(startLambda, endLambda + step, step)

        L_cones = ss.neitz(LambdaMax=self.maxSens['l'], LOG=False,
                            StartWavelength=startLambda,
                            OpticalDensity=OD[0],
                            EndWavelength=endLambda, 
                            resolution=step)
        L_cones /= self.lensMacula
        
        M_cones = ss.neitz(LambdaMax=self.maxSens['m'], LOG=False,
                            StartWavelength=startLambda,
                            OpticalDensity=OD[1],
                            EndWavelength=endLambda, 
                            resolution=step)
        M_cones /= self.lensMacula
        
        S_cones = ss.neitz(LambdaMax=self.maxSens['s'], LOG=False,
                            StartWavelength=startLambda,
                            OpticalDensity=OD[2],
                            EndWavelength=endLambda, 
                            resolution=step)
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

        s, _m, _l, i = self.sRatio * 100, 0, 0, 0
        for m in range(0, 101):
            for l in range(0, 101): 
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
        '''Compute the third stage in the model
        '''

        percentL = self.lRatio / (self.mRatio + self.lRatio)
        percentM = self.mRatio / (self.mRatio + self.lRatio)

        self.ThirdStage = {
            'mCenter': np.zeros(len(self.FirstStage['lambdas'])),
            'lCenter': np.zeros(len(self.FirstStage['lambdas'])),
            }

        for i in self.SecondStage['percent']:

            # surround probabilities:
            lNum = self.SecondStage['percent'][i]['l']
            mNum = self.SecondStage['percent'][i]['m']
            probSur = binom(lNum, lNum + mNum, percentL)
            self.SecondStage['percent'][i]['probSurround'] = probSur

            for num_L, lCenter in self.SecondStage['lmsV_L'][i].iteritems():

                centerProb = binom(num_L, self.center_cones, percentL)
                self.ThirdStage['lCenter'] += (lCenter * probSur) 

            for num_M, mCenter in self.SecondStage['lmsV_M'][i].iteritems():

                centerProb = binom(num_M, self.center_cones, percentM)
                self.ThirdStage['mCenter'] += (mCenter * probSur) 

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


def getLensMaculaFilter(maxLambda=750, age=None, k=1.0):
    '''
    '''
    if age is not None:
        macula = op.filters.stockman(minLambda=390, 
            maxLambda=maxLambda, RETURN_SPECTRUM=False, 
            ONLY_MACULA=True, resolution=1)
        spectrum = np.arange(390, maxLambda + 1, 1)
        lens = op.filters.lens_age_correction(age, spectrum)
        
    else:
        macula = op.filters.stockman(minLambda=390, 
            maxLambda=maxLambda, RETURN_SPECTRUM=False, 
            ONLY_MACULA=True, resolution=1)
        lens = op.filters.stockman(minLambda=390, 
            maxLambda=maxLambda, RETURN_SPECTRUM=False, 
            ONLY_LENS=True, resolution=1)
        
    filters = 10.0 ** (lens + (k * macula))
    return filters


def getCarroll_LMratios():
    '''
    '''
    temp = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
        'data/Carroll2002_lmRatios.txt')
    return np.genfromtxt(temp, 
                delimiter='\t', dtype=None, skip_header=0, 
                names=True)

def getHueScalingData(ConeRatio, maxSens, scale=False, norm=False):
    '''
    '''
    model = colorModel()
    model.genModel(ConeRatio=ConeRatio, maxSens=maxSens)
    # get valence data: returned in OFF config (red+blue)
    FirstStage = model.returnFirstStage()
    ThirdStage = model.returnThirdStage()
    N = len(ThirdStage['lCenter'])
    red = ThirdStage['mCenter'].clip(min=0)
    green = (-1 * ThirdStage['mCenter']).clip(min=0)
    blue = ThirdStage['lCenter'].clip(min=0)
    yellow = (-1 * ThirdStage['lCenter']).clip(min=0)

    #scale to percentage
    if scale and not norm:
        total = red + green + blue + yellow
        red = red / total * 100
        green = green / total * 100
        blue = blue / total * 100
        yellow = yellow / total * 100
    if norm and not scale:
        red /= np.max(red)
        green /= np.max(green)
        blue /= np.max(blue)
        yellow /= np.max(yellow)

    hues = {'red': red, 'green': green,
            'blue': blue, 'yellow': yellow, 
            'lambdas': FirstStage['lambdas']}
    return hues


def optimizeChannel(q, cones, percent, Center):
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

