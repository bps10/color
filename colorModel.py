# -*- coding: utf-8 *-*
from __future__ import division
import os
import numpy as np
from scipy.optimize import fsolve
from math import factorial

from base import spectsens as ss
from base import optics as op


class colorModel():
    '''This class creates the color model described in Schmidt, Neitz, Neitz 2014 in
    Journal of the Optical Society of America A.

    :param center_cones: Change the number of cones assumed in the center \
    of the ganglion cells. This is used during eccentricity analysis. \
    Default = 1
    :type center_cones: float
    :param q: This is the constant that the S cones are multiplied by. \
    It should not be changed as we assume it comes from the ratio of \
    S to (M+L) peak sensitivities (thermal noise).
    :type q: float
    :param age: change the age of the individual. This will be incorporated \
    by the lens filter, which is known to yellow with age. If not set \
    the Stockman lens filter is used. Otherwise, the age parameterizable \
    Pokorny and Smith measurements are used.
    :type age: float
    :param mac_const: Change the magnitude of the macula. This constant \
    will scale the density of the macula. This is used to adjust macula \
    based on individual measurments.

    '''
    def __init__(self, center_cones=1, q=1.300, age=None, mac_const=1.0):
        '''
        '''
        self.test = False
        self._q = q
        self.center_cones = center_cones
        self.lensMacula = getLensMaculaFilter(maxLambda=750, age=age, k=mac_const)

    def findConeRatios(self, fracLvM, fracS=None):
        '''Convert a fraction L / (L+M) into proper L, M, S proportions that sum to 1.

        :param fracLvM: fraction of L:M cones (L / (L+M)). This value must be\
        between 0 and 1.
        :type fracLvM: float
        :param fracS: specify the fraction of S cones in the retina. Default is\
        None and 0.5 is assumed.
        :type fracS: float
        :returns: None. This sets the internal values for self.lRatio, self.mRatio\
        and self.sRatio.

        This method gets called by genModel()
        '''
        if fracS > 1 or fracLvM < 0:
            raise IOError('Fraction of LvM must be between 0 and 1!')

        if fracS is not None:
            self.sRatio = fracS
        self.lRatio = (1 - self.sRatio) * (fracLvM)
        self.mRatio = (1 - self.sRatio) * (1 - fracLvM)

    def genModel(self, ConeRatio={'fracLvM': 0.75, 's': 0.05, },
                 maxSens={'l': 559.0, 'm': 530.0, 's': 417.0, }, OD=None):
        '''This function generates the model. It can be called repeatedly\
        to make changes to the parameters.

        :param ConeRatio: Set the cone ratios with fracLvM and S. Default is\
        {'fracLvM': 0.75, 's': 0.05, }.
        :type ConeRatio: dict
        :param maxSens: Set the peak sensitivity of the L, M and S cones in\
        nm. Default is {'l': 559.0, 'm': 530.0, 's': 417.0, }.
        :type maxSens: dict
        :param OD: Set the optical density assumed for the L, M, and S cones.\
        If not set default values of [0.4, 0.38, 0.33] are assumed.
        :type OD: list
        :returns: None. This method instantiates the model.
        '''
        self.findConeRatios(ConeRatio['fracLvM'], ConeRatio['s'])
        self.maxSens = maxSens
        self.fracLvM = ConeRatio['fracLvM']
        
        self.genFirstStage(OD=OD)
        self.genSecondStage()
        self.genThirdStage()

    def findUniqueHues(self):
        '''Find the unique hues for a given set of parameters. This \
        method will iterate through every value of %L from 0 to 100.

        This method does not return anything. To return the results \
        from this analysis use returnUniqueHues().

        '''
        lambdas = self.FirstStage['lambdas']
        uniqueGreen, uniqueBlue, uniqueYellow = [], [], []
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
                else:
                    zero_cross = np.where(np.diff(np.sign(BY)))[0]
                    uniqueGreen.append(lambdas[zero_cross[0]])

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
            'blue': uniqueBlue,
            'green': uniqueGreen,
            'yellow': uniqueYellow,
            'LMratio': LMratio,
            }

    def get_current_uniqueHues(self):
        '''Find the unique hues for only the specified L:M ratio.

        :returns: a dictionary of the unique hues.
        '''
        lambdas = self.FirstStage['lambdas']
        RG = self.ThirdStage['mCenter']
        BY = self.ThirdStage['lCenter']

        if self.lRatio == 0:
            uniqueGreen = 553
        else:
            zero_cross = np.where(np.diff(np.sign(BY)))[0]
            uniqueGreen = lambdas[zero_cross[0]]

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
            'blue': uniqueBlue,
            'green': uniqueGreen,
            'yellow': uniqueYellow,
            'LMratio': self.fracLvM,
            }
        return hues

    def genFirstStage(self, startLambda=390, endLambda=750, step=1,
                        Out='anti-log', OD=None):
        '''Compute the first stage in the model.

        '''
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
        '''Compute the second stage in the model
        '''
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
                self.ThirdStage['lCenter'] += (lCenter * centerProb * probSur) 

            for num_M, mCenter in self.SecondStage['lmsV_M'][i].iteritems():

                centerProb = binom(num_M, self.center_cones, percentM)
                self.ThirdStage['mCenter'] += (mCenter * centerProb * probSur) 

    def returnFirstStage(self):
        '''Return the data generated by the first stage of the model.

        :returns: Dictionary that contains L, M, S spectral sensitivities,\
        the spectrum and the sampling parameters.
        '''
        return self.FirstStage

    def returnSecondStage(self):
        '''Return the data generated by the second stage of the model.

        :returns: Dictionary that contains lmvV_M and lmsV_L data.
        '''
        return self.SecondStage

    def returnThirdStage(self):
        '''Return the data generated by the third stage of the model.

        :returns: Dictionary that contains the predicted sensitivities of\
        the mCenter and lCenter midget ganglion cells.
        '''
        return self.ThirdStage

    def returnUniqueHues(self):
        '''Return the predicted location of unique hues as a function of\
        L:M cone ratio. The OD, and peak L sensitivies set during genModel()\
        are used, but L:M cone ratio is changed iteratively from 0 to 100.

        :returns: Dictionary containing predicted values for blue, yellow \
        and green.
        '''
        return self.uniqueHues


def getLensMaculaFilter(maxLambda=750, age=None, k=1.0):
    '''Get the lens and macula filters. This function is called\
    during the creation of the class.
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
    '''This function loads and returns the L:M cone rations measured \
    by Joe Carroll in 2002.

    :returns: array with data.
    :rtype: numpy.array
    '''
    temp = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
        'data/Carroll2002_lmRatios.txt')
    return np.genfromtxt(temp, 
                delimiter='\t', dtype=None, skip_header=0, 
                names=True)

def getHueScalingData(ConeRatio, maxSens, scale=False, norm=False):
    '''Return half wave rectified color channels.
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
    '''An internal function that solves for the constant that makes the \
    center and surround of each ganglion cell balanced (i.e. integrate to 0).

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
    '''A binomial function.
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

