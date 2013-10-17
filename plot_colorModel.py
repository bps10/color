#! /usr/bin/env python
# -*- coding: utf-8 *-*
from __future__ import division
import os
import numpy as np
import matplotlib.pylab as plt
from math import factorial

import colorModel as cm
from base import plot as pf


def binomPlot():
    '''
    '''
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)    
    pf.AxisFormat(linewidth=2, markersize=14)
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])

    for probL in range(1, 10):
        probM = 10 - probL
        dist = []
        for percentL in range(0, 101):

            dist.append(cm.binom(percentL, 100, probL / 10))

            color = [probL / 10, probM / 10, 0]
        ax.plot(np.arange(0,101), dist, c=color)

    ax.set_xlabel("%L v M")
    ax.set_ylabel("probability")
    plt.tight_layout()
    plt.show()


def eccentricityAnalysis():
    '''
    '''
    cond = {0: {'percent': 0.40, 'lines': '-'},
                  1: {'percent': 0.60, 'lines': '--'},
                  2: {'percent': 0.80, 'lines': ':'}, }

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)    
    pf.AxisFormat(linewidth=2, markersize=14)
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])

    for c in cond.itervalues():

        uniqueHues = {}
        for center_cones in range(1, 6):

            model = cm.colorModel(q=1.3, center_cones=center_cones)
            model.genModel(ConeRatio={'fracLvM': c['percent'], 's': 0.05, })
            uniqueHues[center_cones] = model.get_current_uniqueHues()

        yellow, green, blue, center_cones = [], [], [], []
        for _cones, hues in uniqueHues.iteritems():
            yellow.append(hues['yellow'])
            green.append(hues['green'])
            blue.append(hues['blue'])
            center_cones.append(_cones)

        ax.plot(center_cones, green, 'go' + c['lines'])
        ax.plot(center_cones, blue, 'bo' + c['lines'])
        ax.plot(center_cones, yellow, 'yo' + c['lines'])

    #ax.set_xlim([0.85, 5.15])
    ax.set_ylim([460, 625])
    ax.set_xlabel("number of center cones")
    ax.set_ylabel("wavelength (nm)")
    plt.tight_layout()
    plt.show()


def MetaAnalysis():
    '''
    '''
    dat = np.genfromtxt('data/hueMeta_Kuehni.txt', delimiter='\t',
        names=True, usecols= (1, 2, 3, 4, 5, 6, 7, 8 , 9))
    for name in dat.dtype.names:
        if name[-1] is not 'N':
            if name[:6] == 'yellow':
                n = 'yellow'
            if name[:5] == 'green':
                n = 'green'
            else:
                n = 'blue'
            total, weight = 0, 0
            for i, study in enumerate(dat[name]):

                if not np.isnan(study):
                    total += dat[n + '_N'][i]
                    if name[-3:] == 'std':
                        weight += (dat[n + '_N'][i] - 1) * (study ** 2)
                    else: 
                        weight += dat[n + '_N'][i] * study

            if name[-3:] == 'std':
                print name + ': ', np.sqrt(weight / (total - 1)), total
            else:
                print name + ': ', weight / total, total


def HueScaling(lPeak=559):
    '''
    '''
    hues = cm.getHueScalingData(
                ConeRatio={'fracLvM': 0.70, 's': 0.05, },
                maxSens={'l': lPeak, 'm': 530.0, 's': 417.0, })
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    pf.AxisFormat(linewidth=3)
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])

    ax.plot(hues['lambdas'], hues['red'], 'r')
    ax.plot(hues['lambdas'], hues['green'], 'g') 
    ax.plot(hues['lambdas'], hues['blue'], 'b') 
    ax.plot(hues['lambdas'], hues['yellow'], 'y')  

    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('percentage')
    ax.set_xlim([390, 750])

    plt.tight_layout()
    plt.show()


def matching(lPeak=559, test=420, match1=485, match2=680):
    '''
    '''
    ratios = np.arange(1, 100) / 100
    prop_match = []
    for ratio in ratios:
        hues = cm.getHueScalingData(
                    ConeRatio={'fracLvM': ratio, 's': 0.05, },
                    maxSens={'l': lPeak, 'm': 530.0, 's': 417.0, })

        test_light = np.where(hues['lambdas'] == test)
        match_1_ind = np.where(hues['lambdas'] == match1)
        match_2_ind = np.where(hues['lambdas'] == match2)

        light = np.array([
            float(hues['blue'][test_light] - hues['yellow'][test_light]),
            float(hues['red'][test_light] - hues['green'][test_light])])

        a1 = float(hues['blue'][match_1_ind] - hues['yellow'][match_1_ind])
        a2 = float(hues['blue'][match_2_ind] - hues['yellow'][match_2_ind])
        b1 = float(hues['red'][match_1_ind] - hues['green'][match_1_ind])
        b2 = float(hues['red'][match_2_ind] - hues['green'][match_2_ind])

        system = np.array([[a1, a2],[b1, b2]])

        match = np.dot(np.linalg.inv(system), light)

        prop_match.append(match[0] / (match.sum()))
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    pf.AxisFormat(linewidth=3)
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])

    ax.plot(ratios * 100, prop_match)

    ax.set_xlabel('%L')
    ax.set_ylabel(('proportion (' + str(match1) + 
        ' / ' + '(' + str(match1) + '+' + str(match2) + ')'))

    plt.tight_layout()
    plt.show()


def match_v2(lPeak=559, test=420, match1=485, match2=680):
    '''
    '''
    get_ind = lambda match: np.where(hues['lambdas'] == match)
    a1 = lambda match: float(hues['blue'][match] - hues['yellow'][match])
    b1 = lambda match: float(hues['red'][match] - hues['green'][match])

    ratios = np.arange(1, 100) / 100
    match = []
    for ratio in ratios:
        hues = cm.getHueScalingData(
                    ConeRatio={'fracLvM': ratio, 's': 0.05, },
                    maxSens={'l': lPeak, 'm': 530.0, 's': 417.0, })

        test_light = get_ind(test)
        match_2_ind = get_ind(match2)

        by_light = float(hues['blue'][test_light] - hues['yellow'][test_light])
        rg_light = float(hues['red'][test_light] - hues['green'][test_light])

        
        a2 = float(hues['blue'][match_2_ind] - hues['yellow'][match_2_ind])      
        b2 = float(hues['red'][match_2_ind] - hues['green'][match_2_ind])

        sol1 = by_light - (0.1 * a2)
        sol2 = rg_light - (0.9 * b2)
        print sol1, sol2

        solution = False
        wvlen = 430
        while solution is not True:
            ind = get_ind(wvlen)
            eq1 = (0.1 * a1(ind)) - by_light
            eq2 = (0.9 * b1(ind)) - rg_light
            print a1(ind), b1(ind), eq1, eq2
            if eq1 + eq2 > 0 and eq1 + eq2 < 0.01:
                solution = True
            elif wvlen < 500:
                wvlen += 1
            else:
                wvlen = 0
                solution = True

        ## Just find the index solution and then interpolate

        match.append(wvlen)
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    pf.AxisFormat(linewidth=3)
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])

    ax.plot(ratios * 100, match)

    ax.set_xlabel('%L')
    ax.set_ylabel('wavelength')

    plt.tight_layout()
    plt.show()


def LMratiosAnalysis(Volbrecht1997=True, returnVals=False, 
                        plot=True, savefigs=False):
    '''
    '''

    model = cm.colorModel(q=1.300)
    
    model.genModel()
    #model.findUniqueHues()
    #uniqueHues = model.returnUniqueHues()
    carroll = cm.getCarroll_LMratios()

    green, yellow, blue = [], [], []
    for i, subject in enumerate(carroll['L']):

        model.genModel(
            maxSens={'l': carroll['lPeak'][i], 'm': 529.0, 's': 417.0, },
            ConeRatio={'fracLvM': carroll['L'][i] / 100.0, 's': 0.05, })
        uniqueHues = model.get_current_uniqueHues()

        green.append(uniqueHues['green'])
        yellow.append(uniqueHues['yellow'])
        blue.append(uniqueHues['blue'])

    print 'hue | \t mean | \t stdv'
    print 'green: ', np.mean(green), np.std(green, ddof=1)
    print 'yellow: ', np.mean(yellow), np.std(yellow, ddof=1)
    print 'blue: ', np.mean(blue), np.std(blue, ddof=1)   
                       
    BINS = np.arange(0, 101, 6)
    if Volbrecht1997:
        BINS_G = np.arange(488, 564, 3)
        temp = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
            'data/Volbrecht1997.txt')
        volb = np.genfromtxt(temp, delimiter='\t',
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
        
        
        BINS, freq = pf.histOutline(freq / np.sum(freq), BINS)
        BINS_Gout, freqGreen = pf.histOutline(freqGreen / np.sum(freqGreen), 
                                           BINS_G)
        BINS_Yout, freqYellow = pf.histOutline(
                                freqYellow / np.sum(freqYellow), BINS_Y)
        BINS_Bout, freqBlue = pf.histOutline(freqBlue / np.sum(freqBlue), 
                                            BINS_B)

        ax1.plot(BINS, freq, 'k', linewidth=3)
        
        if Volbrecht1997:
            binV, freqV = pf.histOutline(volb['025deg'] / np.sum(
                                        volb['025deg']), BINS_G)
            ax2.plot(binV, freqV, c='0.8', linewidth=3,
                             label='Volbrecht 1997')
            ax2.fill_between(binV, freqV, 0, color='0.8')
            
        ax2.plot(BINS_Gout, freqGreen, 'g', linewidth=3,
                 label='predicted')
        ax3.plot(BINS_Yout, freqYellow, 'y', linewidth=3) 
        ax4.plot(BINS_Bout, freqBlue, 'b', linewidth=3) 
          
        ax1.set_xlim([0, 100])
        ax1.set_ylim([-0.002, np.max(freq) + 0.01])
        ax1.set_ylabel('proportion')
        ax1.set_xlabel('% L v M')
        ax1.yaxis.set_label_coords(-0.2, 0.5)
        
        ax2.set_ylabel('proportion')
        ax2.set_xlabel('unique green (nm)')
        ax2.set_xlim([490, 560])
        ax3.set_xlim([514, 590])
        ax4.set_xlim([460, 536])
        if Volbrecht1997:
            ax2.set_ylim([-0.002, np.max(np.max(freqV), np.max(freqGreen)) + 0.01])
        else:
            ax2.set_ylim([-0.002, np.max(freqGreen) + 0.01])
        ax2.yaxis.set_label_coords(-0.2, 0.5)
        
        ax3.set_ylabel('proportion')
        #ax3.set_xlabel('unique yellow (nm)')        
        #ax3.set_xlim([460, 590])
        ax3.tick_params(axis='x', colors='y')
        ax3.set_ylim([-0.005, np.max(np.max(freqBlue), np.max(freqYellow)) + 0.02])
        ax3.yaxis.set_label_coords(-0.2, 0.5)
        
        ax4.tick_params(axis='x', colors = 'b')
        ax3.set_xlabel('unique blue, yellow (nm)')

        #ax4.spines['bottom'].set_visible(True)
        ax3.spines['bottom'].set_visible(False)
        #ax4.set_visible(True)
        ax3.edgecolor  = 'y'
        plt.tight_layout()

        if Volbrecht1997:
            savename = 'uniqueHues_LMcomparison_Volbrecht.eps'
            ax2.legend()
            
        if savefigs:
            plt.savefig(savename)
        plt.show()
    
    if returnVals:
        return freq, freqGreen, freqYellow, (volb['count'] / 
                                                np.sum(volb['count']))
        

def plotModel(plotModel=True, plotCurveFamily=False,
              plotUniqueHues=False, savefigs=False, fracLvM=0.25,
              SHOW=False, age=None, maxSens=None, OD=None):
    """Plot cone spectral sensitivies and first stage predictions.
    """
    
    if plotCurveFamily:
        if maxSens is None:
            maxSens = {'l': 559, 'm': 529.0, 's': 421.0, }
        model = cm.colorModel(age=age)
        model.genModel(ConeRatio={'fracLvM': fracLvM, 's': 0.05, },
            maxSens=maxSens, OD=OD)

        FirstStage = model.returnFirstStage()   
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
            #print SecondStage['percent'][key]
        sortedlist = sorted(sortedlist, key=itemgetter('probSurround'), 
                            reverse=True)
        thresh = sortedlist[15]['probSurround']

        for i in SecondStage['lmsV_L']:
            if i % 2 == 0 or SecondStage['percent'][i][
                    'probSurround'] >= thresh:
                if SecondStage['percent'][i]['probSurround'] >= thresh:
                    print SecondStage['percent'][i]
                    ax1.plot(FirstStage['lambdas'], 
                            SecondStage['lmsV_M'][i][1],
                            c=(1,0,0), linewidth=1, alpha=0.25)
                    ax2.plot(FirstStage['lambdas'], 
                            SecondStage['lmsV_L'][i][1],
                            c=(0,0,1), linewidth=1, alpha=0.25)
                else:
                    ax1.plot(FirstStage['lambdas'], 
                            SecondStage['lmsV_M'][i][1],
                            c=(0,0,0), linewidth=1, alpha=0.10)
                    ax2.plot(FirstStage['lambdas'], 
                            SecondStage['lmsV_L'][i][1],
                            c=(0,0,0), linewidth=1, alpha=0.10)
                
        ax1.set_ylim([-0.4, 0.4])

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

    if plotModel:
        model = cm.colorModel(age=age)
        model.genModel(ConeRatio={'fracLvM': 0.25, 's': 0.05, },
            maxSens=maxSens, OD=OD)

        FirstStage = model.returnFirstStage() 
        ThirdStage = model.returnThirdStage()  

        fig = plt.figure(figsize=(8.5, 11))
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        
        
        
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
        #ax1.set_ylabel('activity')
        ax1.yaxis.set_label_coords(-0.2, 0.5)
        ax1.set_ylim([-0.20, 0.21])
        ax1.text(0.95, 0.95, '25% L', fontsize=16,
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax1.transAxes)

        model.genModel(ConeRatio={'fracLvM': 0.5, 's': 0.05, },
            maxSens=maxSens, OD=OD)
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
        ax2.set_ylabel('sensitivity')
        ax2.yaxis.set_label_coords(-0.2, 0.5)
        ax2.set_ylim([-0.20, 0.21])
        ax2.text(0.95, 0.95, '50% L', fontsize=16, 
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax2.transAxes)


        model.genModel(ConeRatio={'fracLvM': 0.75, 's': 0.05, },
            maxSens=maxSens, OD=OD)
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
        #ax3.set_ylabel('activity')
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
        model = cm.colorModel(age=age)

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        pf.AxisFormat(linewidth=3)
        pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[4, 5])

        style = ['-', '--', '-.']
        i = 0
        for lPeak in [559.0, 557.25, 555.5]:

            model.genModel(ConeRatio={'fracLvM': 0.25, 's': 0.05, },
                maxSens={'l': lPeak, 'm': 529.0, 's': 419.0, }, OD=OD)
            model.findUniqueHues()

            UniqueHues = model.returnUniqueHues()

            ax.plot(UniqueHues['LMratio'], UniqueHues['green'],
                    'g' + style[i], label=str(int(lPeak)))
            ax.plot(UniqueHues['LMratio'], UniqueHues['blue'],
                    'b' + style[i], label=str(int(lPeak)))
            ax.plot(UniqueHues['LMratio'], UniqueHues['yellow'],
                    'y' + style[i], label=str(int(lPeak)))
            i += 1

        ax.set_xlim([20, 100])
        ax.set_ylim([460, 600])
        ax.set_ylabel('wavelength (nm)')
        ax.set_xlabel('% L vs M')

        plt.tight_layout()
        if savefigs:
            firsthalf = '../bps10.github.com/presentations/static/figures/'
            secondhalf = 'colorModel/uniqueHues.png'
            plt.savefig(firsthalf + secondhalf)
        if SHOW:
            plt.show()
        else:
            return ax

def main(args):
    '''
    '''
    if args.LM < 1 and args.LM > 0:
        LMratio = args.LM
    else:
        raise ValueError('LM ratio must be between 0 and 1')

    if args.binom:
        binomPlot()
    
    if args.eccen:
        eccentricityAnalysis()
    
    if args.ratio:
        LMratiosAnalysis(Volbrecht1997=True, savefigs=args.save)

    if args.hues:
        HueScaling()

    if args.match:
        #matching()
        match_v2()

    if args.kuehni:
        MetaAnalysis()

    plotModel( 
            plotModel=args.model,
            plotCurveFamily=args.curve,
            plotUniqueHues=args.unique, 
            savefigs=args.save,
            fracLvM=LMratio)

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description="Color Model: display the \
        Neitz color model")
    
    parser.add_argument("-b", "--binom", action="store_true", 
                        help="plot a series of binomial distributions")
    parser.add_argument("-e", "--eccen", action="store_true",
                        help="plot unique hue as a function of eccentricity")
    parser.add_argument("-r", "--ratio", action="store_true", 
                        help="plot unique hue as a function of LM ratio")
    parser.add_argument("-m", "--model", action="store_true",
                        help="plot the neitz color model")
    parser.add_argument("-c", "--curve", action="store_true",
                        help="plot a family of valence curves")
    parser.add_argument("-u", "--unique", action="store_true",
                        help="plot unique hues")
    parser.add_argument("-q", "--hues", action="store_true",
                        help="plot hue scaling data")
    parser.add_argument("-t", "--match", action="store_true",
                        help="plot moreland match")

    parser.add_argument("-s", "--save", action="store_true",
                        help="save plots - not working right now")
    parser.add_argument("--LM", type=float, default=0.25,
                        help="set L:M ratio, default=0.25")
    parser.add_argument("-k", "--kuehni", action='store_true',
                        help="print Kuehni meta analysis of unique hues")
    # add default save location

    args = parser.parse_args()
    main(args)