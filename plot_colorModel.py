#! /usr/bin/env python
# -*- coding: utf-8 *-*
from __future__ import division
import os
import numpy as np
import matplotlib.pylab as plt
from math import factorial
from operator import itemgetter

from base import plot as pf


def binomPlot(cm):
    '''
    '''
    fig = plt.figure(figsize=(8, 6))
    fig.set_tight_layout(True)
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

    ax.set_xlabel("%L")
    ax.set_ylabel("probability")
    plt.show()


def eccentricityAnalysis(cm):
    '''
    '''
    cond = {0: {'percent': 0.40, 'lines': '-'},
            1: {'percent': 0.60, 'lines': '--'},
            2: {'percent': 0.80, 'lines': ':'}, }

    fig = plt.figure(figsize=(8, 6))
    fig.set_tight_layout(True)
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

    ax.set_ylim([460, 625])
    ax.set_xlabel("number of center cones")
    ax.set_ylabel("wavelength (nm)")
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


def HueScaling(cm, lPeak=559):
    '''
    '''
    hues = cm.getHueScalingData(
                ConeRatio={'fracLvM': 0.70, 's': 0.05, },
                maxSens={'l': lPeak, 'm': 530.0, 's': 417.0, })
    
    ax = pf.get_axes()[0]
    ax.plot(hues['lambdas'], hues['red'], 'r')
    ax.plot(hues['lambdas'], hues['green'], 'g') 
    ax.plot(hues['lambdas'], hues['blue'], 'b') 
    ax.plot(hues['lambdas'], hues['yellow'], 'y')  

    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('percentage')
    ax.set_xlim([390, 750])

    plt.show()


def LMratiosAnalysis(cm, Volbrecht1997=True, returnVals=False, 
                        plot=True, savefigs=False):
    '''
    '''

    model = cm.colorModel()
    
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
        fig.set_tight_layout(True)
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
            ax2.set_ylim([-0.002, np.max(np.max(freqV), 
                                         np.max(freqGreen)) + 0.01])
        else:
            ax2.set_ylim([-0.002, np.max(freqGreen) + 0.01])
        ax2.yaxis.set_label_coords(-0.2, 0.5)
        
        ax3.set_ylabel('proportion')
        ax3.tick_params(axis='x', colors='y')
        ax3.set_ylim([-0.005, np.max(np.max(freqBlue), 
                                     np.max(freqYellow)) + 0.02])
        ax3.yaxis.set_label_coords(-0.2, 0.5)
        
        ax4.tick_params(axis='x', colors = 'b')
        ax3.set_xlabel('unique blue, yellow (nm)')

        ax3.spines['bottom'].set_visible(False)
        ax3.edgecolor  = 'y'

        if Volbrecht1997:
            ax2.legend()
            
        if savefigs:
            savename = 'uniqueHues_LMcomparison_Volbrecht.eps'
            plt.savefig(savename)
        plt.show()
    
    if returnVals:
        return freq, freqGreen, freqYellow, (volb['count'] / 
                                                np.sum(volb['count']))
        

def plotModel(cm, plotModel=True, plotCurveFamily=False,
              plotUniqueHues=False, savefigs=False, 
              fracLvM=0.25, SHOW=True, OD=None, age=None,
              maxSens=None):
    """Plot cone spectral sensitivies and first stage predictions.
    """
    if maxSens is None:
        maxSens = {'l': 559.0, 'm': 530.0, 's': 417.0, }

    if plotCurveFamily:
        model = cm.colorModel(age=age)
        model.genModel(ConeRatio={'fracLvM': fracLvM, 's': 0.05, },
            maxSens=maxSens, OD=OD)

        FirstStage = model.returnFirstStage()   
        SecondStage = model.returnSecondStage()
        
        fig = plt.figure(figsize=(8.5, 8))
        fig.set_tight_layout(True)
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        pf.AxisFormat()

        pf.TufteAxis(ax1, ['left', ], Nticks=[5, 5])
        pf.TufteAxis(ax2, ['left', 'bottom'], Nticks=[5, 5])

        ax1.plot(FirstStage['lambdas'], 
                 np.zeros((len(FirstStage['lambdas']))), 'k', linewidth=1.0)
        ax2.plot(FirstStage['lambdas'], 
                 np.zeros((len(FirstStage['lambdas']))), 'k', linewidth=1.0)
        
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
        
        if savefigs:
            plt.savefig('familyLMS_' + str(int(fracLvM * 100)) + 'L.eps')
        plt.show()

    if plotModel:

        model = cm.colorModel(age=age)
        ax = pf.get_axes(1, 1, nticks=[5, 4])[0]
        ax.spines['left'].set_smart_bounds(False)
        
        style = ['-', '--', '-.']
        for i, LvM in enumerate(np.array([0.0, 0.2, 0.4]) + fracLvM):
            model.genModel(ConeRatio={'fracLvM': LvM, 's': 0.05, },
                           maxSens=maxSens, OD=OD)
            FirstStage = model.returnFirstStage() 
            ThirdStage = model.returnThirdStage()  
        
            # get colors
            blue = ThirdStage['lCenter'].clip(0, 1000)
            yellow = ThirdStage['lCenter'].clip(-1000, 0)
            red = ThirdStage['mCenter'].clip(0, 1000)
            green = ThirdStage['mCenter'].clip(-1000, 0)

            # plot
            ax.plot(FirstStage['lambdas'], blue,
                    'b' + style[i], label=str(int(LvM * 100)) + "%L")
            ax.plot(FirstStage['lambdas'], yellow,
                    'y' + style[i])
            ax.plot(FirstStage['lambdas'], green, 'g' + style[i])
            ax.plot(FirstStage['lambdas'], red, 'r' + style[i])

            # add black line at zero
            ax.plot(FirstStage['lambdas'],
                    np.zeros((len(FirstStage['lambdas']))), 'k', linewidth=2.0)

        ax.set_ylim([-0.28, 0.28])
        ax.set_xlim([FirstStage['wavelen']['startWave'],
                         FirstStage['wavelen']['endWave']])

        ax.legend(loc='upper right', fontsize=18)

        ax.set_ylabel('sensitivity')
        ax.set_xlabel('wavelength (nm)')
        
        if savefigs:
            plt.savefig('percent_L.eps')
            
        plt.show()      
    
    if plotUniqueHues:
        model = cm.colorModel(age=age)
        ax = pf.get_axes(1, 1, nticks=[4, 5])[0]
        ax.spines['bottom'].set_smart_bounds(False)

        style = ['-', '--', '-.']
        i = 0
        for lPeak in [559.0, 557.25, 555.5]:

            model.genModel(maxSens={'l': lPeak, 'm': 530.0, 's': 417.0, }, 
                OD=OD)
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
        ax.set_xlabel('% L')

        if savefigs:
            plt.savefig('unique_hues.eps')
        if SHOW:
            plt.show()
        else:
            return ax

def main(args):
    '''
    '''
    if args.model_type == 'standard':
        import standard_model as cm
    elif args.model_type == 'neitz':
        import colorModel as cm

    if args.LM < 1 and args.LM > 0:
        LMratio = args.LM
    else:
        raise ValueError('LM ratio must be between 0 and 1')

    if args.binom:
        binomPlot(cm)
    
    if args.eccen:
        eccentricityAnalysis(cm)
    
    if args.volbrecht:
        LMratiosAnalysis(cm, Volbrecht1997=True, savefigs=args.save)

    if args.hues:
        HueScaling(cm)

    if args.kuehni:
        MetaAnalysis()

    plotModel(cm,
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
    parser.add_argument("-v", "--volbrecht", action="store_true", 
                        help="plot Carroll/Volbrecht analysis")
    parser.add_argument("-m", "--model", action="store_true",
                        help="plot the neitz color model")
    parser.add_argument("-c", "--curve", action="store_true",
                        help="plot a family of valence curves")
    parser.add_argument("-u", "--unique", action="store_true",
                        help="plot unique hues")
    parser.add_argument("-q", "--hues", action="store_true",
                        help="plot hue scaling data")

    parser.add_argument("-l", "--model_type", default="neitz",
                        help="set model to assume (neitz or standard)")
    parser.add_argument("-s", "--save", action="store_true",
                        help="save plots - not working right now")
    parser.add_argument("--LM", type=float, default=0.25,
                        help="set L:M ratio, default=0.25")
    parser.add_argument("-k", "--kuehni", action='store_true',
                        help="print Kuehni meta analysis of unique hues")

    args = parser.parse_args()
    main(args)
