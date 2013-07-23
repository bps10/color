#! /usr/bin/env python
# -*- coding: utf-8 *-*
from __future__ import division
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
        names=True, usecols= (1, 2, 3, 4, 5, 6, 7))
    print 'hue', 'mean', 'N'
    for name in dat.dtype.names[1:]:
        total, weight = 0, 0
        for i, study in enumerate(dat[name]):

            if not np.isnan(study):
                total += dat['N'][i]
                weight += dat['N'][i] * study

        print name + ': ', weight / total, total

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
            maxSens={'l': carroll['lPeak'][i], 'm': 530.0, 's': 417.0, },
            ConeRatio={'fracLvM': carroll['L'][i] / 100.0, 's': 0.05, })
        uniqueHues = model.get_current_uniqueHues()

        green.append(uniqueHues['green'])
        yellow.append(uniqueHues['yellow'])
        blue.append(uniqueHues['blue'])

    print 'hue | \t mean | \t stdv'
    print 'green: ', np.mean(green), np.std(green)
    print 'yellow: ', np.mean(yellow), np.std(yellow)
    print 'blue: ', np.mean(blue), np.std(blue)   
                       
    BINS = np.arange(0, 101, 6)
    if Volbrecht1997:
        BINS_G = np.arange(488, 564, 3)
        volb = np.genfromtxt('data/Volbrecht1997.txt', delimiter='\t',
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
        

def plotModel(plotModel=True, plotCurveFamily=False,
              plotUniqueHues=False, savefigs=False, fracLvM=0.25):
    """Plot cone spectral sensitivies and first stage predictions.
    """
    
    if plotCurveFamily:
        model = cm.colorModel()
        model.genModel(ConeRatio={'fracLvM': fracLvM, 's': 0.05, })

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
        model = cm.colorModel()
        model.genModel(ConeRatio={'fracLvM': 0.25, 's': 0.05, })

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
        model = cm.colorModel()

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        pf.AxisFormat(linewidth=3)
        pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[4, 5])

        style = ['-', '--', '-.']
        i = 0
        for lPeak in [559.0, 557.0, 555.0]:

            model.genModel(ConeRatio={'fracLvM': 0.25, 's': 0.05, },
                maxSens={'l': lPeak, 'm': 530.0, 's': 417.0, })
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
        plt.show()


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
        LMratiosAnalysis(Volbrecht1997=True)

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
    parser.add_argument("-s", "--save", action="store_true",
                        help="save plots - not working right now")
    parser.add_argument("--LM", type=float, default=0.25,
                        help="set L:M ratio, default=0.25")
    parser.add_argument("-k", "--kuehni", action='store_true',
                        help="print Kuehni meta analysis of unique hues")
    # add default save location

    args = parser.parse_args()
    main(args)