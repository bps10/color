# -*- coding: utf-8 *-*
from __future__ import division
import numpy as np
import matplotlib.pylab as plt
from math import factorial

from spectsens import spectsens
import PlottingFun as pf
from colorModel import colorModel, getCarroll_LMratios

def LMratiosAnalysis(Volbrecht1997=True, returnVals=False, 
                        plot=True, savefigs=False):
    '''Creates a dictionary like object.
    '''

    model = colorModel(q=1.300)
    
    model.genModel()
    model.findUniqueHues()
    uniqueHues = model.returnUniqueHues()
    carroll = getCarroll_LMratios()

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
                       
    BINS = np.arange(0, 101, 6)
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

    LMratiosAnalysis(Volbrecht1997=True)
    #plotModel(plotSpecSens=False, plotCurveFamily=True,
    #          plotUniqueHues=False, savefigs=False)
