#! /usr/bin/env python
import matplotlib.pylab as plt
import numpy as np

from base import optics as op
from base import plot as pf
from stockmanModel import genStockmanAnalysis
from genLMS import genLMS


maxLambda = 770
filters, spectrum = op.filters.stockman(minLambda=380, 
    maxLambda=maxLambda, RETURN_SPECTRUM=True, 
    resolution=1)
Lnorm, Mnorm, Snorm = genLMS(spectrum, filters, 
    fundamental='stockman', LMSpeaks=[559, 530, 421])

stage2, stage3 = genStockmanAnalysis(spectrum, filters, Lnorm,
    Mnorm, Snorm)

def plotStage2Stockman():
    '''
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    for key in stage2:
        ax.plot(spectrum, stage2[key])

    ax.set_xlim([spectrum[0], spectrum[-1]])
    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('sensitivity')
    plt.tight_layout()
    plt.show()


def plotStage3Stockman():
    '''
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    for key in stage3:
        ax.plot(spectrum, stage3[key], c=key)

    ax.set_xlim([spectrum[0], spectrum[-1]])
    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('sensitivity')
    plt.tight_layout()
    plt.show()

def plotUniqueGreenSeries():
    '''
    '''
    # Unique green series plot
    fig1 = plt.figure()
    fig2 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax1, ['left', 'bottom'], Nticks=[5, 5])
    pf.TufteAxis(ax2, ['left', 'bottom'], Nticks=[5, 5])

    for j in [0.04, 0.34, 0.94, 2, 4, 6]:
        stage2, stage3 = genStockmanAnalysis(spectrum, filters, Lnorm,
            Mnorm, Snorm, j)
        ax1.plot(spectrum, stage3['blue'], c='b', alpha=0.7)
        ax2.plot(spectrum, Snorm - (j / 10 * (Lnorm + (0.5 * Mnorm))), 
            c='b', alpha=0.7)

    ax1.plot(spectrum, np.zeros(len(spectrum)), 'k',
        linewidth=1)
    ax2.plot(spectrum, np.zeros(len(spectrum)), 'k',
        linewidth=1)
    ax1.set_xlim([spectrum[0], 650])
    ax1.set_ylim([-0.7, 1.4])
    ax1.set_xlabel('wavelength (nm)')
    ax1.set_ylabel('sensitivity')

    ax2.set_xlim([spectrum[0], 700])
    ax2.set_ylim([-0.9, 1.2])
    ax2.set_xlabel('wavelength (nm)')
    ax2.set_ylabel('sensitivity')

    fig1.tight_layout()
    fig2.tight_layout()
    plt.show()


if __name__ == '__main__':

    #plotStage2Stockman()
    #plotStage3Stockman()
    plotUniqueGreenSeries()


