#! /usr/bin/env python
import matplotlib.pylab as plt
import numpy as np

from base.optics import filters
from base import spectsens as ss
from base import plot as pf
from genLMS import genLMS


maxLambda = 770           
filters, spectrum = filters.stockman(minLambda=380, 
            maxLambda=maxLambda, RETURN_SPECTRUM=True, 
            resolution=1)

def plotCompare(compare=['stockman', 'stockSpecSens', 'neitz']):
    '''
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    style = ['-', '--', '-.']
    for i, condition in enumerate(compare):
        
        Lnorm, Mnorm, Snorm = genLMS(spectrum, filters, 
        fundamental=condition, LMSpeaks=[559, 530, 419])
        
        ax.plot(spectrum, Lnorm, 'r' + style[i], linewidth=2)
        ax.plot(spectrum, Mnorm, 'g' + style[i], linewidth=2)
        ax.plot(spectrum, Snorm, 'b' + style[i], linewidth=2)
    #ax.set_ylim([-0.01, 1.01])
    ax.set_xlim([380, 781])
    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('sensitivity')
    plt.tight_layout()
    plt.show()

def plotFilters():
    '''
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    ax.semilogy(spectrum, filters, 'k', linewidth=2)
    ax.set_ylabel('log density')
    ax.set_xlabel('wavelength (nm)')
    ax.set_xlim([380, 781])
    ax.set_ylim([-10, max(filters)])
    plt.tight_layout()
    plt.show()

def plotSpecSens(plot_norm=False):
    '''
    '''
    Lnorm, Mnorm, Snorm = genLMS(spectrum, filters, 
        fundamental='neitz', LMSpeaks=[559, 530, 419])
    L, M, S = genLMS(spectrum, filters, remove_filters=False,
        fundamental='neitz', LMSpeaks=[559, 530, 419])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])

    ax.plot(spectrum, L, 'r-')
    ax.plot(spectrum, M, 'g-')
    ax.plot(spectrum, S, 'b-')
    
    if plot_norm:
        ax.plot(spectrum, Snorm, 'b', linewidth=2)
        ax.plot(spectrum, Mnorm, 'g', linewidth=2)
        ax.plot(spectrum, Lnorm, 'r', linewidth=2)

    ax.set_ylim([-0.01, 1.01])
    ax.set_xlim([380, 781])
    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('sensitivity')
    plt.tight_layout()
    plt.show()


def main(args):
    '''
    '''
    if args.Compare:
        plotCompare()

    if args.Filters:
        plotFilters()

    if args.SpecSens:
        plotSpecSens()

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description="Color Space: display Neitz or Stockman\
        derived color spaces")

    parser.add_argument("-c", "--Compare", action="store_true", 
                        help="compare stockman and neitz fundamentals")
    parser.add_argument("-f", "--Filters", action="store_true",
                        help="plot lens and macula filters")
    parser.add_argument("-s", "--SpecSens", action="store_true", 
                         help="display spectral sensitivities")
    args = parser.parse_args()
    main(args)