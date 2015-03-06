#! /usr/bin/env python
import matplotlib.pylab as plt
import numpy as np

from base.optics import filters as filt
from base import spectsens as ss
from base import plot as pf
from genLMS import genLMS


maxLambda = 770           
filters, spectrum = filt.stockman(minLambda=390, 
            maxLambda=maxLambda, RETURN_SPECTRUM=True, 
            resolution=1)

def plotCompare(compare=['stockman', 'stockSpecSens', 'neitz'],
    invert=False):
    '''
    '''
    fig = plt.figure()
    fig.set_tight_layout(True)
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

    if invert:
        pf.invert(ax, fig, bk_color='k')

    plt.show()

def plotFilters(invert=False, log=False, lens_only=True, 
                macula_only=False, lens_age=[20, 40, 60, 80], 
                spectrum=spectrum, stiles=True):
    '''
    '''
    if macula_only:
        filters, spectrum = filt.stockman(minLambda=400, 
                                             maxLambda=700, 
                                             RETURN_SPECTRUM=True, 
                                             ONLY_MACULA=True,
                                             resolution=1)
    if lens_only and lens_age is None:
        filters, spectrum = filt.stockman(minLambda=400, 
                                             maxLambda=700, 
                                             RETURN_SPECTRUM=True, 
                                             ONLY_LENS=True,
                                             resolution=1)
    if lens_only and lens_age is not None:
        spectrum = np.arange(400, 700)
        filters = np.zeros((len(spectrum), len(lens_age)))
        for i, age in enumerate(lens_age):
            filters[:, i] = filt.lens_age_correction(age, spectrum, 
                                                     stiles=stiles)

    fig = plt.figure()
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    if log:
        ax.semilogy(spectrum, filters, 'k', linewidth=2)
    else:
        ax.plot(spectrum, filters, 'k', linewidth=2)


    ax.set_xlabel('wavelength (nm)')
    ax.set_xlim([400, 700])
    if log:
        ax.set_ylim([-10, max(filters)])
        ax.set_ylabel('log density')
    else:
        ax.set_ylabel('density')

    if invert:
        pf.invert(ax, fig, bk_color='k')
    plt.show()


def plotSpecSens(plot_norm=True, log=False, invert=False):
    '''
    '''
    Lnorm, Mnorm, Snorm = genLMS(spectrum, filters, 
        fundamental='neitz', LMSpeaks=[559, 530, 419])
    L, M, S = genLMS(spectrum, filters, remove_filters=False,
        fundamental='neitz', LMSpeaks=[559, 530, 419])

    fig = plt.figure()
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    
    if plot_norm:
        if log:
            ax.semilogy(spectrum, Snorm, 'b', linewidth=2)
            ax.semilogy(spectrum, Mnorm, 'g', linewidth=2)
            ax.semilogy(spectrum, Lnorm, 'r', linewidth=2)
        else:
            ax.plot(spectrum, Snorm, 'b', linewidth=2)
            ax.plot(spectrum, Mnorm, 'g', linewidth=2)
            ax.plot(spectrum, Lnorm, 'r', linewidth=2)
    else:
        if log:
            ax.semilogy(spectrum, S, 'b', linewidth=2)
            ax.semilogy(spectrum, M, 'g', linewidth=2)
            ax.semilogy(spectrum, L, 'r', linewidth=2)
        else:
            ax.plot(spectrum, S, 'b', linewidth=2)
            ax.plot(spectrum, M, 'g', linewidth=2)
            ax.plot(spectrum, L, 'r', linewidth=2)

    if log:
        ax.set_ylim([10 ** -5, 1])
    ax.set_xlim([380, 781])
    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('sensitivity')

    if invert:
        pf.invert(ax, fig, bk_color='k')

    plt.show()


def plotRelativeSens(invert=False):
    '''
    '''
    Lnorm, Mnorm, Snorm = genLMS(spectrum, filters, 
        fundamental='neitz', LMSpeaks=[559, 530, 419])
    L, M, S = genLMS(spectrum, filters, remove_filters=False,
        fundamental='neitz', LMSpeaks=[559, 530, 419])

    fig = plt.figure()
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])

    ax.semilogy(spectrum, L / M, 'k')
    ax.semilogy(spectrum, np.ones(spectrum.size), 'k--')
    #ax.set_ylim([-0.01, 1.01])
    ax.set_xlim([380, 781])
    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('L/M sensitivity ratio')
    if invert:
        pf.invert(ax, fig, bk_color='k')
    plt.show()

def main(args):
    '''
    '''
    if args.Compare:
        plotCompare(invert=args.invert)

    if args.Filters:
        plotFilters(invert=args.invert)

    if args.SpecSens:
        plotSpecSens(invert=args.invert, log=args.log, 
                     plot_norm=args.norm)

    if args.Relative:
        plotRelativeSens(invert=args.invert)

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
    parser.add_argument("-r", "--Relative", action="store_true", 
                        help="display relative LM spectral sensitivities")
    parser.add_argument("-i", "--invert", action='store_true',
                        help="invert background (black)")
    parser.add_argument("-l", "--log", action='store_true',
                        help="make y axis a log plot")
    parser.add_argument("-n", "--norm", action='store_true',
                        help="plot normed fundamentals")

    args = parser.parse_args()

    main(args)
