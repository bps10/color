#! /usr/bin/env python
import matplotlib.pylab as plt

import numpy as np

from base import plot as pf
from colorSpace import colorSpace


def plotConeSpace():
    '''
    '''
    space = colorSpace(fundamental='neitz',
                             LMSpeaks=[559.0, 530.0, 421.0])
    space._plotColorSpace(space.Lnorm, space.Mnorm, space.spectrum)
    plt.show()

def plotLUV():
    '''
    '''
    # make sure stim is cie 1931 and fundamentals neitz
    space = colorSpace(fundamental='neitz',
                             LMSpeaks=[559.0, 530.0, 421.0],
                             stim='cie 1931')
    
    space.genXYZ()
    u = 4 * space.X / (-2 * space.X + 12 * space.Y + 3)
    v = 9 * space.Y / (-2 * space.X + 12 * space.Y + 3)
    
    ind1 = np.where(space.spectrum == 420)[0]
    ind2 = np.where(space.spectrum == 700)[0]
    spectrum = space.spectrum[ind1:ind2+1]
    u = u[ind1:ind2]
    v = v[ind1:ind2]
    
    space._plotColorSpace(u, v, spectrum, ee=False, Luv=True, 
                         skipLam=[530, 550])
    space.cs_ax.axes.get_xaxis().set_visible(False)
    space.cs_ax.axes.get_yaxis().set_visible(False)
    plt.axis('off')
    plt.show()
    

def plotCIE():
    '''
    '''        
    space = colorSpace(fundamental='neitz',
                             LMSpeaks=[559.0, 530.0, 417.0])          
    sub1_Neitz, sub2_Neitz, jv = space.genKaiser()
    ## plot confusion lines
    clip_area = Wedge((jv[0, 0], jv[0, 1]), r=10, theta1=0, theta2=360)
    CIEcopunctuals = {'deutan': np.array([1.10, -0.1, 0.1]),
                      'protan': np.array([0.753, 0.247, 0]), 
                      'tritan': np.array([0.17, 0, 0.83]),
                      }
    for deficit in CIEcopunctuals:
        space._plotColorSpace(rVal=jv[:, 0], gVal=jv[:, 1],
                                 spec=space.spectrum)     
        
        print deficit, ': ', CIEcopunctuals[deficit]
        
        if deficit.lower() == 'deutan' or deficit.lower() == 'protan':
            lambdas = [420, 460, 470, 480, 490, 500, 515,]
        elif deficit.lower() == 'tritan':
            lambdas = [420, 460, 480, 500, 520, 535, 545, 555,
                       570, 585, 600, 625, 700]
        
        space.cs_ax.plot(CIEcopunctuals[deficit][0],
                        CIEcopunctuals[deficit][1], 'ko', markersize=8)
        for lam in lambdas:
            R, G, B = jv[:, 0], jv[:, 1], jv[:, 2]
            space.cs_ax.plot([space.find_testLightMatch(lam, 
                                R, G, B)[0],
                             CIEcopunctuals[deficit][0]],

                            [space.find_testLightMatch(lam, 
                                R, G, B)[1],
                             CIEcopunctuals[deficit][1]],
                            'k-', linewidth=1) 
            space.cs_ax.plot([jv[-1, 0], jv[0, 0]], 
                    [jv[-1, 1], jv[0, 1]], 'k-', linewidth=3)
                    
            space.cs_ax.set_clip_path(clip_area)
            
        space.cs_ax.set_ylim([-0.12, 0.9])
        space.cs_ax.set_xlim([-0.05, 1.15])          
        space.cs_ax.set_xlabel('x', fontsize=10)
        space.cs_ax.set_ylabel('y', fontsize=10)
        space.cs_ax.text(0.8, 1, deficit, fontsize=18,
                        horizontalalignment='right',
                        verticalalignment='top',
                        transform=space.cs_ax.transAxes)
        fig.set_tight_layout(True)
        plt.show()        

def plotCMFs():
    '''
    '''
    space = colorSpace()

    fig = plt.figure()
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    ax.plot(space.spectrum, space.CMFs[0, :], 'r', linewidth=2)
    ax.plot(space.spectrum, space.CMFs[1, :], 'g', linewidth=2)
    ax.plot(space.spectrum, space.CMFs[2, :], 'b', linewidth=2)
    ax.set_xlim([space.spectrum[0], space.spectrum[-1]])
    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('sensitivity')
    plt.show()

def plotcoeff():
    '''
    '''
    space = colorSpace()

    fig = plt.figure()
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    ax.plot(space.spectrum, space.rVal, 'r', linewidth=2)
    ax.plot(space.spectrum, space.gVal, 'g', linewidth=2)
    ax.plot(space.spectrum, space.bVal, 'b', linewidth=2)
    ax.set_xlim([space.spectrum[0], space.spectrum[-1]])
    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('coefficients')
    plt.show()


def plotColorSpace(color):
    '''
    '''
    space = colorSpace(fundamental='neitz', LMSpeaks=[559, 530, 421])
    space._plotColorSpace(color=color)
    plt.show()


def plot_dichromatic_system(space, hybrid='ls', clip=True):
    '''
    '''
    space = colorSpace(fundamental='neitz', LMSpeaks=[559, 530, 421])
    space._plotColorSpace()
    
    for x in np.arange(0, 1.1, 0.1):
        if hybrid.lower() == 'ls' or hybrid.lower() == 'sl':
            s = x
            m = 0
            l = -(1.0 - x)
        elif hybrid.lower() == 'lm' or hybrid.lower() == 'ml':
            s = 0
            m = x
            l = -(1.0 - x)
        elif hybrid.lower() == 'ms' or hybrid.lower() == 'sm':
            s = x
            m = -(1.0 - x)
            l = 0
        else:
            raise InputError('hybrid must be ls, lm or ms')

        neutral_points, RG = space.find_spect_neutral([l, m, s], True)
        for neut in neutral_points:
            space.cs_ax.plot([neut[0], RG[0]], [neut[1], RG[1]], 
                             '-o', c=(np.abs(l), np.abs(m), np.abs(s)), 
                             markersize=8, linewidth=2)

    if clip is True:                
        space.cs_ax.set_xlim([-0.4, 1.2])
        space.cs_ax.set_ylim([-0.2, 1.2])
    
    plt.show()
        

def plotConfusionLines(deficit='tritan', clip=True):
    '''add confusion lines
        '''
    space = colorSpace(fundamental='neitz', LMSpeaks=[559, 530, 417])
    space._plotColorSpace()
    space.find_copunctuals()
    print deficit, ': ', space.copunctuals[deficit]
    
    if deficit.lower() == 'deutan' or deficit.lower() == 'protan':
        lambdas = [420, 460, 470, 480, 490, 500, 515,]
    elif deficit.lower() == 'tritan':
        lambdas = [420, 460, 480, 500, 520, 535, 545, 555,
                   570, 585, 600, 625, 700]
    
    space.cs_ax.plot(space.copunctuals[deficit][0],
                    space.copunctuals[deficit][1], 'ko', markersize=8)
    for lam in lambdas:
        space.cs_ax.plot([space.find_testLightMatch(lam)[0],
                         space.copunctuals[deficit][0]],
                        [space.find_testLightMatch(lam)[1],
                         space.copunctuals[deficit][1]],
                        'k-', linewidth=1)   
    
    space.cs_ax.text(0.7, 1, deficit, fontsize=18)
    if clip is True:                
        space.cs_ax.set_xlim([-0.4, 1.2])
        space.cs_ax.set_ylim([-0.2, 1.2])
    plt.show()                 


def main(args):
    '''
    '''
    if args.CMFs:
        plotCMFs()
    
    if args.coeff:
        plotcoeff()
    
    if args.ColorSpace:
        plotColorSpace(color=args.color)
    
    if args.ConfusionLines:
        plotConfusionLines()

    if args.dichromat:
        plot_dichromatic_system(False, hybrid=args.hybrid)
    
    if args.ConeSpace:
        plotConeSpace()
    
    if args.LUV:
        plotLUV()


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description="Color Space: display Neitz or Stockman\
        derived color spaces")
    

    parser.add_argument("-y", "--CMFs", action="store_true",
                        help="plot color matching functions")    
    parser.add_argument("-f", "--coeff", action="store_true",
                        help="plot x,y,z coefficients")
    parser.add_argument("-i", "--ColorSpace", action="store_true",
                        help="plot color space")
    parser.add_argument("-c", "--ConfusionLines", action="store_true",
                        help="plot color space with confusion lines")
    parser.add_argument("-q", "--dichromat", action="store_true",
                        help="plot blue-yellow system on color space")   
    parser.add_argument("--hybrid", type=str, default='ls',
                        help="set hybrid pigment for dichromatic system")
    parser.add_argument("-o", "--ConeSpace", action="store_true",
                        help="displace cone space plot")
    parser.add_argument("-l", "--LUV", action="store_true",
                        help="display best fit LUV space") 
    parser.add_argument("--color", action="store_true",
                        help="add color to Neitz color space")
    
    args = parser.parse_args()
    main(args)
