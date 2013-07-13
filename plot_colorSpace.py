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

def plotKaiser(neitz=False, showBY=True, clip=True,
               showSub1=False, showSub2=True, stockman=False,
               series=True):
    '''
    '''
    space = colorSpace(fundamental='neitz',
                             LMSpeaks=[559.0, 530.0, 421.0])        
    sub1_Neitz, sub2_Neitz, jv = space.genKaiser(neitz)

    if neitz:
        space._plotColorSpace()

    else:
        space._plotColorSpace(rVal=jv[:, 0], gVal=jv[:, 1],
                             spec=space.spectrum)
        space.cs_ax.plot([jv[-1, 0], jv[0, 0]], 
                        [jv[-1, 1], jv[0, 1]], 'k-', linewidth=3)
        space.cs_ax.set_ylim([0, 0.9])
        space.cs_ax.set_xlim([-0.05, 0.8])
    if showSub1:     
        space.cs_ax.plot(sub1_Neitz[:, 0], sub1_Neitz[:, 1], 'ko', 
                        markersize=8, markeredgewidth=2,
                        markerfacecolor='w',
                        linewidth=2)
    if showSub2:
        space.cs_ax.plot(sub2_Neitz[:, 0], sub2_Neitz[:, 1], 'kx',
                        markersize=8, markeredgewidth=2, linewidth=2)

    if showBY:
        if stockman:
            neut2, RG2 = space.BY2lambda(0, 0, 1., True)
            c2 = (1, 0, 0)
            c3 = (0.5, 0.5, 0)
            neut3, RG3 = space.lambda2RG(522, False, True)
        else:
            neut2, RG2 = space.BY2lambda(1, 0, 0, True)
            c2 = (0, 0, 1)
            c3 = (0, 0.5, 0.5)
            neut3, RG3 = space.lambda2BY(522, True)
        neut1, RG1 = space.BY2lambda(0, 1., 0, True)

        c1 = (0, 1, 0)
        # plot green copunctual line
        space.cs_ax.plot([neut1[0], RG1[0]], [neut1[1], RG1[1]], 
                        '-o', c=c1, markersize=8, linewidth=2)  
        # plot red or blue copunctual depending on neitz or stockman
        space.cs_ax.plot([neut2[0], RG2[0]], [neut2[1], RG2[1]], 
                        '-o', c=c2, markersize=8, linewidth=2)  
        # plot 
        space.cs_ax.plot([neut3[0], RG3[0]], [neut3[1], RG3[1]], 
                        '-o', c=c3, markersize=8, linewidth=2)  

    if stockman and series:

        for lam in [500, 505, 510, 515]:
            neut3, RG3 = space.lambda2RG(lam, False, True)
            space.cs_ax.plot([neut3[0], RG3[0]], [neut3[1], RG3[1]], 
                '-o', c=c3, markersize=8, linewidth=2)  


    if clip is True:                
        space.cs_ax.set_xlim([-0.4, 1.2])
        space.cs_ax.set_ylim([-0.2, 1.2])
    
    space.cs_ax.set_xlabel('x', fontsize=10)
    space.cs_ax.set_ylabel('y', fontsize=10)
    
    plt.tight_layout()
    plt.show()

def plotCIE():
    '''
    '''        
    space = colorSpace(fundamental='neitz',
                             LMSpeaks=[559.0, 530.0, 421.0])          
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
        plt.tight_layout()
        plt.show()        

def plotCMFs():
    '''
    '''
    space = colorSpace()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    ax.plot(space.spectrum, space.CMFs[0, :], 'r', linewidth=2)
    ax.plot(space.spectrum, space.CMFs[1, :], 'g', linewidth=2)
    ax.plot(space.spectrum, space.CMFs[2, :], 'b', linewidth=2)
    ax.set_xlim([space.spectrum[0], space.spectrum[-1]])
    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('sensitivity')
    plt.tight_layout()
    plt.show()

def plotcoeff():
    '''
    '''
    space = colorSpace()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    ax.plot(space.spectrum, space.rVal, 'r', linewidth=2)
    ax.plot(space.spectrum, space.gVal, 'g', linewidth=2)
    ax.plot(space.spectrum, space.bVal, 'b', linewidth=2)
    ax.set_xlim([space.spectrum[0], space.spectrum[-1]])
    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('coefficients')
    plt.tight_layout()
    plt.show()

def plotColorSpace():
    '''
    '''
    space = colorSpace(fundamental='neitz', LMSpeaks=[559, 530, 419])
    space._plotColorSpace()
    plt.show()

def plotBYsystem(space, PRINT=False, clip=True):
    '''
    '''
    space = colorSpace(fundamental='neitz', LMSpeaks=[559, 530, 419])

    space._plotColorSpace()
    
    for s in range(0, 11):
        m = (10.0 - s) / 10.0
        s = s / 10.0

        neut, RG = space.BY2lambda(s, m, 0, True)
        if PRINT is True:
            #print RG
            #print neut
            print space.find_testlightFromRG(neut[0], neut[1])
        space.cs_ax.plot([neut[0], RG[0]], [neut[1], RG[1]], 
                        '-o', c=(0, m, s), markersize=8, linewidth=2)
    
    if clip is True:                
        space.cs_ax.set_xlim([-0.4, 1.2])
        space.cs_ax.set_ylim([-0.2, 1.2])
    
    plt.show()
        
def plotRGsystem(PRINT=False, clip=True):
    '''
    '''
    space = colorSpace(fundamental='neitz', LMSpeaks=[559, 530, 419])
    space._plotColorSpace()
    
    for l in range(0, 11):
        m = (10.0 - l) / 10.0
        l = l / 10.0
        
        neut, RG = space.RG2lambda(0, m, l, True)
        
        if PRINT is True:
            #print RG
            #print neut
            print space.find_testlightFromRG(neut[0], neut[1])
        space.cs_ax.plot([neut[0], RG[0]], [neut[1], RG[1]], 
                        '-o', c=(l, m, 0), markersize=8, linewidth=2)
    
    if clip is True:                
        space.cs_ax.set_xlim([-0.4, 1.2])
        space.cs_ax.set_ylim([-0.2, 1.2])
    
    plt.show()

def plotConfusionLines(deficit='tritan', clip=True):
    '''add confusion lines
        '''
    space = colorSpace(fundamental='neitz', LMSpeaks=[559, 530, 419])
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

def genStockmanAnalysis(scale=0.34):
    '''
    '''
    space = colorSpace(fundamental='neitz', LMSpeaks=[559, 530, 419])
    space._plotColorSpace()
    space.genLMS('stockman', [559, 530, 421])

    stage2 = {}
    stage2['L-M'] = space.Lnorm - space.Mnorm
    stage2['M-L'] = space.Mnorm - space.Lnorm
    LM = space.Lnorm + (0.5 * space.Mnorm)
    stage2['S-LM'] = space.Snorm - (0.69 * LM)
    stage2['LM-S'] = (0.69 * LM) - space.Snorm

    stage3 = {}
    stage3['red'] = (2.55 * stage2['L-M']) + stage2['S-LM']
    stage3['green'] = (2.55 * stage2['M-L']) + stage2['LM-S']
    stage3['blue'] = (scale * (2.55 * stage2['M-L']) + 
        (space.Snorm - (scale * 0.69 * LM)))
    stage3['yellow'] = (scale * (2.55 * stage2['L-M']) +
        ((scale * 0.69 * LM) - space.Snorm))

    print 'RG sys: ', np.sum(stage3['red'])
    print 'BY sys: ', np.sum(stage3['blue'])

    return stage2, stage3

def main(args):
    '''
    '''

    if args.CMFs:
        plotCMFs()
    
    if args.coeff:
        plotcoeff()
    
    if args.ColorSpace:
        plotColorSpace()
    
    if args.ConfusionLines:
        plotConfusionLines()

    if args.BYsystem:
        plotBYsystem(False)

    if args.RGsystem:
        plotRGsystem(False)
    
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
    parser.add_argument("-q", "--BYsystem", action="store_true",
                        help="plot blue-yellow system on color space")   
    parser.add_argument("-p", "--RGsystem", action="store_true",
                        help="plot red-green system on color space") 
    parser.add_argument("-o", "--ConeSpace", action="store_true",
                        help="displace cone space plot")
    parser.add_argument("-l", "--LUV", action="store_true",
                        help="display best fit LUV space") 
    
    args = parser.parse_args()
    main(args)
