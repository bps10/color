#! /usr/bin/env python
from __future__ import division
import matplotlib.pylab as plt
import numpy as np

from base import plot as pf
from colorSpace import colorSpace


def plot_abney(hybrid='ms', ratio=0.5):
    '''
    '''
    fig = plt.figure()
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    ax.tick_params(labelsize=24, direction='out',
                   size=14)
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])

    space = colorSpace(fundamental='neitz', LMSpeaks=[559, 530, 421])
    space._plotColorSpace()    
    # hybrid pigment weights
    if hybrid.lower() == 'ls' or hybrid.lower() == 'sl':
        if hybrid.lower()[0] == 'l':
            l = ratio
            s = -(1.0 - ratio)
        else:
            s = ratio
            l = -(1.0 - ratio)
        m = 0
    elif hybrid.lower() == 'lm' or hybrid.lower() == 'ml':
        if hybrid.lower()[0] == 'l':
            l = ratio
            m = -(1.0 - ratio)
        else:
            m = ratio
            l = -(1.0 - ratio)
        s = 0
    elif hybrid.lower() == 'ms' or hybrid.lower() == 'sm':
        if hybrid.lower()[0] == 'm':
            m = ratio
            s = -(1.0 - ratio)
        else:
            s = ratio
            m = -(1.0 - ratio)
        l = 0
    else:
        raise InputError('hybrid must be ls, lm or ms')

    neutral_points, copunct = space.find_spect_neutral([l, m, s], True)
    rs = np.array([neutral_points[0][0], neutral_points[1][0]])
    gs = np.array([neutral_points[0][1], neutral_points[1][1]])
    if not np.all(np.diff(rs) > 0):
        rs = rs[::-1]
        gs = gs[::-1]

    r_vals = np.linspace(rs.min(), rs.max(), 20)
    lms = np.zeros((len(r_vals), 3))
    #ratio = np.zeros((len(r_vals), 1))


    for i, r in enumerate(r_vals):
        g = np.interp(r, rs, gs)
        rgb = np.array([r, g, 1 - (r + g)])
        lms[i, :] = space.rgb_to_lms(rgb)

        space.cs_ax.plot(r, g, 'ko')

        
    k1 = m * lms[0, 1] / (m * lms[0, 1] + -s * lms[0, 2])
    k2 = np.zeros((len(r_vals), 1))
    for i in range(len(r_vals)):
        l = lms[i, 0]
        m = k1 * lms[i, 1]
        s = (1 - k1) * lms[i, 1]
        rgb = space.lms_to_rgb([l, m, s])
        print lms[i, 1] / (lms[i, 1] + lms[i, 2])

    #print ratio
    ax.plot(r_vals, lms[:, 0], 'ro-')
    ax.plot(r_vals, lms[:, 1], 'go-')
    ax.plot(r_vals, lms[:, 2], 'bo-')
    
    ax.set_ylabel('relative excitation', fontsize=24)
    ax.set_xlabel('r coordinate', fontsize=24)
    
    plt.show()

def find_coords(hybrid='ms'):
    # Use purity to search for m + s system with same ratio

    '''wvs = np.linspace(490, 560, len(r_vals))
    lms_p = np.zeros((len(wvs), 3))
    for i, wv in enumerate(wvs):
        r, g, b  = space.find_testLightMatch(testLight=wv)
        dist = np.linalg.norm([r, g] - EEW)
        g_d = ((g - (1 / 3)) * purity) + (1 / 3)
        r_d = ((r - (1 / 3)) * purity) + (1 / 3)
        
        space.cs_ax.plot(r_d, g_d, 'gs')
        lms_p[i, :] = space.rgb_to_lms([r_d, g_d, 1 - (r_d + g_d)])

    # interpolate
    # need to correct for diffs
    ratio = lms_p[:, 1] / (lms_p[:, 1] + lms_p[:, 2])
    wv = np.interp(0.96, ratio, wvs)
    print ratio
    print wv
    ax.plot(wvs, lms_p[:, 0], 'ro-')
    ax.plot(wvs, lms_p[:, 1], 'go-')
    ax.plot(wvs, lms_p[:, 2], 'bo-')'''

    space = colorSpace(fundamental='neitz', LMSpeaks=[559, 530, 421])
    space._plotColorSpace()    

    purities = np.arange(0.1, 1.05, 0.05)
    EEW = np.array([1 / 3, 1 / 3])
    more_pigs = np.arange(0, 0.8, 0.0025)
    for purity in purities:
        pigments = np.arange(0, 1.01, 0.05)
        ratios = np.zeros(len(pigments))
        lms_p = np.zeros((len(pigments), 3))
        for i, p in enumerate(pigments):
            s = p
            m = -(1.0 - p)
            l = 0

            neutral_points, copunct = space.find_spect_neutral([l, m, s], True)
            if purity == purities[0]:
                for neut in neutral_points:
                    space.cs_ax.plot([neut[0], copunct[0]], 
                                     [neut[1], copunct[1]], 
                                     '-o', 
                                     c=(np.abs(l), np.abs(m), np.abs(s)), 
                                     markersize=8, linewidth=2)

            r = neutral_points[0][0]
            g = neutral_points[0][1]
            dist = np.linalg.norm([r, g] - EEW)
            g_d = ((g - (1 / 3)) * purity) + (1 / 3)
            r_d = ((r - (1 / 3)) * purity) + (1 / 3)
            #space.cs_ax.plot(r_d, g_d, 'ks')

            lms_p[i, :] = space.rgb_to_lms([r_d, g_d, 1 - (r_d + g_d)])
            #ratios[i] = lms_p[i, 1] / (lms_p[i, 1] + lms_p[i, 2])
            ratios[i] = lms_p[i, 1] / lms_p[i, 2] * (1 / 32.3)
        if purity > 0.96:
            print purity
            print ratios
        ratios = np.interp(more_pigs, pigments, ratios - 0.967)
        zero_crossing = np.where(np.diff(np.sign(ratios)))[0]
        for cross in zero_crossing:
            pig = more_pigs[cross]
            m_ = -(1 - pig)
            s_ = pig
            print purity
            print m_, s_
            neutral_points, copunct = space.find_spect_neutral([0, m_, s_], 
                                                               True)

            r = neutral_points[0][0]
            g = neutral_points[0][1]
            dist = np.linalg.norm([r, g] - EEW)
            g_d = ((g - (1 / 3)) * purity) + (1 / 3)
            r_d = ((r - (1 / 3)) * purity) + (1 / 3)
            space.cs_ax.plot(r_d, g_d, 'ko')    

    space.cs_ax.set_xlim([-0.4, 1.2])
    space.cs_ax.set_ylim([-0.2, 1.2])
    
    plt.show()




def abney():

    fig = plt.figure()
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    ax.tick_params(labelsize=24, direction='out',
                   size=14)
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])

    space = colorSpace(fundamental='neitz', LMSpeaks=[559, 530, 421])
    space._plotColorSpace()    

    purity = 0.8
    EEW = np.array([1 / 3, 1 / 3])

    pigments = np.arange(0, 1.01, 0.05)
    lms_p = np.zeros((len(pigments), 3))
    lms_w = np.zeros((len(pigments), 3))

    for i, p in enumerate(pigments):
        s = p
        m = -(1.0 - p)
        l = 0

        neutral_points, copunct = space.find_spect_neutral([l, m, s], True)

        for neut in neutral_points:
            space.cs_ax.plot([neut[0], copunct[0]], 
                             [neut[1], copunct[1]], 
                             '-o', 
                             c=(np.abs(l), np.abs(m), np.abs(s)), 
                             markersize=8, linewidth=2)
                
        r = neutral_points[0][0]
        g = neutral_points[0][1]
        dist = np.linalg.norm([r, g] - EEW)
        g_d = ((g - (1 / 3)) * purity) + (1 / 3)
        r_d = ((r - (1 / 3)) * purity) + (1 / 3)
        space.cs_ax.plot(r_d, g_d, 'ks')

        rgb_w = np.array([r, g, 1 - (r + g)])
        rgb_p = np.array([r_d, g_d, 1 - (r_d + g_d)])
        lms_p[i, :] = space.rgb_to_lms(rgb_p)
        lms_w[i, :] = space.rgb_to_lms(rgb_w)
    
    ratio_w = np.log10(lms_w[11, 1] / lms_w[11, 2])
    ratio_p = np.log10(lms_p[:, 1] / lms_p[:, 2])
    new_pigment = np.interp(ratio_w, ratio_p, pigments)
    print ratio_w
    print ratio_p
    print pigments
    print new_pigment


def plot_lms():
    fig = plt.figure()
    fig.set_tight_layout(True)
    ax = fig.add_subplot(111)
    ax.tick_params(labelsize=24, direction='out',
                   size=14)
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])

    space = colorSpace(fundamental='neitz', LMSpeaks=[559, 530, 421])
    space._plotColorSpace()    

    purity = 0.8
    EEW = np.array([1 / 3, 1 / 3])

    pigments = np.arange(0, 1.01, 0.05)
    lms_p = np.zeros((len(pigments), 3))
    lms_w = np.zeros((len(pigments), 3))

    for i, p in enumerate(pigments):
        s = p
        m = -(1.0 - p)
        l = 0

        neutral_points, copunct = space.find_spect_neutral([l, m, s], True)

        for neut in neutral_points:
            space.cs_ax.plot([neut[0], copunct[0]], 
                             [neut[1], copunct[1]], 
                             '-o', 
                             c=(np.abs(l), np.abs(m), np.abs(s)), 
                             markersize=8, linewidth=2)
                
        r = neutral_points[0][0]
        g = neutral_points[0][1]
        dist = np.linalg.norm([r, g] - EEW)
        g_d = ((g - (1 / 3)) * purity) + (1 / 3)
        r_d = ((r - (1 / 3)) * purity) + (1 / 3)
        space.cs_ax.plot(r_d, g_d, 'ks')

        rgb_w = np.array([r, g, 1 - (r + g)])
        rgb_p = np.array([r_d, g_d, 1 - (r_d + g_d)])
        lms_p[i, :] = space.rgb_to_lms(rgb_p)
        lms_w[i, :] = space.rgb_to_lms(rgb_w)

    space.cs_ax.set_xlim([-0.4, 1.2])
    space.cs_ax.set_ylim([-0.2, 1.2])

    ax.plot(pigments, lms_p[:, 0], 'ro-')
    ax.plot(pigments, lms_p[:, 1], 'go-')
    ax.plot(pigments, lms_p[:, 2], 'bo-')

    ax.plot(pigments, lms_w[:, 0], 'ro--')
    ax.plot(pigments, lms_w[:, 1], 'go--')
    ax.plot(pigments, lms_w[:, 2], 'bo--')
    
    ax.set_ylabel('relative activation', fontsize=24)
    ax.set_xlabel('S:M pigment ratio', fontsize=24)

    fig = plt.figure()
    fig.set_tight_layout(True)
    ax1 = fig.add_subplot(111)
    ax1.tick_params(labelsize=24, direction='out',
                   size=14)
    pf.TufteAxis(ax1, ['left', 'bottom'], Nticks=[5, 5])
    ax1.semilogy(pigments, lms_p[:, 1] / lms_p[:, 2], 'ko-')
    ax1.semilogy(pigments, lms_w[:, 1] / lms_w[:, 2], 'ko--')

    ax1.set_ylabel('M/S activation', fontsize=24)
    ax1.set_xlabel('S:M pigment ratio', fontsize=24)
    
    plt.show()

 
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description="Color Space: derived color spaces")
    
    parser.add_argument("-d", "--hybrid", type=str, default='ls',
                        help="set hybrid pigment for dichromatic system")
    parser.add_argument("-r", "--ratio", type=float, default=0.5,
                        help="set hybrid pigment weights")
    
    args = parser.parse_args()

    #find_coords()
    plot_lms()
    #abney()
    #plot_abney(args.hybrid, args.ratio)
