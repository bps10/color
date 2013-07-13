import numpy as np

from genLMS import genLMS


def genStockmanAnalysis(spectrum, filters, Lnorm, Mnorm, Snorm, scale=0.34):
    '''
    '''

    stage2 = {}
    stage2['L-M'] = Lnorm - Mnorm
    stage2['M-L'] = Mnorm - Lnorm
    LM = Lnorm + (0.5 * Mnorm)
    stage2['S-LM'] = Snorm - (0.69 * LM)
    stage2['LM-S'] = (0.69 * LM) - Snorm

    stage3 = {}
    stage3['red'] = (2.55 * stage2['L-M']) + stage2['S-LM']
    stage3['green'] = (2.55 * stage2['M-L']) + stage2['LM-S']
    stage3['blue'] = (scale * (2.55 * stage2['M-L']) + 
        (Snorm - (scale * 0.69 * LM)))
    stage3['yellow'] = (scale * (2.55 * stage2['L-M']) +
        ((scale * 0.69 * LM) - Snorm))

    print 'RG sys: ', np.sum(stage3['red'])
    print 'BY sys: ', np.sum(stage3['blue'])

    return stage2, stage3
