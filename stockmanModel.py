import numpy as np

from base import spectsens as ss

def genLMS():
    '''
    '''
    ind = np.max(self.spectrum)

    sens = ss.stockman(minLambda=390, maxLambda=ind)

    LS = sens[:, 0]
    MS = sens[:, 1]
    SS = sens[:, 2]

    Lresponse = LS / self.filters * self.spectrum
    Mresponse = MS / self.filters * self.spectrum
    Sresponse = SS / self.filters * self.spectrum
    self.Lnorm = Lresponse / np.max(Lresponse)
    self.Mnorm = Mresponse / np.max(Mresponse)
    self.Snorm = Sresponse / np.max(Sresponse)

def genStockmanAnalysis(scale=0.34):
    '''
    '''
    self.genLMS('stockman', [559, 530, 421])

    stage2 = {}
    stage2['L-M'] = self.Lnorm - self.Mnorm
    stage2['M-L'] = self.Mnorm - self.Lnorm
    LM = self.Lnorm + (0.5 * self.Mnorm)
    stage2['S-LM'] = self.Snorm - (0.69 * LM)
    stage2['LM-S'] = (0.69 * LM) - self.Snorm

    stage3 = {}
    stage3['red'] = (2.55 * stage2['L-M']) + stage2['S-LM']
    stage3['green'] = (2.55 * stage2['M-L']) + stage2['LM-S']
    stage3['blue'] = (scale * (2.55 * stage2['M-L']) + 
        (self.Snorm - (scale * 0.69 * LM)))
    stage3['yellow'] = (scale * (2.55 * stage2['L-M']) +
        ((scale * 0.69 * LM) - self.Snorm))

    print 'RG sys: ', np.sum(stage3['red'])
    print 'BY sys: ', np.sum(stage3['blue'])

    return stage2, stage3
