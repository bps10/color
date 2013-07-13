import numpy as np

from base import spectsens as ss
from base.optics import filters


def genLMS(spectrum, filters, fundamental='neitz', LMSpeaks=[559.0, 530.0, 419.0],
			remove_filters=True, normalize=True):
    '''
    '''
    if len(LMSpeaks) != 3:
        print 'LMSpeaks must be length 3! Using defaults: 559, 530, 417nm'
        LMSpeaks = [559.0, 530.0, 421.0]
 
    minLam = np.max(spectrum)
    maxLam = np.min(spectrum) 
    step = spectrum[1] - spectrum[0]

    if fundamental.lower() == 'stockman':
        sens = ss.stockmanfund(minLambda=minLam, maxLambda=maxLam, resolution=step)

        LS = sens[:, 0]
        MS = sens[:, 1]
        SS = sens[:, 2]

        Lresponse = LS * spectrum
        Mresponse = MS * spectrum
        Sresponse = SS * spectrum
        
    elif fundamental.lower() == 'stockspecsens':

        sens = ss.stockman(minLambda=minLam, maxLambda=maxLam, resolution=step)

        LS = sens[:, 0]
        MS = sens[:, 1]
        SS = sens[:, 2]
        
        Lresponse = LS / filters * spectrum
        Mresponse = MS / filters * spectrum
        Sresponse = SS / filters * spectrum
        
    elif fundamental.lower() == 'neitz':

        LS = ss.neitz(LMSpeaks[0], 0.5, False, minLam, 
                                         maxLam, step)
        MS = ss.neitz(LMSpeaks[1], 0.5, False, minLam, 
                                         maxLam, step)
        SS = ss.neitz(LMSpeaks[2], 0.4, False, minLam, 
                                         maxLam, step)
                                        
        Lresponse = LS / filters * spectrum
        Mresponse = MS / filters * spectrum
        Sresponse = SS / filters * spectrum
    
    # normalize to peak at unity
    Lnorm = Lresponse / np.max(Lresponse)
    Mnorm = Mresponse / np.max(Mresponse)
    Snorm = Sresponse / np.max(Sresponse)

    if remove_filters and normalize:
    	return Lnorm, Mnorm, Snorm
    if remove_filters and not normalize:
    	return Lresponse, Mresponse, Sresponse
    if not remove_filters:
    	return LS, MS, SS


def getStockmanFilter(minLambda=390, maxLambda=770, RETURN_SPECTRUM=False):
    '''
    '''
    filters, spectrum = filters.stockman(minLambda=minLambda, 
        maxLambda=maxLambda, RETURN_SPECTRUM=True, 
        resolution=1)
    if RETURN_SPECTRUM:
    	return filters, spectrum
    return filters

