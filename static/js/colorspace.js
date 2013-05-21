var colorSpace = (typeof exports === "undefined")?(function colorSpace() {}):(exports);
if(typeof global !== "undefined") { global.colorSpace = colorSpace; }

colorSpace.version = "1.2.6";

    // 1. run
	colorSpace.params = {'lights': stim.lower, }
	colorSpace.setLights(stim)
	colorSpace.genStockmanFilter()
	colorSpace.genLMS(fundamental, LMSpeaks)
	colorSpace.genConvMatrix()
	
	colorSpace.LMStoCMFs()
	colorSpace.CMFtoEE_CMF()
	colorSpace.EE_CMFtoRGB()
        
    colorSpace.genLMS = function genLMS(fundamental, LMSpeaks) {

        if (LMSpeaks.length != 3) {
            print 'LMSpeaks must be length 3! Using defaults: 559, 530, 417nm'
            LMSpeaks = [559.0, 530.0, 421.0]
        }
        if (fundamental.toLowerCase() == 'stockman') {
            ind = this.spectrum.length
            foo = np.genfromtxt(STATIC_ROOT + 
                                    '/stockman/fundamentals2deg.csv', 
                                 delimiter=',')[::10, :]
            this.Lc = 10.0 ** foo[:ind, 1]
            this.Mc = 10.0 ** foo[:ind, 2]
            this.Sc = 10.0 ** foo[:ind, 3]
    
            Lresponse = this.Lc * this.spectrum
            Mresponse = this.Mc * this.spectrum
            Sresponse = this.Sc * this.spectrum
        }
        else if (fundamental.lower() == 'stockspecsens') {
            ind = len(this.spectrum)
            foo = np.genfromtxt(STATIC_ROOT + 
                                    '/stockman/specSens.csv', 
                                delimiter=',')[::10, :]

            LS = np.log10((1.0 - 10.0 ** -((10.0 ** foo[:, 1]) *
                    0.5)) / (1.0 - 10 ** -0.5))
            MS = np.log10((1.0 - 10.0 ** -((10.0 ** foo[:, 2]) *
                    0.5)) / (1.0 - 10 ** -0.5))
            SS = np.log10((1.0 - 10.0 ** -((10.0 ** foo[:, 3]) *
                    0.4)) / (1.0 - 10 ** -0.4))
          
            this.Lc = 10.0 ** LS[:ind]
            this.Mc = 10.0 ** MS[:ind]
            this.Sc = 10.0 ** SS[:ind]
            
            Lresponse = this.Lc / this.filters * this.spectrum
            Mresponse = this.Mc / this.filters * this.spectrum
            Sresponse = this.Sc / this.filters * this.spectrum
        }    
        else if (fundamental.lower() == 'neitz') {
            minspec = min(this.spectrum)
            maxspec = max(this.spectrum)
            this.Lc = spectsens(LMSpeaks[0], 0.5, 'anti-log', minspec, 
                                             maxspec, 1)[0]
            this.Mc = spectsens(LMSpeaks[1], 0.5, 'anti-log', minspec, 
                                             maxspec, 1)[0]
            this.Sc = spectsens(LMSpeaks[2], 0.4, 'anti-log', minspec, 
                                             maxspec, 1)[0]
                                                         
            Lresponse = this.Lc / this.filters * this.spectrum
            Mresponse = this.Mc / this.filters * this.spectrum
            Sresponse = this.Sc / this.filters * this.spectrum
        }
        #record param
        this.params['fundamentals'] = fundamental
        this.params['LMSpeaks'] = LMSpeaks
        
        this.Lnorm = Lresponse / np.max(Lresponse)
        this.Mnorm = Mresponse / np.max(Mresponse)
        this.Snorm = Sresponse / np.max(Sresponse)
	}
	
    colorSpace.genStockmanFilter =  function genStockmanFilter(maxLambda=770) {

        lens = np.genfromtxt(STATIC_ROOT + '/stockman/lens.csv', 
                             delimiter=',')[::10, :]
        macula = np.genfromtxt(STATIC_ROOT + 
                                '/stockman/macular.csv', 
                                delimiter=',')[::10, :]

        spectrum = lens[:, 0]
        ind = np.where(spectrum == maxLambda)[0]
        this.spectrum = spectrum[:ind+1]
        
        this.filters = 10.0 ** (lens[:ind + 1, 1] +  macula[:ind + 1, 1])
    }  
	
    colorSpace.genConvMatrix = function genConvMatrix(PRINT=False) {

        this.convMatrix = np.array([
            [np.interp(this.lights['l'], this.spectrum, this.Lnorm),
            np.interp(this.lights['m'], this.spectrum, this.Lnorm),
            np.interp(this.lights['s'], this.spectrum, this.Lnorm)],

            [np.interp(this.lights['l'], this.spectrum, this.Mnorm),
            np.interp(this.lights['m'], this.spectrum, this.Mnorm),
            np.interp(this.lights['s'], this.spectrum, this.Mnorm)],

            [np.interp(this.lights['l'], this.spectrum, this.Snorm),
            np.interp(this.lights['m'], this.spectrum, this.Snorm),
            np.interp(this.lights['s'], this.spectrum, this.Snorm)]])

        if PRINT == True:
            print this.convMatrix
	}
	
    colorSpace.setLights = function setLights(stim) {

        if (stim.toLowerCase() != 'wright' && stim.toLowerCase() != 'stiles and burch' 
            && stim.toLowerCase() != 'cie 1931') {
            throw TypeError('Sorry, stim light not understood, using wright')
            stim = 'wright' 
        }
        if (stim.toLowerCase() == 'wright') {
            this.lights = {
                            'l': 650.0,
                            'm': 530.0,
                            's': 460.0,
                            }
		}
        if (stim.toLowerCase() == 'stiles and burch') {
            this.lights = {'l': 645.0, 
                           'm': 526.0, 
                           's': 444.0, }
		}
        if (stim.toLowerCase() == 'cie 1931') {
            this.lights = {'l': 700.0, 
                           'm': 546.1, 
                           's': 435.8, }
		}
	}
	
    colorSpace.TrichromaticEquation = function TrichromaticEquation(r, g, b) {

        rgb = r + g + b
        r_ = r / rgb
        g_ = g / rgb
        b_ = b / rgb
        
        return r_, g_, b_
    }
	
    colorSpace.LMStoCMFs = function LMStoCMFs() {

        
        LMSsens = np.array([this.Lnorm, this.Mnorm, this.Snorm])
        this.CMFs = np.dot(np.linalg.inv(this.convMatrix), LMSsens)

        #save sums for later normalization:            
        Rnorm = sum(this.CMFs[0, :])
        Gnorm = sum(this.CMFs[1, :])
        Bnorm = sum(this.CMFs[2, :])
        this.EEfactors = {'r': Rnorm, 'g': Gnorm, 'b': Bnorm, }
    }
	
    colorSpace.CMFtoEE_CMF = function CMFtoEE_CMF() {

        this.CMFs[0, :], this.CMFs[1, :], this.CMFs[2, :] = this._EEcmf(
                                        this.CMFs[0, :], 
                                        this.CMFs[1, :], 
                                        this.CMFs[2, :])

	}
	
    colorSpace.EE_CMFtoRGB = function EE_CMFtoRGB() {

        this.rVal, this.gVal, this.bVal = this.TrichromaticEquation(
                            this.CMFs[0, :], this.CMFs[1, :], this.CMFs[2, :])
	}
	
    colorSpace.find_copunctuals = function find_copunctuals() {

        protan = this.find_rgb(np.array([1, 0, 0]))
        deutan = this.find_rgb(np.array([0, 1, 0]))
        tritan = this.find_rgb(np.array([0, 0, 1]))
        
        this.copunctuals = {'protan': protan, 
                            'deutan': deutan, 
                            'tritan': tritan, }
	}
        
    colorSpace.find_rgb = function find_rgb(LMS) {

        cmf = np.dot(np.linalg.inv(this.convMatrix), LMS)
        cmf[0], cmf[1], cmf[2] = this._EEcmf(cmf[0], cmf[1], cmf[2])
        out = this.TrichromaticEquation(cmf[0], cmf[1], cmf[2])
        return out
	}
	
    colorSpace._EEcmf = function _EEcmf(r_, g_, b_) {
        
        r_ *= 100. / this.EEfactors['r'] 
        g_ *= 100. / this.EEfactors['g']
        b_ *= 100. / this.EEfactors['b']
        
        return [r_, g_, b_]
	}
}