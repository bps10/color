
function spectsens(lambdaMax, opticalDensity, output,
                startWavelength, endWavelength, step) {
    /*This function returns a photopigment spectral sensitivity curve
    as defined by Carroll, McMahon, Neitz, and Neitz.

    :param LambdaMax: Wavelength peak for photopigment
    :param OpticalDensity: optical density required
    :param OutputType:  log or anti-log.  \n
                        if log, maximum data ouput is 0. \n
                        if anti-log, data output is between 0 and 1. \n
    :param StartWavelength: beginning wavelength
    :param EndWavelength: end wavelength
    :param Resolution: Number of data points

    :returns: array of sensitivity values.
    :rtype: np.array

    .. note::
       Ported from Jim K's Matlab function.
	*/
	typeof lambdaMax == "undefined" && (lambdaMax = 559);
	typeof opticalDensity == "undefined" && (opticalDensity = 0.2);
	typeof output == "undefined" && (output = 'log');
	typeof startWavelength == "undefined" && (startWavelength = 390);
	typeof endWavelength == "undefined" && (endWavelength = 750);
	typeof step == "undefined" && (step = 1);
	
    A_ = 0.417050601;
    B_ = 0.002072146;
    C_ = 0.000163888;
    D_ = -1.922880605;
    E_ = -16.05774461;
    F_ = 0.001575426;
    G_ = 5.11376E-05;
    H_ = 0.00157981;
    I_ = 6.58428E-05;
    J_ = 6.68402E-05;
    K_ = 0.002310442;
    L_ = 7.31313E-05;
    M_ = 1.86269E-05;
    N_ = 0.002008124;
    O_ = 5.40717E-05;
    P_ = 5.14736E-06;
    Q_ = 0.001455413;
    R_ = 4.217640000E-05;
    S_ = 4.800000000E-06;
    T_ = 0.001809022;
    U_ = 3.86677000E-05;
    V_ = 2.99000000E-05;
    W_ = 0.001757315;
    X_ = 1.47344000E-05;
    Y_ = 1.51000000E-05;
    Z_ = opticalDensity;

    A2 = (log10(1.0 / lambdaMax) - log10(1.0 / 558.5));

    vector = array_manip(array_pow(range(startWavelength,
                    endWavelength, step), -1.0), log10);

    con = 1.0 / Math.sqrt(2.0 * Math.PI);

    exTemp = (
		array_add(array_add(array_add(array_subtract(array_subtract(array_subtract(array_add(
			array_add(
		
		array_manip(array_add(array_multiply(array_manip(array_divide(
			array_multiply(-1.0, array_subtract(array_pow(
			10.0, array_subtract(vector, A2)), - F_)), G_), tanh), E_), -E_), log10), D_), 

		array_multiply(A_, array_manip(array_divide(array_multiply(-1, array_subtract(
			array_pow(10.0, array_subtract(vector, A2)), B_)), C_), tanh))),
			
		array_multiply(J_ * I_ * con, array_pow(Math.E, array_multiply(array_pow(
			array_divide(array_subtract(array_pow(10.0, array_subtract(vector, A2)), H_), I_),
			2.0), -0.5)))),
			
		array_multiply(M_ * L_ * con, array_pow(Math.E, array_multiply(array_pow(
			array_divide(array_subtract(array_pow(10.0, array_subtract(vector, A2)), K_), L_),
			2.0), -0.5)))),

		array_multiply(P_ * O_ * con, array_pow(Math.E, array_multiply(array_pow(
			array_divide(array_subtract(array_pow(10.0, array_subtract(vector, A2)), N_), O_),
			2.0), -0.5)))),
			
		array_multiply(S_ * R_ * con, array_pow(Math.E, array_multiply(array_pow(
			array_divide(array_subtract(array_pow(10.0, array_subtract(vector, A2)), Q_), R_),
			2.0), -0.5)))),
			
		array_multiply(V_ * U_ * con, array_pow(Math.E, array_multiply(array_pow(
			array_divide(array_subtract(array_pow(10.0, array_subtract(vector, A2)), T_), U_),
			2.0), -0.5)))),

		array_multiply(Y_ * Z_ * con, array_pow(Math.E, array_multiply(array_pow(
			array_divide(array_subtract(array_pow(10.0, array_subtract(vector, A2)), W_),X_),
			2.0), -0.5))))
		);
					
    ODTemp = array_manip(array_divide(array_subtract(1.0, array_pow(10.0, array_multiply(-1,array_multiply(array_pow(10.0, exTemp), Z_)))), (1.0 - Math.pow(10, -Z_))), log10);

    if (output.toLowerCase() === 'log') {
        extinction = exTemp;
        withOD = ODTemp;
		}
    else {
        extinction = array_pow(10.0, exTemp);
        withOD = array_pow(10.0, ODTemp);
	}
    return withOD, extinction
}

function tanh (arg) {
  return (Math.exp(arg) - Math.exp(-arg)) / (Math.exp(arg) + Math.exp(-arg));
}

function log10(val) {
  return Math.log(val) / Math.LN10;
}

function array_add(term1, term2) {
	var out_array = [],
		typeofterm1 = typeof term1,
		typeofterm2 = typeof term2;  

	if (typeofterm1 != "object" && typeofterm2 != "object") {
		throw TypeError("term1 or term2 must be an array object. Othewise use +");
	}
    if (typeofterm1 == "undefined" || typeofterm2 == "undefined") {
        throw TypeError("Must pass term1 and term2 arguments.");	
	}
	if (typeofterm1 == "object" && typeofterm2 != "object") {
		for (var i=0; i < term1.length; i++) {
			out_array.push(term1[i] + term2);
		}
	}
	else if (typeofterm1 != "object" && typeofterm2 == "object") {
		for (var i=0; i < term2.length; i++) {
			out_array.push(term1 + term2[i]);
		}
	}
	else if (typeofterm1 != "object" && typeofterm2 == "object") {
		if (term1.length != term2.length) {
			throw TypeError("term1 must be same length as term2")
		}
		for (var i=0; i < term1.length; i++) {
			out_array.push(term1[i] + term2[i]);
		}
	}		
	return out_array;	
}

function array_subtract(term1, term2) {
	var out_array = [],
		typeofterm1 = typeof term1,
		typeofterm2 = typeof term2;  

	if (typeofterm1 != "object" && typeofterm2 != "object") {
		throw TypeError("term1 or term2 must be an array object. Othewise use -");
	}
    if (typeofterm1 == "undefined" || typeofterm2 == "undefined") {
        throw TypeError("Must pass term1 and term2 arguments.");	
	}
	if (typeofterm1 == "object" && typeofterm2 != "object") {
		for (var i=0; i < term1.length; i++) {
			out_array.push(term1[i] - term2);
		}
	}
	else if (typeofterm1 != "object" && typeofterm2 == "object") {
		for (var i=0; i < term2.length; i++) {
			out_array.push(term1 - term2[i]);
		}
	}
	else if (typeofterm1 != "object" && typeofterm2 == "object") {
		if (term1.length != term2.length) {
			throw TypeError("term1 must be same length as term2")
		}
		for (var i=0; i < term1.length; i++) {
			out_array.push(term1[i] - term2[i]);
		}
	}		
	return out_array;	
}

function array_multiply(term1, term2) {
	var out_array = [],
		typeofterm1 = typeof term1,
		typeofterm2 = typeof term2;  

	if (typeofterm1 != "object" && typeofterm2 != "object") {
		throw TypeError("term1 or term2 must be an array object. Othewise use *");
	}
    if (typeofterm1 == "undefined" || typeofterm2 == "undefined") {
        throw TypeError("Must pass term1 and term2 arguments.");	
	}
	if (typeofterm1 == "object" && typeofterm2 != "object") {
		for (var i=0; i < term1.length; i++) {
			out_array.push(term1[i] * term2);
		}
	}
	else if (typeofterm1 != "object" && typeofterm2 == "object") {
		for (var i=0; i < term2.length; i++) {
			out_array.push(term1 * term2[i]);
		}
	}
	else if (typeofterm1 != "object" && typeofterm2 == "object") {
		if (term1.length != term2.length) {
			throw TypeError("term1 must be same length as term2")
		}
		for (var i=0; i < term1.length; i++) {
			out_array.push(term1[i] * term2[i]);
		}
	}		
	return out_array;	
}

function array_divide(term1, term2) {
	var out_array = [],
		typeofterm1 = typeof term1,
		typeofterm2 = typeof term2;  

	if (typeofterm1 != "object" && typeofterm2 != "object") {
		throw TypeError("term1 or term2 must be an array object. Othewise use /");
	}
    if (typeofterm1 == "undefined" || typeofterm2 == "undefined") {
        throw TypeError("Must pass term1 and term2 arguments");	
	}
	if (typeofterm1 == "object" && typeofterm2 != "object") {
		for (var i=0; i < term1.length; i++) {
			out_array.push(term1[i] / term2);
		}
	}
	else if (typeofterm1 != "object" && typeofterm2 == "object") {
		for (var i=0; i < term2.length; i++) {
			out_array.push(term1 / term2[i]);
		}
	}
	else if (typeofterm1 != "object" && typeofterm2 == "object") {
		if (term1.length != term2.length) {
			throw TypeError("term1 must be same length as term2")
		}
		for (var i=0; i < term1.length; i++) {
			out_array.push(term1[i] / term2[i]);
		}
	}		
	return out_array;	
}

function array_manip(array, func) {
	var out_array = [],
		typeofarray = typeof array,
		typeoffunc = typeof func;
		
	if (typeofarray != "object") {
		throw TypeError("Array must be an object.");
	}
	if (typeoffunc != "function") {
		throw TypeError("func must be a function.");
	}
	
	for (var i=0; i < array.length; i++) {
		out_array.push(func(array[i]));
	}
	return out_array;
}

function array_pow(base, exponent) {
	var out_array = [],
		typeofbase = typeof base,
		typeofexponent = typeof exponent;
		
	if (typeofbase != "object" && typeofexponent != "object") {
		throw TypeError("Base or exponent must be an array object. Othewise use Math.pow");
	}
	if (typeofbase == "object" && typeofexponent == "object") {
		throw TypeError("Only base or exponent can by an array object");
	}
    if (typeofbase == "undefined" || typeofexponent == "undefined") {
        throw TypeError("Must pass base and exponent arguments.");	
	}
	
	if (typeofbase == "object") {
		for (var i=0; i < base.length; i++) {
			out_array.push(Math.pow(base[i], exponent));
		}
	}
	else if (typeofexponent == "object") {
		for (var i=0; i < exponent.length; i++) {
			out_array.push(Math.pow(base, exponent[i]));
		}
	}
	return out_array;
}

function range(start, end, step) {
	/* from http://stackoverflow.com/questions/3895478/does-javascript-have-a-range-equivalent*/
	
    var range = [];
    var typeofStart = typeof start;
    var typeofEnd = typeof end;

    if (step === 0) {
        throw TypeError("Step cannot be zero.");
    }

    if (typeofStart == "undefined" || typeofEnd == "undefined") {
        throw TypeError("Must pass start and end arguments.");
    } else if (typeofStart != typeofEnd) {
        throw TypeError("Start and end arguments must be of same type.");
    }

    typeof step == "undefined" && (step = 1);

    if (end < start) {
        step = -step;
    }

    if (typeofStart == "number") {

        while (step > 0 ? end >= start : end <= start) {
            range.push(start);
            start += step;
        }

    } else if (typeofStart == "string") {

        if (start.length != 1 || end.length != 1) {
            throw TypeError("Only strings with one character are supported.");
        }

        start = start.charCodeAt(0);
        end = end.charCodeAt(0);

        while (step > 0 ? end >= start : end <= start) {
            range.push(String.fromCharCode(start));
            start += step;
        }

    } else {
        throw TypeError("Only string and number types are supported");
    }

    return range;

}
