

function spectsens(LambdaMax=559, OpticalDensity=0.2, Output='log',
                StartWavelength=380, EndWavelength=780, Res=1000) {
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
	if (arguments.length < 1) {LambdaMax = 559;}
	if (arguments.length < 2) {OpticalDensity = 0.2;}
	if (arguments.length < 3) {Output = 'log';}
	if (arguments.length < 4) {StartWavelength = 380;}
	if (arguments.length < 5) {EndWavelength = 780;}
	if (arguments.length < 6) {Res = 1000;}
		
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
    Z_ = OpticalDensity; //+ 0.00000001;

    A2 = (log10(1.0 / LambdaMax) - log10(1.0 / 558.5));

    vector = log10(np.arange(StartWavelength,
                    EndWavelength + Res, Res) ** -1.0);

    con = 1.0 / Math.sqrt(2.0 * Math.PI);

    exTemp = (log10(-E_ + E_ * tanh(-1 * Math.pow(10.0, vector - A2) - F_) / G_)) + 
				D_ + A_ * tanh((-1 * Math.pow(10.0, vector - A2) - B_) / C_) -
				
              (J_ / I_ * con * Math.pow(Math.E, -0.5 *
              Math.pow((Math.pow(10.0, vector - A2) - H_) / I_, 2.0))) -

              (M_ / L_ * con * Math.pow(Math.E, -0.5 *
              Math.pow((Math.pow(10.0, vector - A2) - K_) / L_, 2.0))) -

              (P_ / O_ * con * Math.pow(Math.E, -0.5 *
              Math.pow((Math.pow(10.0, vector - A2) - N_) / O_, 2.0))) +

              (S_ / R_ * con * Math.pow(Math.E, -0.5 *
              Math.pow((Math.pow(10.0, vector - A2) - Q_) / R_, 2.0))) +
			  
              (V_ / U_ * con * Math.pow(Math.E, -0.5 *
              Math.pow((Math.pow(10.0, vector - A2) - T_) / U_, 2.0)) / 10.0) +

              (Y_ / X_ * con * Math.pow(Math.E, -0.5 *
              Math.pow((Math.pow(10.0, vector - A2) - W_) / X_, 2.0)) / 100.0);
			  
    ODTemp = log10((1.0 - Math.pow(10.0, -1 * (Math.pow(10.0, exTemp) *
                        Z_)) / (1.0 - Math.pow(10, -Z_)));

    if Output.toLowerCase() === 'log' {
        extinction = exTemp;
        withOD = ODTemp;
		}
    else {
        extinction = Math.pow(10.0, exTemp);
        withOD = Math.pow(10.0, ODTemp);
	}
    return withOD, extinction
}

function tanh (arg) {
  return Math.exp(arg) - Math.exp(-arg)) / (Math.exp(arg) + Math.exp(-arg);
}

function log10(val) {
  return Math.log(val) / Math.LN10;
}