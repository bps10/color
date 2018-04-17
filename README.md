# COLOR

*These scripts simulate aspects of human color vision based on the Neitz theory of color opponency.*

## Motivation

A parsimonious account of the distribution of monochromatic lights identified as unique, most notably green, has eluded models of color perception.  The standard model in the field consists of three stages predicated upon opponent signals between L- versus M-cones and S-cones versus L- plus M-cone signals.  This model fails to predict the variability in unique hue matching experiments.  We have developed a three stage zone model that takes into account L/M cone ratio to produce distributions of unique hues with accurate statistical characteristics.  The first stage of our model is the known sensitivity functions of the cone photoreceptors.  The second stage is composed of a center of either an L- or M-cone opposed to a surround of all three cones with varying contribution.  In the third stage these opponent signals are weighted based on the probability of occurrence and summated to construct valence curves for the blue-yellow and red-green mechanisms.  Varying the contribution of L- and M-cones to the surround based on an observed distribution of L:M (Carroll et al 2000), we accurately predict the distributions of unique blue, yellow and green.  Blue and yellow form narrow distributions around 470nm and 578nm, respectively, while green broadly distributes between 495 and 555nm and closely fits the unique green distribution of (Volbrecht et al 1997). We conclude that our model provides a convincing description of color processing and offers insight into the contribution of L/M ratio in color perception. 


This work was described in the publication: 

> [Neurobiological hypothesis of color appearance and hue perception](https://www.osapublishing.org/DirectPDFAccess/59C4123D-D542-852A-7CFA1CDE2D615555_279354/josaa-31-4-A195.pdf?da=1&id=279354&seq=0&mobile=no). Schmidt, Neitz, Neitz. 2014. Journal of the Optical Society of America.

## Installation

### Dependencies

This module is written in the open source programming language [python][py]. Additionally, two libraries are necessary to run this module: [numpy][np] and [matplotlib][mpl]. Both of these can be easily installed with the `pip` command if it is installed on your computer. Additionally, numerous scientific computing bundles of python come with these packages. These bundles are free and come with easy to use gui installers. The [enthought][enth] suite or the [pythonxy][pyxy] (windows and linux machines only) are two excellent examples.

[enth]: https://www.enthought.com/
[pyxy]: https://code.google.com/p/pythonxy/
[np]: http://scipy.org
[mpl]: http://matplotlib.org
[py]: http://python.org

### Download

To install these scripts [git](http://git-scm.com/) must be installed on your computer. Then from the command line run

```bash
git clone https://github.com/bps10/color
```

This will download the files to a directory called color. You then need to initialize the submodule, base, which contains some general plotting routines and the spectral sensitivity functions.

```bash
git submodule init
git submodule update
```

You can then run

```bash
python plot_colorModel.py -h
```

To display options for plotting the color model. This same command can be used to display options for `plot_colorSpace.py`, `plot_sensitivity.py` and `plot_stockmanModel.py`. 


