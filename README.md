# COLOR

*These scripts simulate aspects of human color vision based on the Neitz theory of color opponency.*

## Motivation

De Valois and De Valois [Vis. Res. 33, 1053 (1993)] showed that to explain hue appearance, S-cone signals have to
be combined with M versus L opponent signals in two different ways to produce red–green and yellow–blue axes,
respectively. Recently, it has been shown that color appearance is normal for individuals with genetic mutations
that block S-cone input to blue-ON ganglion cells. This is inconsistent with the De Valois hypothesis in which
S-opponent konio-geniculate signals are combined with L-M signals at a third processing stage in cortex. Instead,
here we show that color appearance, including individual differences never explained before, are predicted by a
model in which S-cone signals are combined with L versus M signals in the outer retina. 

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


