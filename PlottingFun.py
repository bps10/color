import numpy as np
import matplotlib.pylab as plt
import matplotlib.ticker as tic

## To DO:
## 1. Help files.
## 2. More functions.


### Plotting cosmetics ###
##########################


def histOutline(histIn,binsIn):
    """

    This will take a histogram and return an outline of the histogram.

    :param histIn: histogrammed data.
    :type histIn: np.array
    :param binsIn: bins used to histogram data.
    :type binsIn: np.array

    :returns: bins and data for plotting a histogram outline.


    """
    stepSize = binsIn[1] - binsIn[0]

    bins = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    data = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    for bb in range(len(binsIn)):
        bins[2*bb + 1] = binsIn[bb]
        bins[2*bb + 2] = binsIn[bb] + stepSize
        if bb < len(histIn):
            data[2*bb + 1] = histIn[bb]
            data[2*bb + 2] = histIn[bb]

    bins[0] = bins[1]
    bins[-1] = bins[-2]
    data[0] = 0
    data[-1] = 0
    return (bins, data)


def simpleaxis(ax):
    """

    Change the axes of a plot to have simple axes.

    :param ax: handle to the axes.

    :returns: altered axes.

    .. note::
       Taken from here: `stackoverflow
       <http://stackoverflow.com/questions/925024/how-can-i-remove-the-top-and-right-axis-in-matplotlib>`_


    """
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def centerAxes(ax):
    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_smart_bounds(True)
    ax.spines['bottom'].set_smart_bounds(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

def TufteAxis(ax, spines, Nticks = None, integer='on'):
    """

    Change the axes of a plot to have axes that are offset, inspired by Edward Tufte.

    :param ax: handle to the axes.
    :param spines: list of spines you wish to display. Others will be hidden.
    :param Nticks: number of ticks you wish to show.
    :type Nticks: list of int
    :param integer: plot only integers.

    :returns: altered axes.

    .. note::
       * Adapted from a `Matplotlib Example`_
       * This function should work if called prior to plot input.


    **Usage**

    .. code-block:: python
       :emphasize-lines: 4

       fig = plt.figure(figsize=(8,6))
       ax = fig.add_subplot(111)

       pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5,7])

       ax.plot(np.arange(0,100), np.arange(0,100), 'b-', linewidth =2)

       plt.tight_layout()
       plt.show()


    **Exmaple**

    *simple plot without TufteAxis call*

    .. plot:: pyplots/TufteAxis1.py
       :include-source:

    *simple plot with TufteAxis call*

    .. plot:: pyplots/TufteAxis2.py
       :include-source:

    .. _Matplotlib Example:
       http://matplotlib.org/examples/pylab_examples/spine_placement_demo.html
    """
    for loc, spine in ax.spines.iteritems():
        if loc in spines:
            spine.set_position(('outward',10)) # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none') # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
        if integer.lower == 'on':
            if Nticks == None:
                ax.yaxis.set_major_locator(tic.MaxNLocator(integer=True))
            else:
                ax.yaxis.set_major_locator(tic.MaxNLocator(Nticks[1] + 1,
                                                           integer=True))
        else:
            if Nticks == None:
                ax.yaxis.set_major_locator(tic.MaxNLocator(integer=False))
            else:
                ax.yaxis.set_major_locator(tic.MaxNLocator(Nticks[1] + 1,
                                                           integer=False))
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
        if integer.lower == 'on':
            if Nticks == None:
                ax.xaxis.set_major_locator(tic.MaxNLocator(integer=True))
            else:
                ax.xaxis.set_major_locator(tic.MaxNLocator(Nticks[0] + 1,
                                                           integer=True))
        else:
            if Nticks == None:
                ax.xaxis.set_major_locator(tic.MaxNLocator(integer=False))
            else:
                ax.xaxis.set_major_locator(tic.MaxNLocator(Nticks[0] + 1,                                                           integer=False))
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])

    if 'top' in spines:
        ax.xaxis.set_ticks_position('top')
        if integer.lower == 'on':
            if Nticks == None:
                ax.xaxis.set_major_locator(tic.MaxNLocator(integer=True))
            else:
                ax.xaxis.set_major_locator(tic.MaxNLocator(Nticks[0] + 1,
                                                           integer=True))
        else:
            if Nticks == None:
                ax.xaxis.set_major_locator(tic.MaxNLocator(integer=False))
            else:
                ax.xaxis.set_major_locator(tic.MaxNLocator(Nticks[0] + 1,
                                                           integer=False))




def SciNoteAxis(gca_handle,spines):
    """

    Force scientific notation.

    :param gca_handle: plt.gca handle
    :param spines: indicate what spines [x and/or y] you wish to force sci note on.
    :type spines: list of strings

    .. note::
       Not thoroughly tested though it does work in a few isolated cases.

    """
    if 'y' in spines:

        gca_handle.yaxis.major.formatter.set_powerlimits((-1, 0))
        #t.ticklabel_format(style='sci', axis='y')

    if 'x' in spines:

        gca_handle.xaxis.major.formatter.set_powerlimits((-1, 0))
        #t.ticklabel_format(style='sci', axis='x')



def AxisFormat(FONTSIZE = 22, TickSize = 10, TickDirection = 'out'):
    """

    Format axes to standard design.

    :param FONTSIZE: desired fontsize of all fonts.
    :type FONTSIZE: int
    :param TickSize: size of ticks in pxls.
    :type TickSize: int
    :param TickDirection: decide whether to plot ticks in or out
    :type TickDirection: str

    :returns: altered axes.

    .. note::
       * This function should work if called prior to plot input.

    **Usage**

    .. code-block:: python
       :emphasize-lines: 3

       fig = plt.figure(figsize=(8,6))
       ax = fig.add_subplot(111)
       pf.AxisFormat()

       ax.plot(np.arange(0,100), np.arange(0,100), 'b-', linewidth =2)

       plt.tight_layout()
       plt.show()

    **Exmaple**


    *simple plot without AxisFormat call*

    .. plot:: pyplots/AxisFormatDemo1.py
       :include-source:

    *simple plotnwith AxisFormat call*

    .. plot:: pyplots/AxisFormatDemo2.py
       :include-source:

    """
    font = {'weight': 'norm', 'size': FONTSIZE}
    legend = {'frameon': False}
    ticks = {'direction': TickDirection, 'major.size': TickSize,
             'minor.size': TickSize - 2}

    plt.rc('font', **font)
    plt.rc('legend', **legend)
    plt.rc('xtick', **ticks)
    plt.rc('ytick', **ticks)
