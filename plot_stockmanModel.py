

def plotStockmanAnalysis(space):
    '''
    '''
    stage2, stage3 = space.genStockmanAnalysis()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    for key in stage2:
        ax.plot(space.spectrum, stage2[key])

    ax.set_xlim([space.spectrum[0], space.spectrum[-1]])
    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('sensitivity')
    plt.tight_layout()
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    for key in stage3:
        ax.plot(space.spectrum, stage3[key], c=key)

    ax.set_xlim([space.spectrum[0], space.spectrum[-1]])
    ax.set_xlabel('wavelength (nm)')
    ax.set_ylabel('sensitivity')
    plt.tight_layout()
    plt.show()

    # Unique green series plot
    fig1 = plt.figure()
    fig2 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax1, ['left', 'bottom'], Nticks=[5, 5])
    pf.TufteAxis(ax2, ['left', 'bottom'], Nticks=[5, 5])

    for j in [0.04, 0.34, 0.94, 2, 4, 6]:
        stage2, stage3 = space.genStockmanAnalysis(j)
        ax1.plot(space.spectrum, stage3['blue'], c='b', alpha=0.7)
        ax2.plot(space.spectrum, space.Snorm -
            j / 10 * (space.Lnorm + (0.5 * space.Mnorm)), 
            c='b', alpha=0.7)

    ax1.plot(space.spectrum, np.zeros(len(space.spectrum)), 'k',
        linewidth=1)
    ax2.plot(space.spectrum, np.zeros(len(space.spectrum)), 'k',
        linewidth=1)
    ax1.set_xlim([space.spectrum[0], 650])
    ax1.set_ylim([-0.7, 1.4])
    ax1.set_xlabel('wavelength (nm)')
    ax1.set_ylabel('sensitivity')

    ax2.set_xlim([space.spectrum[0], 700])
    ax2.set_ylim([-0.9, 1.2])
    ax2.set_xlabel('wavelength (nm)')
    ax2.set_ylabel('sensitivity')

    fig1.tight_layout()
    fig2.tight_layout()
    plt.show()

def main(args):
    '''
    '''
    if args.Stockman:
        plotStockmanAnalysis()


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description="Color Space: display Neitz or Stockman\
        derived color spaces")

    parser.add_argument("-m", "--Stockman", action="store_true",
                        help="plot Stockman model.")
    args = parser.parse_args()
    main(args)

