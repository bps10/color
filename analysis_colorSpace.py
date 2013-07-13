import numpy as np

from colorSpace import colorSpace

def trichromaticAnalysis(self, Lmax=560, Smax=417):
    '''
    '''
    M_lamMax = []
    volume = []
    for i in range(420, 561):
        M_lamMax.append(i)
        self.genLMS('Neitz', [Lmax, i, Smax])
        self.genConvMatrix()
        volume.append(np.linalg.det(self.convMatrix))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    
    ax.plot(M_lamMax, volume, 'k', linewidth=3)
    
    ind = np.where(np.array(M_lamMax) == 530)[0]
    ax.plot(M_lamMax[ind], volume[ind], 'g+', markersize=15, 
            markeredgewidth=3)
    
    ax.set_ylabel('color space volume')        
    ax.set_xlabel('$M_{\lambda_{max}}$') 
    
    ax.text(0.5, 0.5, 
        ('$L_{\lambda_{max}}$: ' + str(Lmax) + 'nm\n$S_{\lambda_{max}}$: '
        + str(Smax) + 'nm'), 
        fontsize=20, 
        horizontalalignment='center',
        verticalalignment='top',
        transform=ax.transAxes)
    
    plt.tight_layout()
    plt.show()

def tetrachromaticAnalysis(self, Lmax=560, Mmax=530, Smax=417):
    '''
    '''
    X_lamMax = []
    volume = []
    self.genLMS('Neitz', [Lmax, Mmax, Smax])
    
    for i in range(300, 850):
        X_lamMax.append(i)
        
        tetraSystem = self.genTetraConvMatrix(i)
        volume.append(abs(np.linalg.det(tetraSystem)))
        '''
        if i % 10 == 0:
            print i, 'nm, rank: ', np.linalg.matrix_rank(tetraSystem)
        '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[5, 5])
    
    ax.plot(X_lamMax, volume, 'k', linewidth=3)
    
    ind = np.where(np.array(X_lamMax) == 495)[0]
    ax.semilogy(X_lamMax[ind], volume[ind], 'g+', markersize=15, 
            markeredgewidth=3)
    
    ax.set_ylabel('color space volume')        
    ax.set_xlabel('$X_{\lambda_{max}}$') 
    
    ax.text(0.2, 0.95, 
        ('$L_{\lambda_{max}}$: ' + str(Lmax) + 'nm\n$M_{\lambda_{max}}$: '
        + str(Mmax) + 'nm\n$S_{\lambda_{max}}$: ' + str(Smax) + 'nm'), 
        fontsize=20, 
        horizontalalignment='center',
        verticalalignment='top',
        transform=ax.transAxes)
    
    plt.tight_layout()
    plt.show()
