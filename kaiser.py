import matplotlib.pylab as plt
import numpy as np

from colorSpace import colorSpace


def genKaiser(neitz=True):
    '''
    '''
    space = colorSpace()

    kaiser = np.genfromtxt('data/kaiser1987.csv', delimiter=",")
    # sort by subject:
    subj1 = np.where(kaiser[:, 3] == 1)
    subj2 = np.where(kaiser[:, 3] == 2)
    s1_xy = kaiser[subj1, :2][0]
    s2_xy = kaiser[subj2, :2][0]
    # solve for z:
    s1_z = 1.0 - (s1_xy[:, 0] + s1_xy[:, 1])
    s2_z = 1.0 - (s2_xy[:, 0] + s2_xy[:, 1])
    # combine:
    s1_xyz = np.zeros((len(s1_xy), 3))
    s1_xyz[:, :2] = s1_xy
    s1_xyz[:, 2] = s1_z
    s2_xyz = np.zeros((len(s2_xy), 3))
    s2_xyz[:, :2] = s2_xy
    s2_xyz[:, 2] = s2_z

    jv = space._genJuddVos()

    if neitz:
        JV_Neitz_transMatrix = space._genJuddVos2Neitz(jv)

        sub1_Neitz = np.dot(JV_Neitz_transMatrix, s1_xyz.T).T 
        sub2_Neitz = np.dot(JV_Neitz_transMatrix, s2_xyz.T).T
    else:
        sub1_Neitz = s1_xyz                     
        sub2_Neitz = s2_xyz 
                    
    return sub1_Neitz, sub2_Neitz, jv   


def plotKaiser(neitz=True, showBY=True, clip=True,
               showSub1=False, showSub2=True, stockman=False,
               series=True):
    '''
    '''
    space = colorSpace(fundamental='neitz',
                             LMSpeaks=[559.0, 530.0, 421.0])        
    sub1_Neitz, sub2_Neitz, jv = genKaiser()

    if neitz:
        space._plotColorSpace()

    else:
        space._plotColorSpace(rVal=jv[:, 0], gVal=jv[:, 1],
                             spec=space.spectrum)
        space.cs_ax.plot([jv[-1, 0], jv[0, 0]], 
                        [jv[-1, 1], jv[0, 1]], 'k-', linewidth=3)
        space.cs_ax.set_ylim([0, 0.9])
        space.cs_ax.set_xlim([-0.05, 0.8])
    if showSub1:     
        space.cs_ax.plot(sub1_Neitz[:, 0], sub1_Neitz[:, 1], 'ko', 
                        markersize=8, markeredgewidth=2,
                        markerfacecolor='w',
                        linewidth=2)
    if showSub2:
        space.cs_ax.plot(sub2_Neitz[:, 0], sub2_Neitz[:, 1], 'kx',
                        markersize=8, markeredgewidth=2, linewidth=2)

    if showBY:
        if stockman:
            neut2, RG2 = space.BY2lambda(0, 0, 1., True)
            c2 = (1, 0, 0)
            c3 = (0.5, 0.5, 0)
            neut3, RG3 = space.lambda2RG(522, False, True)
        else:
            neut2, RG2 = space.BY2lambda(1, 0, 0, True)
            c2 = (0, 0, 1)
            c3 = (0, 0.5, 0.5)
            neut3, RG3 = space.lambda2BY(522, True)
        neut1, RG1 = space.BY2lambda(0, 1., 0, True)

        c1 = (0, 1, 0)
        # plot green copunctual line
        space.cs_ax.plot([neut1[0], RG1[0]], [neut1[1], RG1[1]], 
                        '-o', c=c1, markersize=8, linewidth=2)  
        # plot red or blue copunctual depending on neitz or stockman
        space.cs_ax.plot([neut2[0], RG2[0]], [neut2[1], RG2[1]], 
                        '-o', c=c2, markersize=8, linewidth=2)  
        # plot 
        space.cs_ax.plot([neut3[0], RG3[0]], [neut3[1], RG3[1]], 
                        '-o', c=c3, markersize=8, linewidth=2)  

    if stockman and series:

        for lam in [500, 505, 510, 515]:
            neut3, RG3 = space.lambda2RG(lam, False, True)
            space.cs_ax.plot([neut3[0], RG3[0]], [neut3[1], RG3[1]], 
                '-o', c=c3, markersize=8, linewidth=2)  


    if clip is True:                
        space.cs_ax.set_xlim([-0.4, 1.2])
        space.cs_ax.set_ylim([-0.2, 1.2])
    
    space.cs_ax.set_xlabel('x', fontsize=10)
    space.cs_ax.set_ylabel('y', fontsize=10)
    
    plt.show()



if __name__ == '__main__':

    plotKaiser()