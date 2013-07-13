import numpy as np

from colorSpace import colorSpace

def genKaiser(self, neitz=False):
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

    
def main(args):
    '''
    '''
    if args.Kaiser:
        plotKaiser(neitz=True, stockman=True)  


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description="Color Space: display Neitz or Stockman\
        derived color spaces")
    
    parser.add_argument("-k", "--Kaiser", action="store_true",
                        help="plot Kaiser data in Neitz or CIE space")   
    args = parser.parse_args()
    main(args)