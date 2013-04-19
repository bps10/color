import numpy as np
import matplotlib as mpl

def clip_rgb_color (rgb_color):
    '''Convert a linear rgb color (nominal range 0.0 - 1.0), into a displayable
    irgb color with values in the range (0 - 255), clipping as necessary.

    The return value is a tuple, the first element is the clipped irgb color,
    and the second element is a tuple indicating which (if any) clipping processes were used.
    '''
    clipped_chromaticity = False
    clipped_intensity = False

    rgb = rgb_color.copy()

    # add enough white to make all rgb values nonnegative
    # find max negative rgb (or 0.0 if all non-negative), we need that much white
    rgb_min = np.min(0.0, np.min(rgb))
    # get max positive component
    rgb_max = np.max(rgb)
	# get scaling factor to maintain max rgb after adding white
    scaling = 1.0
    if rgb_max > 0.0:
        scaling = rgb_max / (rgb_max - rgb_min)
    # add enough white to cancel this out, maintaining the maximum of rgb
    if rgb_min < 0.0:
        rgb [0] = scaling * (rgb [0] - rgb_min)
        rgb [1] = scaling * (rgb [1] - rgb_min)
        rgb [2] = scaling * (rgb [2] - rgb_min)
        clipped_chromaticity = True

    # clip intensity if needed (rgb values > 1.0) by scaling
    rgb_max = np.max(rgb)
    # we actually don't overflow until 255.0 * intensity > 255.5, so instead of 1.0 use ...
    intensity_cutoff = 1.0 + (0.5 / 255.0)
    if rgb_max > intensity_cutoff:
        # must scale intensity, so max value is intensity_cutoff
        scaling = intensity_cutoff / rgb_max
        rgb *= scaling
        clipped_intensity = True

    # gamma correction
    #for index in xrange (0, 3):
    #    rgb [index] = display_from_linear_component (rgb [index])

    # ensure that values are in the range 0-1
    thresh = rgb > 1
    rgb[thresh] = 1
    thresh = rgb < 0
    rgb[thresh] = 0
    
    return (rgb, (clipped_chromaticity, clipped_intensity))    


def rgb_patch_plot (
    ax,
    rgb_colors,
    color_names,
    #title,
    #filename = None,
    patch_gap = 0.05,
    num_across = 60):
    '''Draw a set of color patches, specified as linear rgb colors.'''
    
    def draw_patch (x0, y0, color, patch_gap):
        '''Draw a patch of color.'''
        # patch relative vertices
        m = patch_gap
        omm = 1.0 - m
        poly_dx = [m, m, omm, omm]
        poly_dy = [m, omm, omm, m]
        # construct vertices
        poly_x = [ x0 + dx_i for dx_i in poly_dx ]
        poly_y = [ y0 + dy_i for dy_i in poly_dy ]
        ax.fill (poly_x, poly_y, mpl.colors.rgb2hex(color))


    # make plot with each color with one patch
    #plt.clf()
    num_colors = rgb_colors.shape[1]
    for i in xrange (0, num_colors):
        (iy, ix) = divmod (i, num_across)
        
        # get color as a displayable string
        draw_patch (float (ix), float (-iy), rgb_colors[:,i], patch_gap)
    #plt.axis ('off')
    #plt.title (title)
    #if filename is not None:
    #    print 'Saving plot %s' % str (filename)
    #    plt.savefig (filename)
        
def shark_fin_plot (ax):
    '''Draw the 'shark fin' CIE chromaticity diagram of the pure spectral lines (plus purples) in xy space.'''
    # get array of (approximate) colors for the boundary of the fin
    import math
    xyz_list = ciexyz.get_normalized_spectral_line_colors (brightness=1.0, num_purples=200, dwl_angstroms=2)
    # get normalized colors
    xy_list = xyz_list.copy()
    (num_colors, num_cols) = xy_list.shape
    for i in xrange (0, num_colors):
        colormodels.xyz_normalize (xy_list [i])
    # get phosphor colors and normalize
    red   = colormodels.PhosphorRed
    green = colormodels.PhosphorGreen
    blue  = colormodels.PhosphorBlue
    white = colormodels.PhosphorWhite
    colormodels.xyz_normalize (red)
    colormodels.xyz_normalize (green)
    colormodels.xyz_normalize (blue)
    colormodels.xyz_normalize (white)

    def get_direc_to_white (xyz):
        '''Get unit vector (xy plane) in direction of the white point.'''
        direc = white - xyz
        mag = math.hypot (direc [0], direc [1])
        if mag != 0.0:
            direc /= mag
        return (direc[0], direc[1])

    # plot
    #pylab.clf ()
    # draw color patches for point in xy_list
    s = 0.025     # distance in xy plane towards white point
    for i in xrange (0, len (xy_list)-1):
        x0 = xy_list [i][0]
        y0 = xy_list [i][1]
        x1 = xy_list [i+1][0]
        y1 = xy_list [i+1][1]
        # get unit vectors in direction of white point
        (dir_x0, dir_y0) = get_direc_to_white (xy_list [i])
        (dir_x1, dir_y1) = get_direc_to_white (xy_list [i+1])
        # polygon vertices
        poly_x = [x0, x1, x1 + s*dir_x1, x0 + s*dir_x0]
        poly_y = [y0, y1, y1 + s*dir_y1, y0 + s*dir_y0]
        # draw (using full color, not normalized value)
        color_string = colormodels.irgb_string_from_rgb (
            colormodels.rgb_from_xyz (xyz_list [i]))
        ax.fill (poly_x, poly_y, color_string, edgecolor=color_string)
    # draw the curve of the xy values of the spectral lines and purples
    ax.plot (xy_list [:,0], xy_list [:,1], color='#808080', linewidth=3.0)
    # draw monitor gamut and white point
    ax.plot ([red  [0], green[0]], [red  [1], green[1]], 'o-', color='k')
    ax.plot ([green[0], blue [0]], [green[1], blue [1]], 'o-', color='k')
    ax.plot ([blue [0], red  [0]], [blue [1], red  [1]], 'o-', color='k')
    ax.plot ([white[0], white[0]], [white[1], white[1]], 'o-', color='k')
    # label phosphors
    dx = 0.01
    dy = 0.01
    ax.text (red   [0] + dx, red   [1], 'Red',   ha='left',   va='center')
    ax.text (green [0], green [1] + dy, 'Green', ha='center', va='bottom')
    ax.text (blue  [0] - dx, blue  [1], 'Blue',  ha='right',  va='center')
    ax.text (white [0], white [1] + dy, 'White', ha='center', va='bottom')
    # titles etc
    #ax.axis ([0.0, 0.85, 0.0, 0.85])
    ax.xlabel (r'CIE $x$')
    ax.ylabel (r'CIE $y$')
    ax.title (r'CIE Chromaticity Diagram')
    filename = 'ChromaticityDiagram'
    print 'Saving plot %s' % (str (filename))
    ax.savefig (filename)