################################################################################
# History:
#      v1 - A. Eigenbrot Jan. 2015
#      v1.1 - Feb. 2015 - Added Fits image capabilities
#
################################################################################
"""
************
GradPak_plot
************

This module contains an API for plotting data taken with GradPak. At
it's heart it contains hardcoded information about the location of all
the GradPak fibers on the sky, which can then be used to make nice
plots with pyplot's Patch Collection object.

User interaction should be limited to the plot*() functions, all of
which, at their most basic level, take in an array of length 109 (the
number of GradPak fibers) and return an pyplot Axes object containing
a plot that can be further modified by the user. Specific details
about the plot*() functions and their copious advanced options can be
found below.

Useful Example
--------------

An extremely basic example can be used to plot the GradPak fibers
color coded by fiber number::

 >>> ax = GradPak_plot.plot(np.arange(109))
 >>> ax.figure.show()

More options, including synergy with :mod:`GradPak_bin`, FITS image
overplotting, rotations, and more can be found in :func:`plot`, which is
really the only function you should ever use.

Functions
---------
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pyfits
import scipy.interpolate as spi
import scipy.ndimage.interpolation as spndi
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection

plt.ioff() #Fuck you, new versions of pyplot
tau = 2 * np.pi #Fuck you pi, your days are numbered

def GradPak_patches():
    '''
    Produce the "data" of where the GradPack fibers live. Because we're all
    about plotting, these data are stored as pyplot Patches. Fiber locations
    are given in arcsec relative to the central, bottom row fiber (fiber 10).

    Returns a tuple of fiber number and the corresponding Patch
    '''

    raw_patches = [
        Circle((48.1542, 85.5171), radius=0.9374),
        Circle((-48.1542, 85.5171), radius=0.9374),
        Circle((15.8139, 0.0000), radius=0.9374),
        Circle((13.5548, 0.0000), radius=0.9374),
        Circle((11.2957, 0.0000), radius=0.9374),
        Circle((9.0365, 0.0000), radius=0.9374),
        Circle((6.7774, 0.0000), radius=0.9374),
        Circle((4.5183, 0.0000), radius=0.9374),
        Circle((2.2591, 0.0000), radius=0.9374),
        Circle((0.0000, 0.0000), radius=0.9374),
        Circle((-2.2591, 0.0000), radius=0.9374),
        Circle((-4.5183, 0.0000), radius=0.9374),
        Circle((-6.7774, 0.0000), radius=0.9374),
        Circle((-9.0365, 0.0000), radius=0.9374),
        Circle((-11.2957, 0.0000), radius=0.9374),
        Circle((-13.5548, 0.0000), radius=0.9374),
        Circle((-15.8139, 0.0000), radius=0.9374),
        Circle((-42.3236, 85.5171), radius=0.9374),
        Circle((42.3236, 85.5171), radius=0.9374),
        Circle((48.9201, 90.8997), radius=1.4061),
        Circle((15.7343, 3.1159), radius=1.4061),
        Circle((12.2378, 3.1159), radius=1.4061),
        Circle((8.7413, 3.1159), radius=1.4061),
        Circle((5.2448, 3.1159), radius=1.4061),
        Circle((1.7483, 3.1159), radius=1.4061),
        Circle((-1.7483, 3.1159), radius=1.4061),
        Circle((-5.2448, 3.1159), radius=1.4061),
        Circle((-8.7413, 3.1159), radius=1.4061),
        Circle((-12.2378, 3.1159), radius=1.4061),
        Circle((-15.7343, 3.1159), radius=1.4061),
        Circle((-48.9201, 90.8997), radius=1.4061),
        Circle((41.5578, 90.8997), radius=1.4061),
        Circle((15.7343, 6.6124), radius=1.4061),
        Circle((12.2378, 6.6124), radius=1.4061),
        Circle((8.7413, 6.6124), radius=1.4061),
        Circle((5.2448, 6.6124), radius=1.4061),
        Circle((1.7483, 6.6124), radius=1.4061),
        Circle((-1.7483, 6.6124), radius=1.4061),
        Circle((-5.2448, 6.6124), radius=1.4061),
        Circle((-8.7413, 6.6124), radius=1.4061),
        Circle((-12.2378, 6.6124), radius=1.4061),
        Circle((-15.7343, 6.6124), radius=1.4061),
        Circle((-41.5578, 90.8997), radius=1.4061),
        Circle((49.7300, 97.7793), radius=1.8748),
        Circle((15.9124, 10.8720), radius=1.8748),
        Circle((11.3660, 10.8720), radius=1.8748),
        Circle((6.8196, 10.8720), radius=1.8748),
        Circle((2.2732, 10.8720), radius=1.8748),
        Circle((-2.2732, 10.8720), radius=1.8748),
        Circle((-6.8196, 10.8720), radius=1.8748),
        Circle((-11.3660, 10.8720), radius=1.8748),
        Circle((-15.9124, 10.8720), radius=1.8748),
        Circle((-49.7300, 97.7793), radius=1.8748),
        Circle((40.7478, 97.7793), radius=1.8748),
        Circle((15.9124, 15.4184), radius=1.8748),
        Circle((11.3660, 15.4184), radius=1.8748),
        Circle((6.8196, 15.4184), radius=1.8748),
        Circle((2.2732, 15.4184), radius=1.8748),
        Circle((-6.8196, 15.4184), radius=1.8748),
        Circle((-11.3660, 15.4184), radius=1.8748),
        Circle((-15.9124, 15.4184), radius=1.8748),
        Circle((-40.7478, 97.7793), radius=1.8748),
        Circle((45.2389, 82.7434), radius=2.3435),
        Circle((16.6201, 20.6997), radius=2.3435),
        Circle((11.0801, 20.6997), radius=2.3435),
        Circle((5.5400, 20.6997), radius=2.3435),
        Circle((-0.0000, 20.6997), radius=2.3435),
        Circle((-5.5400, 20.6997), radius=2.3435),
        Circle((-11.0801, 20.6997), radius=2.3435),
        Circle((-16.6201, 20.6997), radius=2.3435),
        Circle((-45.2389, 82.7434), radius=2.3435),
        Circle((16.6201, 26.2397), radius=2.3435),
        Circle((11.0801, 26.2397), radius=2.3435),
        Circle((5.5400, 26.2397), radius=2.3435),
        Circle((-0.0000, 26.2397), radius=2.3435),
        Circle((-5.5400, 26.2397), radius=2.3435),
        Circle((-11.0801, 26.2397), radius=2.3435),
        Circle((-16.6201, 26.2397), radius=2.3435),
        Circle((45.2389, 88.2909), radius=2.3435),
        Circle((16.6201, 31.7797), radius=2.3435),
        Circle((11.0801, 31.7797), radius=2.3435),
        Circle((5.5400, 31.7797), radius=2.3435),
        Circle((-0.0000, 31.7797), radius=2.3435),
        Circle((-5.5400, 31.7797), radius=2.3435),
        Circle((-11.0801, 31.7797), radius=2.3435),
        Circle((-16.6201, 31.7797), radius=2.3435),
        Circle((-45.2389, 88.2909), radius=2.3435),
        Circle((45.2389, 94.4215), radius=2.8122),
        Circle((16.7092, 38.1297), radius=2.8122),
        Circle((10.0255, 38.1297), radius=2.8122),
        Circle((3.3418, 38.1297), radius=2.8122),
        Circle((-3.3418, 38.1297), radius=2.8122),
        Circle((-10.0255, 38.1297), radius=2.8122),
        Circle((-16.7092, 38.1297), radius=2.8122),
        Circle((-45.2389, 94.4215), radius=2.8122),
        Circle((16.7092, 44.8133), radius=2.8122),
        Circle((10.0255, 44.8133), radius=2.8122),
        Circle((3.3418, 44.8133), radius=2.8122),
        Circle((-3.3418, 44.8133), radius=2.8122),
        Circle((-10.0255, 44.8133), radius=2.8122),
        Circle((-16.7092, 44.8133), radius=2.8122),
        Circle((-45.2389, 101.1361), radius=2.8122),
        Circle((16.7092, 51.4970), radius=2.8122),
        Circle((10.0255, 51.4970), radius=2.8122),
        Circle((-16.7092, 51.4970), radius=2.8122),
        Circle((-3.3418, 51.4970), radius=2.8122),
        Circle((-10.0255, 51.4970), radius=2.8122),
        Circle((3.3418, 51.4970), radius=2.8122),
        Circle((45.2389, 101.1361), radius=2.8122)]

    patch_list = [[i+1,p] for i,p in enumerate(raw_patches)]

    return np.array(patch_list)

def get_binned_patches(header):
    """Produce patches from the header of a GradPak multispec data
    file produced by :mod:`GradPak_bin`.

    The output is interchangeable with the output of :func:`GradPak_patches`
    """
    
    patches = GradPak_patches()[:,1]
    patch_list = []

    for i in range(109):

        try:
            fibers = header['BIN{:03}F'.format(i+1)]
            pos = header['BIN{:03}P'.format(i+1)]
        except KeyError:
            break

        r = patches[int(fibers.split(' ')[0]) - 1].get_radius()
        xpos = float(pos.split(' ')[0])
        ypos = float(pos.split(' ')[1])
        patch_list.append([i+1,Circle((xpos,ypos), radius=r)])

    return np.array(patch_list)

def fill_fibers_with_bins(header, values):
    """Given a FITS header containing information about the fibers in
    each bin return a list of values where each individual GradPak
    fiber is assigned the value of its corresponding bin.

    This function is useful if you want to plot every GradPak fiber,
    but fill them with binned data values.
    """
    
    newvalues = np.zeros(109) + np.nan
    for i in range(109):
        
        try:
            fibers = header['BIN{:03}F'.format(i+1)]
        except KeyError:
            break

        for f in fibers.split():
            newvalues[int(f)-1] = values[i]

    return newvalues

def exclude_bins(header, exclude):
    """Convert an exclusion vector using binned apertures to an
    exclusion vector using raw GradPak fibers.

    When using a binheader the user will specify excluded data by
    *aperture* number. This function converts these numbers to the
    range of underlying GradPak *fiber* numbers, which is what
    :func:`prep_patches` and its ilk expect.
    """

    new_exclude = []
    
    for ap in exclude:
        try:
            fibers = header['BIN{:03}F'.format(ap)]
        except KeyError:
            print "WARNING: Tried to exclude ap {}, but this can't be done".format(int(ap))
            continue
            
        for f in fibers.split():
            new_exclude.append(int(f))

    return new_exclude

def get_bin_boxes(header, patches, pval):
    """Given a FITS header containing binning information return a
    list of pyplot Polygon objects that can be used to plot boxes
    around all the fibers in a given bin.
    """
    
    boxlist = []
    bval = np.array([])
    for i in range(109):
        try:
            fibers = header['BIN{:03}F'.format(i+1)]
        except KeyError:
            break

        fiblist = fibers.split(' ')
        if len(fiblist) == 1: continue

        exist_fib = False
        for f in fiblist:
            if int(f) in patches[:,0]:
                exist_fib = True
                idx1 = np.where(patches[:,0] == int(fiblist[0]))[0][0]
                break

        if not exist_fib:
            print 'Bin {} was not found (excluded?), skipping'.format(i+1)
            continue

        # A special case becase 108 and 105 are out of order.
        # 105 is on the edge and should therefore always be x2, y2
        if '105' in fiblist:
            idx2 = np.where(patches[:,0] == 105)[0][0]
        else:
            idx2 = np.where(patches[:,0] == int(fiblist[-1]))[0][0]

        bval = np.r_[bval, pval[idx1]]

        x1, y1 = patches[idx1,1].center
        x2, y2 = patches[idx2,1].center
        r = patches[idx1,1].get_radius()
        theta = np.arctan((y2 - y1)/(x2 - x1))
        ex = np.cos(theta)*r
        ey = np.sin(theta)*r
        dx = np.sin(theta)*r
        dy = np.cos(theta)*r

        ex = 0
        ey = 0

        p1 = [x1 - ex + dx, y1 - ey - dy]
        p2 = [x2 + ex + dx, y2 + ey - dy]
        p3 = [x2 + ex - dx, y2 + ey + dy]
        p4 = [x1 - ex - dx, y1 - ey + dy]
        
        boxlist.append(Polygon(np.array([p1,p2,p3,p4]), closed=True))
    
    return boxlist, bval

def transform_patches(patches, pa=0, center=[0,0], reffiber=105, scale=1.,
                      refpatches=None):
    """
    Rotate and shift the centers of the GradPak Patches.

    Parameters
    ----------
    
    patches : list of matplotlib.patches.Circle
        GradPak patches produced by GradPak_patches()

    pa : float
        Position angle, in degrees, of GradPak on sky.

    center : list 
        Tuple or length 2 list containing the position, in decimal 
        degrees, of the fiber specified in **reffiber**

    reffiber : int
        The fiber to place at center

    scale : float
        Radial scale factor, in px/arcsec. This is needed because 
        Patches are plotted in pixel space.

    refpatches : None or list of matplotlib.patches.Circle
        A full set of all GradPak patches that are used to compute the
        transformed locations of binned patches.
        
    Returns
    -------

    patches : list of pyplot.Patches
        The shifted, rotated, and scaled patches

    """
    if refpatches is None:
        refpatches = np.copy(patches)
        recurse = False
    else:
        recurse = True
    refcenter = np.array(refpatches[reffiber - 1,1].center) #in arcsec
    parad = -1.*pa*tau/360.      #Radians
    decrad = center[1]*tau/360.  #Radians
    rotrefx = refcenter[0]*np.cos(parad) - refcenter[1]*np.sin(parad)
    rotrefy = refcenter[0]*np.sin(parad) + refcenter[1]*np.cos(parad)
    for c in patches[:,1]:
        startx = c.center[0] #in arcsec
        starty = c.center[1] 
        
        #Rotate
        rotx = startx*np.cos(parad) - starty*np.sin(parad)
        roty = startx*np.sin(parad) + starty*np.cos(parad)

        #Shift
        shiftx = (rotx - rotrefx)/(np.cos(decrad)*3600.) + center[0]
        shifty = (roty - rotrefy)/3600. + center[1]

        c.center = (shiftx, shifty)
        c.radius *= scale

    if recurse:
        return patches, transform_patches(refpatches, pa=pa, center=center,
                                          reffiber=reffiber, scale=scale,
                                          refpatches=None)[0]
    else:
        return patches, False

def wcs2pix(patches, header):
    """
    Given a pyfits header object, convert the centers of the GradPak patches
    from WCS coordinates (decimal degrees) to pixel coordinates on the
    corresponding FITs file. Radial scaling is done as part of
    transform_patches().
    """
    header_wcs = pywcs.WCS(header)

    for c in patches[:,1]:
        c.center= tuple(header_wcs.wcs_sky2pix([c.center],0)[0])

    return patches

def prep_axis(fitsfile = None, invert = True, sky = False, imrot = False,
              wcsax = True, figure = None, figsize = (8,8), geometry = (1,1,1)):
    """Create a pyplot Axes object and get it ready to receive some GradPack
    patches. This includes the creation of a colorbar and setting reasonable
    limits on the plot. It is also possible to provide a FITS image that will
    be displayed on the axes with WCS axis labels.

    Parameters
    ----------
    
    fitsfile : str
        The name of a fits file to plot. This 
        file MUST have a well formed, useful header.

    invert : bool
        If a fitsfile is provided, do you want to invert the colormap?

    sky : bool
        Do you want to show sky fibers in the plot?

    imrot : float
        The angle (in degrees, E of N) to rotate the fitsfile by

    wcsax : bool
        Display axis labels in WCS coordinates?

    figsize : list or tuple
        Width x height of figure, in inches

    figure : matplotlib.pyplot.Figure instance 
        If not None, the axis object created by will be placed into
        this figure with the geometry specified by **geometry**.

    geometry : int or tup
        A pyplot-style geometry argument for placing the axis in its
        parent figure. Can be a 3 digit integer (e.g., 111) or a tuple
        (e.g., (1,1,1)

    Returns
    -------

    ax : matplotlib.axes.Axes
        The requested axis, all setup and ready for plotting

    hdu : None or pyfits.PrimaryHDU
        The HDU used to create the WCS axis

    """

    if fitsfile:
        try:
            import pywcs
            global pywcs
        except ImportError:
            raise Exception("Could not find pywcs! Without it you cannot use the fitsfile option :(")
        
        hdu = pyfits.open(fitsfile)[0]
        imdata = hdu.data
        if imrot:
            imdata = spndi.rotate(imdata,-1*imrot,reshape=False)
            hdu.header.update('CROTA2',imrot)
        if wcsax:
            try:
                import pywcsgrid2 as wcsgrid
                global wcsgrid
            except ImportError:
                raise Exception("Could not fine pywcsgrid2! Withoutit you cannot use the wcsax option :(")
            axistype = (wcsgrid.Axes, dict(header=hdu.header))
        else:
            axistype = None
    else:
        hdu = None
        axistype = None

    if not figure:
        fig = plt.figure(figsize=figsize)

    grid = ImageGrid(figure, geometry,
                     nrows_ncols = (1,1),
                     cbar_mode = 'each',
                     cbar_location = 'top',
                     cbar_pad = '1%',
                     axes_class = axistype)
    ax = grid[0]
    
    if fitsfile:
        if invert:
            imdata = -1*(imdata - np.max(imdata))

        if wcsax:
            ax.set_display_coord_system('fk5')
                        
            #labdict = {'nbins':10}
            ax.set_ticklabel_type('dms','hms')
                                  # ,center_pixel=centpx,
                                  # labtyp1_kwargs=labdict,
                                  # labtyp2_kwargs=labdict)
            ax.add_compass(1)

        ax.imshow(imdata,
                  cmap = plt.get_cmap('gray'),
                  origin = 'lower', aspect = 'equal')
        
    else:
        ax.set_xlabel('arcsec')
        ax.set_ylabel('arcsec')
        if sky:
            ax.set_xlim(59,-59)
            ax.set_ylim(-2,116)
        else:
            ax.set_ylim(-2,60)
            ax.set_xlim(31,-31)
        
    return ax, hdu

def prep_patches(values,
                 binheader = None, plotbins = False,
                 hdu = None, pa = 0, center = [0,0], reffiber = 105,
                 sky = False, exclude = []):
    """
    Generate GradPak patches and prepare them for plotting. This function is
    used to transform the patches and remove any unwanted fibers (sky or user
    defined).

    Parameters
    ----------
    values : numpy.ndarray
        Should be length 109 and contain the data value of each 
        fiber, in order.
        
    binheader : pyfits.header.Header 
        A header (probably produced by :mod:`GradPak_bin`)
        containing information on what fibers went into what
        apertures

    plotbins : bool
        Somewhat strangely named; if True plot the raw GradPak
        fibers, regardless of aperture assignment (the actual
        fiber values will still be taken on an
        aperture-by-aperture basis).
            
    hdu : pyfits.PrimaryHDU
        Contains a header with a WCS solution that is used to 
        transform the patches.
        
    pa, center, reffiber : float
        Parameters that describe how to transform the
        patches. See :func:`transform_patches` for more info.

    sky : bool
        If True, plot sky fibers

    exclude : list
        List of fiber numbers to be excluded from the plot. If a
        **binheader** is provided this list should contain
        aperture numbers, not fiber numbers.

    Returns
    -------
    patches : list of matplotlib.patches.Circle
        Transformed and culled list of pyplot Patch objects

    pvals : numpy.ndarray
        Culled list of values so that pval[i] is the value of patch[i]

    refcenter : list of float 
        The pixel center of the reference fiber
    """
    if binheader is None or plotbins is False:
        patches = GradPak_patches()
        refpatches = None
    else:
        patches = get_binned_patches(binheader)
        refpatches = GradPak_patches()

    skyidx = [0,1,17,18,19,30,31,42,43,52,53,61,62,70,78,86,87,94,101,108]
    if hdu:
        scale = 2./((np.abs(hdu.header['CDELT1']) + \
                     np.abs(hdu.header['CDELT2']))*
                    3600.) #px/arcsec
        patches, refpatches = transform_patches(patches,
                                                reffiber = reffiber,
                                                pa = pa,
                                                center = center,
                                                scale = scale,
                                                refpatches = refpatches)
        patches = wcs2pix(patches, hdu.header) # now in px
        if refpatches is False:
            refpatches = patches
        else:
            refpatches = wcs2pix(refpatches, hdu.header)

    if binheader is None or plotbins is False:
        refpatches = patches
    refcenter = refpatches[reffiber - 1,1].center # in px
    
    if binheader is not None and plotbins is False:
        values = fill_fibers_with_bins(binheader, values)
        exclude = exclude_bins(binheader, exclude)

    if not sky and (binheader is None or plotbins is False):
        exclude = np.r_[skyidx,np.array(exclude)-1] #-1 needed because fiber
                                                    #numbers start at 1
    else:
        exclude = np.array(exclude) - 1

    exclude = np.array(exclude)
    exclude = np.unique(exclude) #in case the user specified some sky fibers
    patches = np.delete(patches, exclude, axis=0)
    values = np.delete(values, exclude)
    patches = patches[values == values]
    pval = values[values == values]

    return patches, pval, refcenter

def plot(values, binheader = None, plotbins = False,
         ax = None, figsize = (8,8), figure = None, geometry = (1,1,1), 
         fitsfile = None, imrot = False, wcsax = True, invert=True,
         pa = 0, center = [0,0], reffiber = 105, alpha=1.0,
         clabel = '', cmap = 'gnuplot2', minval = None, maxval = None,
         labelfibers = True, sky = False, exclude = []):
    """Generate a spatial plot of the GradPack IFU fibers with each fiber colored
    based on user-supplied values. This is one of the main user-level
    functions in this module, and returns an Axes object for integration into
    whatever higher-level plotting the user is doing.

    It is of the utmost importance that the values input variable is of length
    109 and is ordered by fiber number.

    A quick example of how to use this function:

    >>> GradPak_plot.plot(np.arange(109)).figure.show()

    This will show you a simple plot with the fibers colored by their fiber
    number and presented on a relative arcsec scale. More advanced usage can
    be achieved with the following options:

    Parameters
    ----------

    values : numpy.ndarray
        The data values associated with each fiber. The length of this array 
        should be equal to either 109 or the number of bins specified in **binheader**. 

    binheader : pyfits.header.Header
        A FITS header that contains keywords BINFXXX and BINPXXX. This header will be used
        to group fibers together during plotting.

    plotbins : bool
        An extremely poorly named parameter. If True, plot every *fiber* (not aperture) with
        a value given by the value of its aperture assignment as specified by **binheader**.
        If False the apertures specified by **binheader** will be plotted as contiguous regions.
        
    alpha : float
        The desired opacity of the plotted patches.
        
    ax : matplotlib.axes.Axes
        If supplied, the GradPack patches will be plotted
        in this axis. This is very useful for plotting multiple pointings on
        the same plot. Setting this option causes **fitsfile**, **imrot**, **invert**,
        and **wcsax** to be ignored.

    figure : matplotlib.pyplot.Figure instance
        If not None, the axis object created by :func:`plot` will be
        placed into this figure with the geometry specified by
        **geometry**.

    figsize : tup
        The size of the figure, in inches. Passed directly to
        plt.figure(). This option is ignored if **figure** is not None
    
    geometry : int or tup
        A pyplot-style geometry argument for placing the axis in its
        parent figure. Can be a 3 digit integer (e.g., 111) or a tuple
        (e.g., (1,1,1)

    fitsfile : str - 
        The name of a FITS image to draw on the plot. The
        FITS header must contain WCS parameters in the CDELT, CRVAL, CRPIX
        format.

    imrot : float
        Rotation of fits image in relation to the axes. This
        is useful for, e.g., aligning a galaxy's major axis along the x
        axis. This option is ignored if **fitsfile** = None or **ax** != None

    wcsax : bool
        If True the axis labels will be in Fk5 WCS
        coordinates. This option is ignored if **fitsfile** = None or **ax** != None

    invert : bool
        If True, the colormap of the fits image will be
        inverted. This option is ignored if **fitsfile** = None or **ax** != None

    pa : float
        Position angle of GradPak in decimal degrees. This
        angle is measured East of North and should be whatever you told the
        telescope operator.

    center : tup or list
        Length 2 tuple or list containing the coordinates of the GradPak array. 
        The units should be decimal Ra and Dec. This is probably the 
        coordinates you had listed in your cache.

    reffiber : int
        The IFU fiber placed at the coordinate given in
        center. Default is fiber 105.

    clabel : str
        The label of the colorbar. This is typically a
        description of the values being plotted.

    cmap : str
        The name of a matplotlib colormap that will be applied
        to the data values. This is passed directly to plt.get_cmap()
        Don't use 'jet'.
        
    minval/maxval : float
        The lower and upper limits of the colorbar,
        respectively. These are passed directly to
        matplotlib.colors.Normalize()

    labelfibers : bool
        If True, each fiber will be labeled with its
        fiber number.
    
    sky : bool
        If True, sky fibers will be plotted and the axes limits
        expanded to view the sky fibers.

    exclude : list
        A list of fiber numbers to be excluded from
        plotting. These patches are simply deleted. 

    Returns
    -------

    ax : matplotlib.axes.Axes
        The Axes containing all the plotting requested.

    """
    tmpax, hdu = prep_axis(fitsfile = fitsfile, 
                           invert = invert, 
                           sky = sky, 
                           imrot = imrot, 
                           wcsax = wcsax,
                           figsize = figsize,
                           geometry = geometry,
                           figure = figure)

    if not ax:
        ax = tmpax
  
    patches, pval, refcenter = prep_patches(values, binheader=binheader,
                                            plotbins = plotbins,
                                            hdu = hdu, pa = pa, 
                                            center = center,
                                            reffiber = reffiber,
                                            sky = sky, exclude = exclude)

    if hdu is not None:
        xdelt = 1.5/(60. * hdu.header['CDELT1'])
        ydelt = 1.5/(60. * hdu.header['CDELT2'])
        ax.set_xlim(refcenter[0] + xdelt, refcenter[0] - xdelt)
        ax.set_ylim(refcenter[1] - ydelt, refcenter[1] + ydelt)
    else:
        xmin, xmax = ax.get_xlim()
        ax.set_xlim(xmax, xmin)

    if labelfibers:
        for c in patches:
            ax.text(c[1].center[0],c[1].center[1],c[0],fontsize=7,
                    ha='center',va='center')
    
    if minval is None:
        minval = pval.min()
    if maxval is None:
        maxval = pval.max()

    collection = PatchCollection(patches[:,1],
                                 cmap=plt.get_cmap(cmap),
                                 norm=matplotlib.colors.Normalize(
                                     vmin=minval,vmax=maxval),
                                 edgecolor = 'none',
                                 alpha=alpha)
    collection.set_array(pval)
    ax.add_collection(collection)

    cbar = plt.colorbar(collection,cax=ax.cax,orientation='horizontal')
    cbar.set_label(clabel)
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')

    if not plotbins and binheader is not None:
        boxes, bval = get_bin_boxes(binheader, patches, pval)
        boxColl = PatchCollection(boxes, 
                                  cmap=plt.get_cmap(cmap),
                                  norm=matplotlib.colors.Normalize(
                                      vmin=minval,vmax=maxval),
                                  edgecolor = 'none',
                                  alpha=alpha)
        boxColl.set_array(bval)
        ax.add_collection(boxColl)

    return ax

def plot_img(values, 
             ax = None, figsize = (8,8),
             fitsfile = None, imrot = False, invert = True,
             pa=0, center=[0,0], reffiber = 105,
             clabel='', cmap='gnuplot2', minval = None, maxval = None,
             numpoints=500, method='nearest',
             sky = False, exclude=[]):
    """
    Generate an interplotaed image of the GradPak IFU using user supplied
    values. This is one of the main user-level functions in this module, and
    returns an Axes object for integration into whatever higher-level plotting
    the user is doing.

    It is of the utmost importance that the values input variable is of length
    109 and is ordered by fiber number.

    A quick example of how to use this function:

    >>> GradPak_plot.plot_img(np.arange(109)).figure.show()

    This will show you a simple image of GradPack colored by fiber number and
    presented on a relative arcsec scale. More advanced usage can be achieved
    with the following options:

    Input Options:

        o ax (pyplot.Axes) - If supplied, the GradPack image will be plotted
          in this axis. This is very useful for plotting multiple pointings on
          the same plot. Setting this option causes fitsfile, imrot, invert,
          and wcsax to be ignored.

        o figsize (tup) - The size of the figure, in inches. Passed directly
          to plt.figure()
    
        o fitsfile (str) - The name of a FITS image to draw on the plot. The
          FITS header must contain WCS parameters in the CDELT, CRVAL, CRPIX
          format.

        o imrot (float) - Rotation of fits image in relation to the axes. This
          is useful for, e.g., aligning a galaxy's major axis along the x
          axis. This option is ignored if fitsfile = None or ax != None

        o wcsax (bool) - If True the axis labels will be in Fk5 WCS
          coordinates. This option is ignored if fitsfile = None or ax != None

        o invert (bool) - If True, the colormap of the fits image will be
          inverted. This option is ignored if fitsfile = None or ax != None

        o pa (float) - Position angle of GradPak in decimal degrees. This
          angle is measured East of North and should be whatever you told the
          telescope operator.

        o center (tup or list) - Length 2 tuple or list containing the
          coordinates of the GradPak array. The units should be decimal Ra and
          Dec. This is probably the coordinates you had listed in your cache.

        o reffiber (int) - The IFU fiber placed at the coordinate given in
          center. Default is the lower left fiber (viewed on Wifoe), which is
          fiber 105.

        o clabel (str) - The label of the colorbar. This is typically a
          description of the values being plotted.

        o cmap (str) - The name of a matplotlib colormap that will be applied
          to the data values. This is passed directly to plt.get_cmap()

        o minval/maxval (float) - The lower and upper limits of the colorbar,
          respectively. These are passed directly to
          matplotlib.colors.Normalize()
    
        o numpoints (int) - The number of points to use when interpolating the
          GradPak IFU. This is actually the number of points to interpolate
          over in each dimension, so the total number of points will be
          numpoints**2

        o method (str) - The interpolation method used. This is passed
          directly to spi.griddata(). The available options are:

            'nearest', 'linear', 'cubic'

        o sky (bool) - If True, sky fibers will be plotted and the axes limits
          expanded to view the sky fibers.

        o exclude (list) - A list of fiber numbers to be excluded from
          plotting. These patches are simply deleted

    Output:
        ax (pyplot.Axes) - The Axes containing all the plotting requested.

    """
    if not ax:
        ax, hdu = prep_axis(fitsfile, invert, sky, imrot, figsize)
        
    patches, pval, refcenter = prep_patches(values,
                                            hdu = hdu, pa = pa, 
                                            center = center,
                                            reffiber = reffiber,
                                            sky = sky, exclude = exclude)
    if hdu is not None:
        xdelt = 2./(60. * hdu.header['CDELT1'])
        ydelt = 2./(60. * hdu.header['CDELT2'])
        ax.set_xlim(refcenter[0] + xdelt, refcenter[0] - xdelt)
        ax.set_ylim(refcenter[1] - ydelt, refcenter[1] + ydelt)
    else:
        xmin, xmax = ax.get_xlim()
        ax.set_xlim(xmax, xmin)

    x = np.array([c.center[0] for c in patches[:,1]])
    y = np.array([c.center[1] for c in patches[:,1]])

    xi = np.linspace(x.min(),x.max(),numpoints)
    yi = np.linspace(y.min(),y.max(),numpoints)

    vi = spi.griddata((x,y),
                      pval,
                      (xi[None,:],yi[:,None]),
                      method=method,
                      fill_value=np.inf)
    if minval is None:
        minval = pval.min()
    if maxval is None:
        maxval = pval.max()

    im = ax.imshow(vi, cmap=cmap, origin='lower', 
                   extent=(xi.min(),xi.max(),yi.min(),yi.max()),
                   vmin=minval,vmax=maxval)
    cbar = ax.cax.colorbar(im)
    cbar.set_label_text(clabel)

    return ax

def plot_rows(values, binheader = None,
              ylabel='', label='',
              ax = None, fullout = False,
              weights=None, err=False,
              kpc_scale=None, zcorr=0,
              **plot_kwargs):
    """
    Bin values by GradPak row and produce a plot of the results.  Each row's
    value is the weighted average of the fibers in that row. This is one of
    the main user-level functions in this module, and returns an Axes object
    for integration into whatever higher-level plotting the user is doing.

    It is of the utmost importance that the values input variable is of length
    109 and is ordered by fiber number.

    A quick example of how to use this function:

    >>> GradPak_plot.plot_img(np.arange(109)).figure.show()

    This will show you a simple image of GradPack colored by fiber number and
    presented on a relative arcsec scale. More advanced usage can be achieved
    with the following options:

    Input Options:

        o ylabel (str) - Y axis labels. If ax != None this option is ignored.
    
        o label (str) - Label applied to plotted lines. Useful for any future
          calls to ax.legend()

        o ax (pyplot.Axes) - If supplied, the GradPack patches will be plotted
          in this axis. This is very useful for plotting multiple pointings on
          the same plot. Setting this option causes fitsfile, imrot, invert,
          and wcsax to be ignored.

        o weights (length 109 numpy array) - Weights to be applied to
          individual fibers when averaging rows together. If None then all
          fibers are assumed to have equal weight.

        o err (bool) - If True the standard error is used for error bars. If
          False then the weighted standard deviation is used for plot error
          bars.

        o kpc_scale (float) - Kpc/arcsec scale used to convert height (x axis)
          to something meaningful.

        o fullout (bool) - If True the binned values and associated errors are
          returned along with the Axes.

        o **plot_kwargs (dict) - Plotting keywords passed directly to
          ax.errorbar

    Outputs:
    
        o ax (pyplot.Axes) - The Axes containing all the plotting requested.
    
        o abcissa (ndarray) - The height of each binned row. If kpc_scale was
          set then the units of abcissa are kpc, otherwise they are arcsec.

        o binned_vals (ndarray) - The binned (weighted avg.) values of each
          GradPak row.

        o binned_errs (ndarray) - The standard error of each row

        o binned_stds (ndarray) - The weighted standard deviation of each row

    """
    if binheader:
        y_values = np.array([c.center[1] for c in get_binned_patches(binheader)[:,1]])
    else:
        y_values = np.array([c.center[1] for c in GradPak_patches()[:,1]])
    row_pos = np.unique(y_values)
    binned_vals = np.array([])
    binned_errs = np.array([])
    binned_stds = np.array([])
    abcissa = np.array([])
    if weights is None:
        weights = np.ones(y_values.size)

    for i, row in enumerate(row_pos):
        if row < 80: #Don't consider sky fibers
            idx = np.where(y_values == row)[0]
            abcissa = np.append(abcissa, row)
            mean = np.sum(values[idx]*weights[idx]/np.sum(weights[idx]))
            std = np.sqrt(
                np.sum(weights[idx]*(values[idx] - mean)**2) /
                ((idx.size - 1.)/(idx.size) * np.sum(weights[idx])))
            std /= np.sqrt(idx.size)
            err = mean*np.sqrt(1./np.sum(weights[idx]**2))
            binned_vals = np.append(binned_vals,mean)
            binned_errs = np.append(binned_errs,err)
            binned_stds = np.append(binned_stds,std)

    if kpc_scale is not None:
        abcissa *= kpc_scale
        abcissa += zcorr
        xlabel = 'Height [kpc]'
    else:
        xlabel = 'Height [arcsec]'

    if ax is None:
        ax = plt.figure(figsize=(8,8)).add_subplot(111)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    if err:
        ax.errorbar(abcissa, binned_vals, 
                    yerr = binned_errs,label = label, **plot_kwargs)
    else:
        ax.errorbar(abcissa, binned_vals, 
                    yerr = binned_stds,label = label, **plot_kwargs)
        
    if fullout:
        return ax, abcissa, binned_vals, binned_errs, binned_stds
    else:
        return ax
    

def format_tpl(tpl):
    """Take in a ds9 .tpl file and print out corresponding pyplot Patch
    definitions.

    WARNING: It is almost certainly the case that ds9 template files are NOT
    listed in fiber order. It is therefore NOT SAFE to assume that the output
    of this function can be used to produce something like
    GradPak_patches(). Manual confirmation of fiber order is crucial.
    """
    
    with open(tpl,'r') as f:
        lines = f.readlines()
    
    for line in lines:
        if line[0:6] == 'circle':
            data = line.split('(')[1].split(')')[0].replace('"','')
            x,y,r = [float(i) for i in data.split(',')]
            print 'Circle(({:6.4f}, {:6.4f}), radius={:6.4f}),'.format(
                x*-3600.,y*3600,r)

    return
