#! /usr/bin/python
#
"""
***********
GradPak_bin
***********

This module provides useful functions for binning GradPak fibers together to
achieve higher Signal to Noise (SNR). Given a desired SNR the algorithm will
bin fibers together until either the SNR or the end of a row is reached. The
current main limitation with this algorithm is that it will not smartly add a
single fiber at the end of a row to the preceeding bin even if that single
fiber has a low SNR.

There are many ways to compute the SNR, this module does it simply:

SNR = <signal/noise>

where signal are the values in the datafile and noise are the values in the errfile.

In addition to generating binned spectra this module provides a method for
computing the location of each bin, in arcsec or kpc, in a coordinate system
defined by the user (usually relative to some galaxy's center).

Useful Example
--------------

For my NGC 891 work I use::

 >>> GPB.bin('NGC_891_P3.ms_rfs_lin.fits','NGC_891_P3.mse_rfs_lin.fits',30,'NGC_891_P3_bin30',exclude=[92,93])
 >>> GPB.create_locations('NGC_891_P3_bin30.ms.fits',ifucenter=[35.636325,42.37935])

Functions
---------
"""

import time
import numpy as np
import pyfits
import GradPak_plot as GPP
import matplotlib.pyplot as plt

def bin(datafile, errfile, SNR, outputfile, waverange=None, exclude=[], 
        logfile=None, logfits=None):
    """Bin GradPak fibers along a row until the desired SNR is achieved.

    The way fibers are grouped together is somewhat simplistic, but the
    underlying mathematics is correct. 

    The first bin starts at the smallest non-sky, non-excluded fiber (usually
    fiber #3) and just keeps adding the next fiber until the desire SNR is
    reached. Fibers are not binned across rows, which ensures that each bin is
    made up only of fibers of one size, but also limits the amount of binning
    that can be done. The major limitation of this method is that the ends of
    rows can often be left with a single, low SNR fiber. For example, if a bin
    is constructed containing fibers 14-16 and fiber 17 (at the end of the
    row) is below the threshold it will not be added to the previous bin and
    instead will make up a single bin with low a low SNR. This should be easy
    to fix but in practice it doesn't cause that much of an issue so I just
    haven't done it.

    The spectrum for each bin is an average of all the contributing fibers,
    weighted by each individual's SNR^2. The corresponding error bin is a sum
    in quadrature, weighted in the same way.

    The history of each bin is recorded in the output FITS header using two keywords:

     **BINXXXF** The fibers that went into bin XXX

     **BINXXXP** The position (in arcsec) of the center of bin XXX. This is
       the *unweighted* average of the centers of the contributing fibers.

    Parameters
    ----------

    datafile : str
        The name of a multispec FITS file that you want to bin

    errfile : str
        The name of a multispec FITS file containing error vectors corresponding to the datafile

    SNR : float
        The minimum SNR in each bin

    outputfile : str 
        The base name of the binned spectral and error arrays. The results
        will be called OUTPUTFILE.ms.fits and OUTPUTFILE.me.fits.

    waverange : list (default: None)
        A 2 element list containing the minimum and maximum wavelengths
        to consider *for the SNR calculations*. Regardless of this parameter
        the final output will span the same wavelength range as the input
        files.
    
    exclude : list (default: [])
        A list containing fibers to exlcude from binning. Sky fibers are
        automatically excluded and fiber numbers start at 1.

    logfile : str
        Name of a log file that will contain, for each aperture, a
        list of the individual fibers and the total S/N.

    logfits : str
        Name of a FITS file that will contain a separate HDU for each
        aperture. Each HDU will contain an array with two dimensions
        that records the fiber numbers and associated weights for that
        aperture.
        
    Returns
    -------
    
    finalf : numpy.ndarray
        Array with shape NBINS x NWAVE containing the binned data

    finatle : numpy.ndarray
        Array with same shape as finalf containing the errors on the binned data

    fibdict : dict
        Dictionary of questionable usefulness that contains information about
        which fibers went into each bin. The keys are ROW_BIN and the entries
        are a list of the fibers that went into that bin (e.g., '1_2' for the
        second bin in the first row).

    """
    hdu = pyfits.open(datafile)[0]
    data = hdu.data
    err = pyfits.open(errfile)[0].data

    if waverange is not None:
        wave = (np.arange(data.shape[1]) + hdu.header['CRPIX1'] - 1)\
               *hdu.header['CDELT1'] + hdu.header['CRVAL1']
        waveidx = np.where((wave >= waverange[0]) & (wave <= waverange[1]))[0]
    else:
        waveidx = None

    data *= 1e17
    err *= 1e17

    y_values = np.array([c.center[1] for c in GPP.GradPak_patches()[:,1]])
    x_values = np.array([c.center[0] for c in GPP.GradPak_patches()[:,1]])
    fibnums = np.arange(109) + 1
    row_pos = np.unique(y_values)

    finalf = np.zeros(data.shape[1])
    finale = np.zeros(data.shape[1])

    fibdict = {}
    binnum = 1

    if logfile is not None:
        lf = open(logfile,'w')
        
    if logfits is not None:
        log_HDUs = [] #This will hold [fiber_list, weight_list] for each aperture
        
    for i in range(row_pos.size):

        if row_pos[i] > 80:
            continue

        idx = np.where(y_values == row_pos[i])[0]
        b = 0
        n = 0

        while fibnums[idx[n]] in exclude:
            print 'Skipping fiber {}'.format(fibnums[idx[n]])
            n += 1

        while n < len(idx):
            try:
                while fibnums[idx[n]] in exclude:
                    print 'Skipping fiber {}'.format(fibnums[idx[n]])
                    n += 1
            except IndexError: #The rest of the fibers in the row are excluded
                break

            fstack = data[idx[n]][None,:]
            estack = err[idx[n]][None,:]
            tmp = compute_SN(data[idx[n]], err[idx[n]], waveidx)
            snstack = np.array([[tmp]])
            fibers = [fibnums[idx[n]]]
            xpos = [x_values[idx[n]]]
            ypos = [y_values[idx[n]]]

            while tmp < SNR:
                n += 1
                print 'fibers: {}, SNR: {}'.format(fibers, tmp)
                if n > len(idx) - 1:
                    print "WARNING, SN threshold not met in row {}, bin {}".\
                        format(i,b)
                    break

                if fibnums[idx[n]] in exclude:
                    print 'Skipping fiber {}'.format(fibnums[idx[n]])
                    continue
                    
                fstack = np.vstack((fstack, data[idx[n]]))
                estack = np.vstack((estack, err[idx[n]]))
                snstack = np.vstack((snstack, compute_SN(data[idx[n]], 
                                                         err[idx[n]], 
                                                         waveidx)))
                tmpf, tmpe = create_bin(fstack, estack, snstack)
                tmp = compute_SN(tmpf, tmpe, waveidx)
                fibers.append(fibnums[idx[n]])
                xpos.append(x_values[idx[n]])
                ypos.append(y_values[idx[n]])

            binf, bine = create_bin(fstack, estack, snstack)
            binsn = compute_SN(binf, bine, waveidx)
            print 'binned aperture {}: {}, SNR: {}'.format(binnum,fibers, binsn)
            if logfile is not None:
                lf.write('binned aperture {}: {}, SNR: {}\n'.format(binnum,fibers, binsn))
            if logfits is not None:
                log_HDUs.append(pyfits.ImageHDU(np.vstack([fibers, snstack.T[0]])))

            bin_x_pos = np.mean(xpos)
            bin_y_pos = np.mean(ypos)
            fibstr = [str(i) for i in fibers]
            hdu.header.update('BIN{:03}F'.format(binnum),' '.join(fibstr))
            hdu.header.update('BIN{:03}P'.format(binnum),' '.\
                              join([str(bin_x_pos),str(bin_y_pos)]))

            finalf = np.vstack((finalf,binf))
            finale = np.vstack((finale,bine))
            fibdict['{}_{}'.format(i,b)] = fibers
            b += 1
            n += 1
            binnum += 1

    finalf = finalf[1:]/1e17
    finale = finale[1:]/1e17
    pyfits.PrimaryHDU(finalf, hdu.header).\
        writeto('{}.ms.fits'.format(outputfile),clobber=True)
    pyfits.PrimaryHDU(finale, hdu.header).\
        writeto('{}.me.fits'.format(outputfile),clobber=True)

    if logfits is not None:
        lP = pyfits.PrimaryHDU()
        lP.header.update('HDUDIM','APERTURE')
        lP.header.update('AXIS1','FIBERS')
        lP.header.update('AXIS2','0: FIBER NUMBER 1: WEIGHT')
        pyfits.HDUList([lP] + log_HDUs).writeto(logfits, clobber=True)
        
    return finalf, finale, fibdict

def create_bin(fstack, estack, snstack):
    """Weight and combine data and errors into new bins.

    The weights are simply SNR^2.

    Parameters
    ----------

    fstack : numpy.ndarray
        A NxM array where of data values N is the number of fibers to add to
        the bin and M is the number of wavelength channels. Generally the 1st
        element of this array is the current bin, not an individual fiber.

    estack : numpy.ndarray
        The errors corresponding to fstack. Should be the same shape as fstack

    snstack : numpy.ndarray
        Array containing the SNR of each of the N entries in fstack

    Returns
    -------
    
    newbin : numpy.ndarray
        Array of length M containing the weighted average of fstack
    
    newerr : numpy.ndarray
        Array of length M containing the errors on newbin. Errors are combined
        in weighted quadrature.

    """
    sumW = np.sum(snstack**2)

    newbin = np.sum(fstack*snstack**2,axis=0)/sumW
    newerr = np.sqrt(np.sum((estack*snstack**2/sumW)**2,axis=0))

    return newbin, newerr

def compute_SN(signal, noise, idx=None):
    """Given signal and noise vectors, compute the SNR (possibly over a specific range).

    Parameters
    ----------

    signal : numpy.ndarray
        A 1D vector of data values

    noise : numpy.ndarray
        A 1D vector of error values

    idx : numpy.ndarray (default=None)
        A 1D vector containing the index values of a sub region used to
        compute the SNR

    Returns
    -------

    SNR : float
        The signal to noise. Computed simply as <signal/noise>

    """
    zidx = np.where(noise[idx] != 0)

    return np.mean(signal[idx][zidx]/(noise[idx][zidx]))

def create_locations(binfile, galcenter=[35.637962,42.347629], 
                     ifucenter=[35.637962,42.347629], reffiber=105,
                     galpa=293.3, ifupa=295.787, kpc_scale=0.0485):
    """Given a binned data file, creat a list of bin locations, relative to a galactic center.

    The function is designed to give physical meaning to the location of each
    bin created by :func:`bin`. For each bin found in the FITS header it
    computes physical coordinates in both arcsec and kpc. This method was
    written with a 2D projection of cylindrical coordinates in mind (r,z), but
    the user is free to fully define the location and meaning of the
    coordinate system through keyword options.

    The output is a text file containig the location of each bin in the
    desired coordinate system.

    The current defaults are for NGC 891 and WIYN proposal 14B-0456. You
    should probably change them.

    Parameters
    ----------
    
    binfile : str
        The name of a FITS file produced by :func:`bin`
    
    galcenter : list, optional
        A list containing the [RA,DEC] coordinates of the center of the
        coordinate system. The units are decimal degrees.

    ifucenter : list, optional 
        A list containing the center of the fiber specified in
        **reffiber**. Units are decimal degrees. This is usually the
        coordinates of your GradPak pointing.

    reffiber : int, optional
        The *fiber* number whoes coordinates are given in **ifucenter**. This
        is generally the fiber you used to align GradPak during observations.

    galpa : float, optional
        The position angle of the desired coordinate system. This angle
        defines rotation of the horizontal axis relative to the sky.
    
    ifupa : float, optional 
        The position angle of GradPak. This is defined as an absolute PA,
        relative to North, *not* to the defined coordinate system. This is
        generally the PA you entered during GradPak observations.
    
    kpc_scale : float, optional
        kpc/arcsec of the object in question

    Returns
    -------

    None :
        The result is file with the suffix "locations.dat" that contains
        information about the location in the user defined coordinates of each
        aperture (bin). Coordinates are given in arcsec and kpc. Note that the
        reported "size" is the size of the fibers that went into each bin, not
        the size of the bin itself. This is admittedly confusing.

    """

    hdu = pyfits.open(binfile)[0]
    numaps = hdu.data.shape[0]
    binhead = hdu.header
    
    patches = GPP.get_binned_patches(binhead)
    refpatches = GPP.GradPak_patches()
    
    patches, refpatches = GPP.transform_patches(patches,refpatches=refpatches,
                                                pa=ifupa, center=ifucenter, 
                                                reffiber=reffiber)

    decrad = galcenter[1]*2*np.pi/360.
    parad = galpa*2*np.pi/360.
    
    f = open('{}_locations.dat'.format(binfile.split('.ms.fits')[0]),'w')
    f.write("""# Generated on {}
# Inpute file: {}
#
""".format(time.asctime(),binfile))
    f.write('# {:4}{:>10}{:>10}{:>10}{:>10}{:>10}\n#\n'.format('Apnum',
                                                              'size (")',
                                                              'r (")',
                                                              'z (")',
                                                              'r (kpc)',
                                                              'z (kpc)'))

    for i, p in enumerate(patches[:,1]):
    
        fibers = binhead['BIN{:03}F'.format(i+1)]
        radius = refpatches[int(fibers.split(' ')[0]) - 1][1].get_radius()

        ra_diff = 3600*(galcenter[0] - p.center[0])*np.cos(decrad)
        dec_diff = 3600*(galcenter[1] - p.center[1])

        r_diff = ra_diff*np.cos(parad) - dec_diff*np.sin(parad)
        z_diff = -1*(ra_diff*np.sin(parad) + dec_diff*np.cos(parad))

        print i+1, p.center, ra_diff, dec_diff, r_diff, z_diff

        f.write(str('{:7n}'+5*'{:10.3f}'+'\n').format(i,
                                                             radius,
                                                             r_diff,
                                                             z_diff,
                                                             r_diff*kpc_scale,
                                                             z_diff*kpc_scale))
    f.close()
    return

def plot_test():

    import matplotlib.pyplot as plt
    for SN in [0, 20, 60]:

        flux, err, _ = bin('../NGC_891_P1_final.ms_rfsz_lin.fits',
                           '../NGC_891_P1_final.me_rfz_lin.fits',
                           SN, 'NGC_891_P1_bin{}'.format(SN), 
                           waverange=[5450,5550])

        ax = plt.figure().add_subplot(111)
        ax.plot(np.arange(flux.shape[1]), flux[3])
        ax.set_title('SN > {}'.format(SN))
        ax.set_ylim(-0.2e-15,0.8e-15)
        ax.figure.show()



