#! /usr/bin/python

"""
************
GradPak_scat
************

A quick and non-so-dirty scattered light correction designed for GradPak
data. You'll need this because IRAF's scatter light routines cannot deal with
the very small gaps between the fibers.

The algorithm is simple. For each wavelength channel the "scattered light" is
just a linear interpolation of the average of the dark regions beyond the
edges of the slit. Essentially it finds the light scattered by the edge fibers
and interpolates it across the whole slit. It's not perfect, but it's better
than nothing and it works.

The calling syntax is:
::

    > python GradPak_scat.py input.lst output.lst

where input.lst and output.lst contain the names of the input and output fits
files, one on each line. If you've ever used @lists in IRAF you will
understand what is required here.

.. warning:: The pixel locations of the off-slit dark regions is currently
    hardcoded based on a 2x2 binning scheme. If you didn't use 2x2 binning
    this will not work as expected.

Functions
---------

"""

import sys
import pyfits
import numpy as np


def do_single(input_image, output_image):
    """Perform a GradPak scattered light correction on a single image.
    
    First compute the mean "gutter" (the regions beyond the edges of the slit)
    spectrum, then interpolated between the two gutters at each
    wavelength. Finally, subtract this interpolated image from the input and
    write the results to the output file.

    Parameters
    ----------

    input_image : str
        Name of the input FITS image

    output_image : str
        Name of the file to contain the scattered light corrected FITS image

    Returns
    -------

    None

    Notes
    -----

    Currently hardwired for 2x2 binning. Will likely work in any case
    where binning in the spatial dimension is 2, but this has not been tested.

    """
    h = pyfits.open(input_image)[0]
    d = h.data
    sr = np.zeros(d.shape,dtype=d.dtype)
    sr[:,0] = np.mean(d[:,0:70], axis=1)
    sr[:,-1] = np.mean(d[:,1260:], axis=1)
    
    x = np.arange(sr.shape[1])
    xp = np.array([x[0],x[-1]])
    fp = np.array([sr[:,0], sr[:,-1]])
    
    l = [np.interp(x,xp,fp[:,i]).astype(d.dtype) for i in range(fp.shape[1])]
    scat = np.vstack(l)
    
    pyfits.PrimaryHDU(d - scat, h.header).writeto(output_image)

    return

def do_multi(inputs, outputs):
    """Read input and output files and run :meth:`GradPak_scat.do_single` on each pair

    Very simple; just read the lines and run the main function. The only
    little hiccup is that readlines() keeps the newline characters around in the
    names so we need to nuke those.

    Parameters
    ----------

    inputs : str
        Name of the file containing input image names, one per line

    outputs : str
        Name of the file containing output image names, one per line

    Returns
    -------
 
    None

    """
    with open(inputs,'r') as fi:
        with open(outputs,'r') as fo:
            input_list = fi.readlines()
            output_list = fo.readlines()
            for i, o in zip(input_list, output_list):
                i = i.replace('\n','')
                o = o.replace('\n','')
                print i, o
                do_single(i,o)

    return

if __name__ == '__main__':
    do_multi(sys.argv[1],sys.argv[2])
