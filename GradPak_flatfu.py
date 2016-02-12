#! /usr/bin/python
#
# History:
#      v1 - A. Eigenbrot Nov. 2014
#      v1.1 - A. Eigenbrot Dec. 2014
#      v2.0 - A. Eigenbrot Aug. 2015
#               Corrected implimentation of fitflat+
#      v2.1 - A. Eigenbrot Oct. 2015
#               Added shutter lag correction
#
####################################################
"""
GradPak_flatfu
==============

This script combines up multiple flats so they can be used as a single input to dohydra. You will want to do this because an exposure that puts enough signal in the 200 micron fibers will push the larger fibers into the non-linear regime. The calling syntax is::
 
    > python GradPak_flatfu.py Flat1 Flat2... Flatn pivot1 pivot2... pivotn-1

Where the pivots define the aperture at which to cut the flat,
inclusive. For example, if the call is::

    > python GradPak_flatfu.py Flat1.fits Flat2.fits 70

Then the resulting master flat will have apertures 1 - 70 from
Flat1.fits and 71 - 109 from Flat2.fits.

By default, Flat1 is the flat used for the aperture extraction of the entire run, but this can be changed with the *-t* flag (see below).

All the flats are scaled by their exposure times and then cut up as specified by the user. The final, master flat lives only as a .ms multispec file. For this reason it is imperative that you set the "flat" parameter to EXACTLY "dFlat_master.fits". Nothing else will work. Even when you do this dohydra will complain that it can't find dFlat_master.fits. This warning is OK; dohydra will still do everything you want it to.

Options
-------

-t        Use the next flat as the flat for tracing apertures. **NOTE** the resulting aperture definitions will be used for all of your data, so make sure it's good.
-r file   Use the flat in *file* for fiber throughput determinations. Equivalent to *throughput=file* in **dohydra**
-nf       Do NOT divide the final flat by its average spectrum. Normally this division is done so that the spectral signature of the flat lamps is not divided into the data frames, but if you want to keep that spectrum this option is for you. Similar to *fitflat-* in **dohydra**.

Useful Example
--------------

I have found that often the lower exposure flat fields produce better aperture traces (because the large fibers are not saturated), and that typically it is a good idea to pivot around fiber 43 so that the 2'' and 3'' fibers use one flat and the 4'', 5'', and 6'' use another. With this in mind I almost always type::

    > python GradPak_flatfu.py dFlat_4s.fits -t dFlat_1s.fits 43

Functions
---------
"""

import os
import sys
import pyfits
import numpy as np

#Load the IRAF packages we'll need
try:
    current_dir = os.getcwd()
    if os.getlogin() == 'Arthur':
            os.chdir('/Users/Arthur/Ureka/iraf/local')
    from pyraf import iraf
    os.chdir(current_dir)
    iraf.imred(_doprint=0)
    iraf.hydra(_doprint=0)
except Exception as e:
    print "Failure: could not find pyraf/iraf"
    sys.exit(1)

if os.getlogin() == 'Arthur':
    APIDTABLE = '/Users/Arthur/Documents/School/MetaPak/gradpak_sizes.iraf'
else:
    APIDTABLE = '/usr/users/eigenbrot/research/Pak/gradpak_sizes.iraf'

def scale_images(hdulist):
    """Take in a list of fits HDUS and scale the data in all of them to the exposure time of the first HDU

    I think this function is no longer used, so I won't document it.
    """
    exptimes = [h.header['EXPTIME'] for h in hdulist]
    scales = exptimes[0]/np.array(exptimes)
    
    print 'Scaling flats...'
    for h, scale in zip(hdulist, scales):
        print '\t{} {}'.format(h.header['OBJECT'],scale)
        h.data *= scale
        h.header['EXPTIME'] *= scale

    return hdulist

def scale_spectra(imagelist):
    """Scale spectra based on an overlapping group of fibers.

    I think this function is no longer used, so I won't document it.
    """
    fiber1 = 20
    fiber2 = 43

    hdulist = [pyfits.open(image)[0] for image in imagelist]
    means = [np.mean(h.data[fiber1 - 1:fiber2 - 1,:]) for h in hdulist]
    scales = means[0]/np.array(means)

    print 'Scaling extracted flats...'
    outputnames = ['{}_scale.ms.fits'.format(image.split('.ms.fits')[0]) for image in imagelist]
    for h, scale, name in zip(hdulist, scales, outputnames):
        print '\t{} {}'.format(name, scale)
        h.data *= scale
        h.writeto(name,clobber=True)

    return outputnames

def shutter_correction(imagelist):
    """Apply a correction to account for spectral signatures caused by lagging shutter.

    The WIYN bench has linear shutter that has asymmetric opening and closing times. At exposure times less than ~7 seconds this asymmetry can introduce an artificial spectral signature at the few % level. Furthermore, the magnitude of this artificial spectrum varies with exposure time. This becomes a problem when trying to stitch together two flats taken at short, but not equal exposure times; each flat will have a different "shutter spectrum".

    To account for this we look at a set of fibers expected to have good signal in both flats and use them to normalize both the spectral shapes and total flux to a common value.
    
    Parameters
    ----------

    imagelist : list of str
        The names of IRAF multispec files to perform the correction on

    Returns
    -------
    outputnames : list of str
        The main result is to create and write new files with the shutter correction applied. These are giving the suffix "_shut.ms.fits" and their names are returned in a list.

    Notes
    -----

    The current default is to use the 3'' fibers as the overlapping set of fibers on which to base the correction. This might not be always true for all GradPak usage cases. Watch out.

    """

    fiber1 = 20
    fiber2 = 43

    hdulist = [pyfits.open(image)[0] for image in imagelist]
    avgs = [np.mean(h.data[fiber1 - 1:fiber2 - 1,:],axis=0) for h in hdulist]
    corrections = avgs[0]/np.array(avgs)

    pyfits.PrimaryHDU(corrections,hdulist[0].header).writeto('corrections.ms.fits',clobber=True)

    outputnames = ['{}_shut.ms.fits'.format(image.split('.ms.fits')[0]) for image in imagelist]
    print 'Correcting for shutter lag...'
    for h, corr, name in zip(hdulist, corrections, outputnames):
        print '\t{} {}'.format(name,np.mean(corr))
        h.data *= corr
        h.writeto(name,clobber=True)

    return outputnames

def setup_files(flatlist):
    """Prepare flats for dohydra run

    Take in a list of names of flat fields, scale the data by exposure
    time, and then write them out to appropriatly named files.
    Also create a dummy spectrum that will be "reduced" by dohydra.
    
    This function might not be used anymore.
    """
    hdulist = [pyfits.open(i)[0] for i in flatlist]
    hdulist = scale_images(hdulist)

    scalednames = []
    print 'Writing flats...'
    for h, flat in zip(hdulist, flatlist):
        newname = '{}_scale.fits'.format(flat.split('.fits')[0])
        print '\t{}'.format(newname), np.mean(h.data[:,3:70]/hdulist[0].data[:,3:70])
        h.writeto(newname,clobber=True)
        scalednames.append(newname)
        
    make_tmp(hdulist[0])

    return scalednames

def make_tmp(hdu):
    """Make the dummy spectrum to be reduced by **dohydra**. 

    This spectrum is passed as the 'object' spectrum to **dohydra** and only exists to allow the routine to run properly. It is set to all ones so that the affect of the flat scaling can be seen.

    Parameters
    ----------

    hdu : pyfits.PrimaryHDU object
        HDU that is characteristic of the flats and data used during reduction. This header will be copied to the temporary file and the object changed to 'FlatFu tmp'.

    Returns
    -------
    None
        The result is 'ffTmp.fits', which will passed as the object to internal **dohydra** runs.

    """
    tmpdata = np.ones(hdu.data.shape)
    tmphdu = pyfits.PrimaryHDU(tmpdata,hdu.header)
    tmphdu.header['OBJECT'] = 'FlatFu tmp'

    print 'Writing tmp spectrum'
    tmphdu.writeto('ffTmp.fits',clobber=True)

    return

def initial_run(scalednames,traceflat,throughput=''):
    """Use **dohydra** to extract .ms multispec files from the input flat fields. 
    
    This is done so that we can stitch them back together in aperture space rather than pixel space. To save time we don't bother with any sort of wavelength solution or other garbage. Essentially all we do is aperture extraction, but in a way that will enable future **dohydra** runs on actual data.

    Parameters
    ----------

    scalednames : list of str
        The names of flat images that have had shutter correction applied
    
    traceflat : str
        Name of the flat to be used for tracing apertures
    
    throughput : str, optional
        Name of the image to use for throughput correction

    Returns
    -------
    
    outputnames : list of str
        The names of the extracted flat field multispec images
    
    outputscales : list of float 
        List of scaling applied to each flat by **dohydra**. This is returned so that it can be undone later.

    Notes
    -----
    **dohydra** automatically scales each extracted flat so that the mean is 1. This messes up our careful shutter correction/scaling so we need to keep track of these scales for future use.
    
    """
    
    oldlog = iraf.hydra.logfile
    print 'Doing initial flat aperture extraction...'
    idtable = os.path.basename(iraf.dohydra.apidtable)
    print '\tusing idtable {}'.format(idtable)
    outputnames = []
    outputscales = []
    for flat in scalednames:
        print '\trunning {}'.format(flat)
        iraf.hydra.logfile = '{}.log'.format(flat)
        try:
            iraf.dohydra('ffTmp.fits',
                         apref=traceflat,
                         flat=flat,
                         through=throughput,
                         readnoise=3.9,
                         gain=0.438,
                         fibers=109,
                         width=5,
                         minsep=1,
                         maxsep=10,
                         apidtable=APIDTABLE,
                         scatter=False,
                         fitflat=False,
                         clean=False,
                         dispcor=False,
                         savearc=False,
                         skyalign=False,
                         skysubt=False,
                         skyedit=False,
                         savesky=False,
                         splot=False,
                         redo=False,
                         update=False,
                         batch=False,
                         listonl=False,
                         Stdout=1)
        except iraf.IrafError:
            '''For some reason, if you run this from the shell
            (non-interactively) dohydra will always fail the first
            time when trying to copy the database file. If/when this
            happens we'll just return false so main() knows just to
            try it again.
            '''
            print 'Fucked up, trying again'
            return False, False
        f = open('{}.log'.format(flat),'r')
        o = f.readlines()
        for dump in o:
            if 'bscale' in dump:
                scale = dump.split()[-1]
                outputscales.append(float(scale))
            if 'Create the normalized response' in dump:
                outname = '{}.fits'.format(dump.split()[-1])
                if outname not in outputnames:
                    outputnames.append('{}.fits'.format(dump.split()[-1]))
        os.remove('ffTmp.ms.fits')
        f.close()

    print "Extracted flats:"
    for out,scale in zip(outputnames,outputscales):
        print '\t{}  {}'.format(out,scale)
    iraf.hydra.logfile = oldlog
    return outputnames, outputscales

def stitch_flats(outputnames,pivots,outstring):
    """Take a list of multispec files and a list of pivots and stitch together a master flat.

    The pivot values are inclusive, so if pivots = [72] then the master flat will contain fibers 1 - 72 from flat #1 and 73 - 109 from flat #2.

    Parameters
    ----------
    
    outputnames : list of str
        Names of the multispec files to stitch together
    
    pivots : list of int
        The fibers that form the borders of the stitch. If outputnames is length N, pivots must be length N - 1
    
    outstring : str
        The special, IRAF scrunch string used to identify intermediate files associated with a **dohydra** run

    Returns
    -------
    
    mastername : str
        The name of the stitched master multispec flat

    Notes
    -----
    The specifics of outstring depend on system and IRAF distribution. See :meth:`GradPak_flatfu.get_scrunch`

    """
    pivots = [0] + pivots + [109]

    tmpfiles = []
    print 'Extracting flat apertures...'
    for i, flat in enumerate(outputnames):
        print '\ttaking {} from {} to {}'.format(flat,pivots[i]+1,pivots[i+1])
        name = 'tmp{}'.format(flat)
        iraf.scopy(flat,name,
                   apertur='{}-{}'.format(pivots[i]+1,pivots[i+1]),
                   w1='INDEF',
                   w2='INDEF',
                   format='multispec',
                   verbose=False)
        tmpfiles.append(name)


    mastername = 'dFlat_master{}.ms.fits'.format(outstring)
    print 'Stitching together master flat {}'.format(mastername)    
        
    iraf.scombine(','.join(tmpfiles),mastername,
                  apertur='',
                  group='apertures',
                  first=True,
                  w1='INDEF',
                  w2='INDEF',
                  dw='INDEF',
                  nw='INDEF',
                  log=False,
                  scale='none',
                  zero='none',
                  weight='none',
                  logfile='flatfu.log')

    for tmp in tmpfiles:
        os.remove(tmp)
    return mastername

def get_scrunch(flatname, msname):
    """Take in two strings and return the difference between them, modulo the difference between .ms.fits and .fits. 

    This is used to figure out what scrunching scheme IRAF used to name the aperture-extracted flats. This function exists because IRAF has different behavior when it comes to scrunching on different systems.

    Parameters
    ----------
    
    flatname : str
        Name of the raw, un-extracted file

    msname : str
        Name of the same file, after extraction by **dohydra**

    Returns
    -------

    scrunchstring : str
        The scrunch string used by this systems's IRAF

    """
    basefname = flatname.split('.fits')[0]
    basemsname = msname.split('.ms.fits')[0]
    
    return basemsname.replace(basefname,'')

def mean_scale(mslist,scalelist):
    """Take in a list of fits file names and scale each file by the corresponding value in the scalelist.

    When constructing the aperture-extracted flat that will be applied to all data apertures IRAF's last step is to normalize the entire multispec file to a unity mean. Because we will be stitching together a few different flats we need to undo this scaling so the relative strengths of the individual flats is maintained. This assumes that these flats were already properly scaled (see :meth:GradPak_flatfu.shutter_correction).

    Parameters
    ----------
    
    mslist : list of str
        Names of files to have **dohydra** scaling removed
    
    scalelist : list of float
        Scales applied to flats by **dohydra**

    Returns
    -------
    None : 
        The input files are overwritted with **dohydra**'s scaling removed

    """
    print 'Undoing IRAF flat scaling...'
    mshdulist = [pyfits.open(i)[0] for i in mslist]
    for h, name, scale in zip(mshdulist, mslist, scalelist):
        print '\t{:} scale: {:5.4f}'.format(name,scale)
        h.data *= scale
        h.writeto(name,clobber=True)

    return

def fit_flat(mastername):
    """Remove average spectra signature form multispec flat field
    
    This mimics the *fitflat* parameter in **dohydra**. The average spectrum is computed across all fibers and then divided into each fiber.

    Parameters
    ----------

    mastername : str
        Name of the flat to fit

    Returns
    -------
    None : 
        The input image is overwritted with the spectrally normalized version of itself.

    Notes
    -----
    Because of the different fiber sizes it is very likely that there will be some systematic residual spectral signature in the master flat. This is OK because you'll get rid of that during flux calibration.
    
    """
    print 'Dividing out average spectrum'
    h = pyfits.open(mastername)[0]
    h.data /= np.mean(h.data,axis=0)
    h.writeto(mastername,clobber=True)

    return

def normalize(mastername):
    """Normalize a fits file so the mean over all pixels is 1.

    This function exists to keep reduction as close as possible to IRAF's own routines. In msresp1d the last step is to normalize the flat so the mean is one. In :meth:GradPak_flatfu.mean_scale we undid this scaling so that the flats could be stitched together. This function recomputes this normalization to the stitched together flat.

    Parameters
    ----------

    mastername : str
        Name of the flat to be normalized

    Returns
    -------

    None : None
        The input file is overwritten with a version of itself with total flux = 1

    """
    print 'Renormalizing final flat'
    h = pyfits.open(mastername)[0]
    h.data /= np.mean(h.data)
    h.writeto(mastername,clobber=True)

    return

def parse_input(inputlist):
    """Parse arguments given to this module by the shell into useful function arguments and options

    Parameters
    ----------
    
    inputlist : list of str
        Probably sys.argv

    Returns
    -------
    
    flat_list : list of str
        Names of flats to be stitched together

    pivot_list : list of int
        Pivot points about which to stitch the flats

    traceflat : str
        Flat used for aperture extraction

    throughput : st
        Flat used for throughput correction. Can be an empty string (no special throughput file is used).

    fitflat : bool
        If True, fit out the average spectrum from the final flat

    """
    flat_list = []
    pivot_list = []
    fitflat = True
    traceflat = False
    throughput = ''

    for i, token in enumerate(inputlist):
        if token == '-nf':
            fitflat = False
        elif token == '-t':
            traceflat = inputlist[i+1]
            flat_list.append(inputlist[i+1])
        elif token == '-r':
            throughput = inputlist[i+1]
            del inputlist[i+1]
        elif '.fits' in token:
            if token not in flat_list:
                flat_list.append(token)
        else:
            try:
                pivot_list.append(int(token))
            except ValueError:
                print 'Could not parse option {}'.format(token)
    if not traceflat:
        traceflat = flat_list[0]
        
    return flat_list, pivot_list, traceflat, throughput, fitflat

def print_help():
    """Print usage help.
    """
    print """Calling syntax:
    
> GradPak_flatfu.py flat1.fits flat2.fits ... flatn.fits pivot

    Options:
            -t   Use the next flat as the trace flat. Default is 
                 the first flat given.

            -r   Use the next flat as the throughput (sky) image.
                 This flat is not used to make the master flat.

            -nf  Do not fit a average spectral function to the flats.
                 Equivalent to fitflat- in IRAF.

    Common Example:
             > GradPak_flatfu.py dFlat_4s.fits -t dFlat_1s.fits -r sFlat.fits -nf 43

             This will create a master flat consisting of apertures 1
             - 43 from dFlat_4s.fits and apertures 44 - 109 from
             dFlat_1s.fits. No functions will be fit to the flats but
             a fiber-by-fiber scaling will be applied based on the
             total counts in each aperture extracted from sFlat.fits
    """
    return

def main():
    """Parse the user inputs, check that there are the right number of pivot points, and run through the script.

    The steps are:
    
    #. Make dummy ffTmp.fits file
    #. Extract input flats using **dohydra**
    #. Undo the **dohydra** normalization
    #. Perform shutter correction
    #. Stitch flats together
    #. Optionally remove average spectral signature
    #. Normalize master flat to mean flux of 1
    #. Set external **dohydra** parameters so the user can call it from IRAF directly

    """
    flat_list, pivot_list, traceflat, throughput, fitflat = parse_input(sys.argv[1:])

    print 'Flat list is {}'.format(flat_list)
    print 'Pivots are {}'.format(pivot_list)
    print 'Tracing using {}'.format(traceflat)
    print 'Throughput file is {}'.format(throughput)
    
    if len(pivot_list) != len(flat_list) - 1:
        print "There must be one less pivot than flats"
        return 1

    '''Run the script'''
    # sl = setup_files(flat_list)
    make_tmp(pyfits.open(flat_list[0])[0])
    msl, scales = initial_run(flat_list, traceflat, throughput)
    if not msl:
        '''Here is where we catch IRAF being bad'''
        msl, scales = initial_run(msl, traceflat, throughput)
    outstring = get_scrunch(flat_list[0],msl[0])
    mean_scale(msl,scales)
    #msl = scale_spectra(msl)
    msl = shutter_correction(msl)
    master = stitch_flats(msl,pivot_list,outstring)
    if fitflat:
        fit_flat(master)
    normalize(master)

    '''Finally, set some IRAF parameters to make it easier to run dohydra with
    the master flat'''
    pd = iraf.dohydra.getParDict()
    pd['apref'].set(traceflat)
    pd['flat'].set('dFlat_master.fits')
    pd['through'].set(throughput)
    iraf.dohydra.saveParList()

    return 0
            
if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "The request was made but it was not good"
        sys.exit(1)
    elif sys.argv[1] == '-h':
        sys.exit(print_help())
    # try:
    #     sys.exit(main())
    # except Exception as e:
    #     print "The request was made but it was not good"
    #     sys.exit(1)
    sys.exit(main())
