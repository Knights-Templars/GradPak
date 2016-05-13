The files in this directory are very useful when planning and executing a
GradPak program. They are:

 GradPak_fibermap.pdf - A figure showing the loction of all the active GradPak
 fibers. Fibers are labeled with the canonical numbering scheme.

 gradpak_w_sky_comp.tpl - A ds9 template file of the GradPak IFU. Very useful
 for planning your pointings.

 gradpak_w_sky_astrometry_table.txt - Very similar to gradpak_w_sky_comp.tpl,
 but in a slightly different format. Eric H. needed this at some point so
 maybe you will too.

 gradpak.iraf - Aperture ID table for use with IRAF reduction
 tasks. Specifically aperture identification and extraction. This version of
 the ApID table is only for people too boring/afraid to unleash the power of
 gradpak_sizes.iraf.

 gradpak_sizes.iraf - Aperture ID table for use with IRAF reduction tasks. In
 this table each fiber size is assigned a "beam" with number equal to the
 fiber size (in microns) divided by 100. The sky fibers are identified as the
 corrseponding beam number times 11. For example, a 3" fiber is assigned to
 beam 3 and a 4" sky fiber is assigned to beam 44. This makes sky subtraction
 a breeze. In fact, if you're going to use GradPak_skysub you need to use
 gradpak_sizes.iraf.
