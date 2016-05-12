************
Installation
************

Hopefully this is the easiest part of the whole process. All of the special
GradPak stuff lives in a GitHub repo at https://github.com/eigenbrot/GradPak
so you can just clone it to wherever you want::

 > git clone https://github.com/eigenbrot/GradPak.git

This will result in a folder called GradPak in your current directory. The
structure is::

 GradPak
 |-- docs/
 |   |-- Source code for these document pages. Can be ignored.
 |
 |-- extras/
 |   |-- Things like fiber maps, ds9 templates, etc.
 |
 |-- GradPak code

All of the main python code you'll need for data reduction and analysis lives
in the root directory. The items in the extra directory will be very usefull
for planning and executing your observation run.

That's it! You are now all setup to enter the beautiful world of GradPak. If
you want, shoot an email over to eigenbrot at astro.wisc.edu and I'll keep you
in the loop on any updates. Please also send me any bug reports/feature
requests.


Set up PATH's
=============

It is perfectly reasonable to just call the GradPak code from the main repo
directory created above, but if you want to use the more advanced analysis
modules like :mod:`GradPak_plot` or :mod:`GradPak_bin` you'll need them in
your PYTHONPATH, so let's just set them up now.

C Shell
-------

Bust open your .cshrc file and add/modify the following lines::

 setenv PYTHONPATH EXISTING_PATH:PATH_TO_GRADPAK_REPO
 set path = ( PATH_TO_GRADPAK_REPO $path )

Of course, fill in values that make sense for your setup.

Bash
----

Add the following lines to your .bashrc::

 export PYTHONPATH="${PYTHONPATH}:PATH_TO_GRADPAK_REPO
 PATH=PATH_TO_GRADPAK_REPO:$PATH
