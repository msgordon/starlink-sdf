starlink.sdf
============

starlink.sdf is a Python extension to Starlink that allows
for conversion of Starlink HDS (.sdf, .ndf) format files.
It does not require the Starlink software suite to be installed.
This project is based on `starlink-pyhds <https://github.com/sfgraves/starlink-pyhds>`_, a Python C-extension that liberates NDF file access from Starlink.

Note that this package is not well-tested or mature.  It only provides access
to JCMT/SCUBA .sdf files, and the conversion is not complete.  The NDF HISTORY
records are not preserved, though the main data, header, and WCS information
are extracted correctly.  Currently, only the DATA_ARRAY, VARIANCE, and QUALITY maps are pulled from the HDS file.

I am not associated with Starlink or JCMT.  I just wanted a quick, light-weight, pythonic interface to the NDF/SDF file format that was free of the entire Starlink architecture.


Installation
************

The package can be cloned from git and installed via::
  
  python setup.py install

or::
  
  python setup.py install --prefix<installation dir>

It can also be installed directly from this repository (recommended) as::

  pip install git+https://github.com/msgordon/starlink-sdf.git

The necessary dependencies (astropy, `starlink-pyast <https://github.com/timj/starlink-pyast>`_, `starlink-pyhds <https://github.com/sfgraves/starlink-pyhds>`_) will be downloaded and installed from pypi.  Note that the Starlink AST and HDS dependencies may take some time to compile; however, a full Starlink software installation is not required.
  
Usage
*****
The interface is designed to match closely with astropy.io.fits::

  from starlink import sdf
  hdul = sdf.open(sdf_image_filename)

In this case, a fits.HDUList object is returned, in the same manner as
fits.open(). Similar convenience functions exist as well::

  hdr = sdf.getheader(sdf_image_filename)  # get default HDU, i.e. primary HDU's header
  hdr = sdf.getheader(sdf_image_filename, 0)  # get primary HDU's header
  hdr = sdf.getheader(sdf_image_filename, 1)  # get first HDU's header
  hdr = sdf.getheader(sdf_image_filename, 'VARIANCE')  # get HDU header with name 'VARIANCE'

Note that for the last example, only PRIMARY (DATA_ARRAY), VARIANCE, and QUALITY
from JCMT files are implemented.

The function getdata() gets the data of an HDU. Similar to getheader(), it only requires the input .sdf file name while the extension is specified through the optional arguments. It does have one extra optional argument, header. If header is set to True, this function will return both data and header, otherwise only data is returned::

  data = getdata(sdf_image_filename, 'VARIANCE')
  # get 1st extension's data AND header:
  data, hdr = getdata(sdf_image_filename, 1, header=True)

To get all the supported extensions in one shot, use sdf.open() as shown above.

Installation of this module also adds a command-line program to your Python bin path, sdf2fits.  This program takes an input the .sdf filename, and optionally an output fits file.  For more usage information, invoke it with the `-h` switch


   sdf2fits -h
   
   usage: sdf2fits [-h] [-c] filename [outfile]

   Convert sdf/ndf to fits
