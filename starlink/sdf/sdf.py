from starlink import Atl,Ast,hds
from astropy.io import fits
from astropy.wcs import WCS,FITSFixedWarning
import numpy as np
from pathlib import Path
from tempfile import NamedTemporaryFile
import argparse
#import warnings
#warnings.filterwarnings('ignore',category=FITSFixedWarning)

### These are default options for JCMT Starlink hds files
DEF_PATH_MAP = {0:
                {'data_path':('DATA_ARRAY','DATA'),
                 'header_path':('MORE','FITS'),
                 'wcs_path':('WCS','DATA'),
                 'origin_path':('DATA_ARRAY','ORIGIN'),
                 'label_path':('LABEL',),
                 'unit_path':('UNITS',)},
                1:
                {'data_path':('VARIANCE','DATA'),
                 'header_path':None,
                 'wcs_path':('WCS','DATA'),
                 'origin_path':('VARIANCE','ORIGIN'),
                 'label_path':('LABEL',),
                 'unit_path':('UNITS',)},
                2:
                {'data_path':('QUALITY','QUALITY','DATA'),
                 'header_path':None,
                 'wcs_path':('WCS','DATA'),
                 'origin_path':('QUALITY','QUALITY','ORIGIN'),
                 'label_path':('LABEL',),
                 'unit_path':None},
                3:
                {'data_path':('MORE','SMURF','EXP_TIME','DATA_ARRAY','DATA'),
                 'header_path':None,
                 'wcs_path':('WCS','DATA'),
                 'origin_path':('MORE','SMURF','EXP_TIME','DATA_ARRAY','ORIGIN'),
                 'label_path':('MORE','SMURF','EXP_TIME','LABEL'),
                 'unit_path':('MORE','SMURF','EXP_TIME','UNITS')},
                4:
                {'data_path':('MORE','SMURF','WEIGHTS','DATA_ARRAY','DATA'),
                 'header_path':None,
                 'wcs_path':('WCS','DATA'),
                 'origin_path':('MORE','SMURF','WEIGHTS','DATA_ARRAY','ORIGIN'),
                 'label_path':('MORE','SMURF','WEIGHTS','LABEL'),
                 'unit_path':('MORE','SMURF','WEIGHTS','UNITS')},
                }
NAME_PATH_MAP = {'DATA':DEF_PATH_MAP[0],
                 'DATA_ARRAY':DEF_PATH_MAP[0],
                 'PRIMARY':DEF_PATH_MAP[0],
                 '0':DEF_PATH_MAP[0],
                 'VARIANCE':DEF_PATH_MAP[1],
                 '1':DEF_PATH_MAP[1],
                 'QUALITY':DEF_PATH_MAP[2],
                 '2':DEF_PATH_MAP[2],
                 'EXP_TIME':DEF_PATH_MAP[3],
                 '3':DEF_PATH_MAP[3],
                 'WEIGHTS':DEF_PATH_MAP[4],
                 '4':DEF_PATH_MAP[4]
                 }

def hds_path(root, path):
    """Recursively find on root"""
    if isinstance(path,str):
        # this is the end
        if path == '.':
            return root
        return root.find(path)
    elif len(path) == 1:
        return root.find(path[0])
    else:
        if path[0] == '.':
            node = root
        else:
            node = root.find(path[0])
        return hds_path(node,path[1:])


def _getext(*args, ext=None, extname=None, extver=None):
    # This is copied from astropy.io.fits
    err_msg = ('Redundant/conflicting extension arguments(s): {}'.format(
        {'args': args, 'ext': ext, 'extname': extname,
         'extver': extver}))

    # This code would be much simpler if just one way of specifying an
    # extension were picked.  But now we need to support all possible ways for
    # the time being.
    if len(args) == 1:
        # Must be either an extension number, an extension name, or an
        # (extname, extver) tuple
        if fits.util._is_int(args[0]) or (isinstance(ext, tuple) and len(ext) == 2):
            if ext is not None or extname is not None or extver is not None:
                raise TypeError(err_msg)
            ext = args[0]
        elif isinstance(args[0], str):
            # The first arg is an extension name; it could still be valid
            # to provide an extver kwarg
            if ext is not None or extname is not None:
                raise TypeError(err_msg)
            extname = args[0]
        else:
            # Take whatever we have as the ext argument; we'll validate it
            # below
            ext = args[0]
    elif len(args) == 2:
        # Must be an extname and extver
        if ext is not None or extname is not None or extver is not None:
            raise TypeError(err_msg)
        extname = args[0]
        extver = args[1]
    elif len(args) > 2:
        raise TypeError('Too many positional arguments.')

    if (ext is not None and
            not (fits.util._is_int(ext) or
                 (isinstance(ext, tuple) and len(ext) == 2 and
                  isinstance(ext[0], str) and fits.util._is_int(ext[1])))):
        raise ValueError(
            'The ext keyword must be either an extension number '
            '(zero-indexed) or a (extname, extver) tuple.')
    if extname is not None and not isinstance(extname, str):
        raise ValueError('The extname argument must be a string.')
    if extver is not None and not fits.util._is_int(extver):
        raise ValueError('The extver argument must be an integer.')

    if ext is None and extname is None and extver is None:
        ext = 0
    elif ext is not None and (extname is not None or extver is not None):
        raise TypeError(err_msg)
    elif extname:
        if extver:
            ext = (extname, extver)
        else:
            ext = (extname, 1)
    elif extver and extname is None:
        raise TypeError('extver alone cannot specify an extension.')

    if isinstance(ext,(int,str)):
        return ext
    if ext in ('',None):
        return 0
    if len(ext) == 2:
        try:
            ext = int(ext[0])
        except ValueError:
            # ignore version and just return string
            return ext[0]

    return ext

def extract_wcs(data, wcs_path=('WCS','DATA'),**kwargs):
    """Extract WCS information from data.
    data can be a filename, Path, or starlink hds object.
    path is an iterable of hds paths to search.
       default wcs_path=('WCS','DATA') yields:
          data.find('WCS').find('DATA')
    """

    if isinstance(data,(str,Path)):
        # open data as hds
        data = hds.open(str(data), 'r')

    wcsnode = hds_path(data, wcs_path)
    wcsdat = wcsnode.get().astype(str).tolist()

    # prune linebreak '+'
    plus_idx = [i for i,line in enumerate(wcsdat) if line[0] == '+']
    for idx in plus_idx:
        wcsdat[idx-1] += wcsdat[idx][1:]
    wcsdat = np.delete(np.array(wcsdat),plus_idx)
    wcsstring = '\n'.join(wcsdat)

    # open tempfile for Ast.Channel
    with NamedTemporaryFile('w') as f:
        f.write(wcsstring)
        c = Ast.Channel()
        c.set('SourceFile=%s' % f.name)
        fs = c.read()

    # create hdu for FitsChan
    hdu = fits.PrimaryHDU()
    Atl.writefitswcs(fs, hdu)
    wcs = WCS(hdu.header)
    return wcs

def extract_header(data, header_path=('MORE','FITS'),
                   wcs_path=('WCS','DATA'),
                   origin_path=('DATA_ARRAY','ORIGIN'),
                   label_path=('LABEL',),
                   unit_path=('UNITS',),
                   with_wcs=True,**kwargs):
    """Extract HEADER information from data.
    data can be a filename, Path, or starlink hds object.
    path is an iterable of hds paths to search.
       default header_path=('MORE','FITS') yields:
          data.find('MORE').find('FITS')
    if with_wcs == True, return header with WCS prepended
    """
    if isinstance(data,(str,Path)):
        # open data as hds
        data = hds.open(str(data), 'r')
        
    if header_path:      
        hnode = hds_path(data, header_path)
        hdat = hnode.get().astype(str)
        hstr = '\n'.join(hdat)
        h = fits.Header.fromstring(hstr, sep='\n')
    else:
        h = fits.Header()

    # add LABEL and BUNIT
    meta = {}
    if label_path:
        lnode = hds_path(data, label_path)
        label = lnode.get().decode()
        meta['LABEL'] = (label,'Label of the primary array')
    if unit_path:
        bnode = hds_path(data, unit_path)
        bunit = bnode.get().decode()
        meta['BUNIT'] = (bunit,'Units of the primary array')
    if kwargs.get('data_path'):
        extname = kwargs['data_path'][0]
        meta['EXTNAME'] = extname
        if extname == 'VARIANCE' and meta.get('BUNIT'):
            meta['BUNIT'] = ('(%s)**2'%bunit,'Units of the Variance array')
        if extname == 'DATA_ARRAY':
            meta['EXTNAME'] = 'PRIMARY'
        if extname == 'MORE':
            meta['EXTNAME'] = kwargs['data_path'][2]
    h.update(meta)

    if with_wcs:
        wcs = extract_wcs(data, wcs_path)
        hwcs = wcs.to_header()

        # add origin keywords
        onode = hds_path(data, origin_path)
        odata = onode.get().tolist()
        ohdr = {'LBOUND%i'%i:(x,'Pixel origin along axis %i'%i) \
                for i,x in enumerate(odata,start=1)}
        hwcs.update(ohdr)

        # finally, append original header
        hwcs.update(h)
        h = hwcs
        
    return h

def extract_data(data, data_path=('DATA_ARRAY','DATA'),
                 header_path=('MORE','FITS'),
                 wcs_path=('WCS','DATA'),
                 origin_path=('DATA_ARRAY','ORIGIN'),
                 label_path=('LABEL',),
                 unit_path=('UNITS',),
                 as_hdu=False,
                 with_header=False,
                 with_wcs=True):
    """Extract the data array.
    data can be a filename, Path, or starlink hds object.
    path is an iterable of hds paths to search.
       default data_path=('DATA_ARRAY','DATA') yields:
          data.find('DATA_ARRAY').find('DATA')
    if as_hdu == True, return ImageHDU, else return ndarray.
    """

    if isinstance(data,(str,Path)):
        # open data as hds
        data = hds.open(str(data), 'r')

    dnode = hds_path(data, data_path)
    data_array = dnode.get()

    # find bad values if any
    dtype = dnode.type
    badvalue = hds.getbadvalue(dtype)

    #mask = data_array == badvalue
    #masked_data = np.ma.MaskedArray(data_array, mask=mask)
    try:
        data_array[data_array == badvalue] = np.nan
    except ValueError:
        # do not assign badvalue
        pass
        
    if with_header:
        header = extract_header(data, header_path=header_path,
                                wcs_path=wcs_path,
                                origin_path=origin_path,
                                label_path=label_path,
                                unit_path=unit_path,
                                with_wcs=with_wcs)
        if as_hdu:
            extname = header.get('EXTNAME') if header.get('EXTNAME') \
                      else data_path[0]
            if extname == 'DATA_ARRAY':
                extname = 'PRIMARY'
            if extname == 'MORE':
                extname = data_path[2]
            hdu = fits.ImageHDU(data=data_array,header=header,name=extname)
            return hdu
        else:
            return data_array, header
    else:
        if as_hdu:
            header = extract_header(data, header_path=header_path,
                                    wcs_path=wcs_path,
                                    origin_path=origin_path,
                                    label_path=label_path,
                                    unit_path=unit_path,
                                    with_wcs=with_wcs)
            extname = header.get('EXTNAME') if header.get('EXTNAME') \
                      else data_path[0]
            if extname == 'DATA_ARRAY':
                extname = 'PRIMARY'
            if extname == 'MORE':
                extname = data_path[2]
            hdu = fits.ImageHDU(data=data_array,header=header,name=extname)
            return hdu
        else:
            return data_array
    

def _getopts(ext=None):
    ext = _getext(ext)

    if isinstance(ext,int):
        opts = DEF_PATH_MAP.get(ext)
    else:
        opts = NAME_PATH_MAP.get(ext)

    if opts is None:
        raise InputError('Extension %s not understood. Only DATA, VARIANCE, QUALITY, EXP_TIME, and WEIGHTS are implemented' % str(ext))

    return opts
        
def getdata(filename, *args, header=None, as_hdu=False):
    """High-level function to get data from sdf (and optionally the header)"""
    opts = _getopts(*args)

    if header:
        opts['with_header'] = True
    if as_hdu:
        opts['as_hdu'] = True

    return extract_data(filename, **opts)

def getheader(filename, *args):
    """High-level function to get header from sdf"""
    opts = _getopts(*args)
    return extract_header(filename, **opts)

def sdfopen(filename):
    """High-level function open sdf as hdulist"""
    exts = ('DATA','VARIANCE','QUALITY','EXP_TIME','WEIGHTS')

    hd = hds.open(filename,'r')

    hdlist = [getdata(hd,ext,as_hdu=True) for ext in exts]
    hdlist[0] = fits.PrimaryHDU(data=hdlist[0].data,
                                header=hdlist[0].header)
    hdulist = fits.HDUList(hdus=hdlist)
    return hdulist


def main():
    parser = argparse.ArgumentParser(description='Convert sdf/ndf to fits')
    parser.add_argument('filename',type=str,help='Input sdf file')
    parser.add_argument('outfile',type=str,nargs='?',default=None,
                        help='Optional output file. If not specified, file will have input name with .fits extension')
    parser.add_argument('-c','--clobber',dest='clobber',action='store_true',
                        help='Clobber/overwrite output file.')
    args = parser.parse_args()

    if args.outfile is None:
        args.outfile = str(Path(args.filename).with_suffix('.fits'))

    hdu = sdfopen(args.filename)
    hdu.info()
    try:
        hdu.writeto(args.outfile,overwrite=args.clobber)
        print('-> %s' % args.outfile)
    except OSError:
        print("ERROR:  File %s already exists. Use --clobber to overwrite." % args.filename)

if __name__ == '__main__':
    main()
    
