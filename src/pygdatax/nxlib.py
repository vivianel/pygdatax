#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:37:16 2019

@author: achennev
"""

import nexusformat.nexus as nx
from nexusformat.nexus.tree import NeXusError
import datetime
import os
import numpy as np
import fabio
from scipy.io import loadmat, matlab
NXREAD_VERSION = '0.0'


def build_nexus_from_mat(filename):
    # TO FINISH
    def _check_keys(di):
        """
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        """
        for key in di:
            if isinstance(di[key], matlab.mio5_params.mat_struct):
                di[key] = _todict(di[key])
            elif isinstance(di[key], np.ndarray):
                subDict = []
                for element in di[key]:
                    subDict.append(_todict(element))
                di[key] = subDict
        return di

    def _todict(matobj):
        """
        A recursive function which constructs from matobjects nested dictionaries
        """
        di = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, matlab.mio5_params.mat_struct):
                di[strg] = _todict(elem)
            else:
                di[strg] = elem
        return di

    data = loadmat(filename, struct_as_record=False, squeeze_me=True)
    d = _check_keys(data)
    d.pop('__version__')
    d.pop('__globals__')
    d.pop('__header__')
    c = d.pop('c')
    for i, data in enumerate(c):
        key = 'c{%i}' % (i)
        d[key] = data

    root = nx.NXroot()
    for key, entry_n in zip(d, range(len(d))):
        entry_key = 'entry' + str(entry_n)
        entry = nx.NXentry(name=entry_key, attrs={'default': 'data'})
        m = d[key]['m']
        data = nx.NXdata(name='data')
        leg = d[key]['leg']
        if d[key]['et']['type'] == 'XY':
            data.nxsignal = nx.NXfield(m, name='m')
        elif d[key]['et']['type'] == 'Y':
            data.nxsignal = nx.NXfield(m[:, 1], name='y', units=leg['y'])
            data.nxerrors = m[:, 2]
            data.nxaxes = nx.NXfield(m[:, 0], name='x', units=leg['x'])
            data.insert(nx.NXfield(m[:, 3], name='x_errors', units=leg['x']))
        else:
            pass
        entry.insert(data)
        root.insert(entry)

    root.attrs['default'] = get_last_entry_key(root)
    new_name = filename.split('.')[0]
    new_name += '.nxs'
    root.save(new_name, mode='w')
    return root


def build_nexus_from_txt(filename, fields_name=None, comments='#', skiprows=0, delimiter='\t'):
    mat = np.loadtxt(filename)  # , comments=comments, skiprows=skiprows, delimiter=delimiter)
    root = nx.NXroot()
    root.attrs['default'] = 'entry0'
    entry = nx.NXentry()
    entry.attrs['default'] = 'data'
    entry.attrs['version'] = '1.0'
    entry.data = nx.NXdata()
    entry.data.nxsignal = nx.NXfield(mat[:, 1], name='y')
    entry.data.nxaxes = nx.NXfield(mat[:, 0], name='x', attrs={'uncertainties': 'x_errors'})
    if mat.shape[1] > 2:
        entry.data.nxerrors = mat[:, 2]
    if mat.shape[1] > 3:
        entry.data.insert(mat[:, 3], name='x_errors')
    root.entry0 = entry
    new_name = filename.split('.')[0]
    new_name += '.nxs'
    root.save(new_name, mode='w')
    return root


def build_nexus_from_edf(filename):
    fileObj = fabio.open(filename)
    header = fileObj.header
    root = nx.NXroot()
    root.attrs['default'] = 'entry0'
    root.attrs['NX_class'] = b'NXroot'
    entry = nx.NXentry()
    entry.attrs['default'] = 'data'
    entry.attrs['version'] = '1.0'
    entry.title = nx.NXfield(value=header['title'])
    entry.run = nx.NXfield(value=header['run'])
    entry.definition = nx.NXfield(definition='NXsas')
    # entry.date=nx.NXfield(date=)
    # building instrument
    instrument = nx.NXinstrument(description='Xeuss')
    y_gap = float(header['s2bot']) + float(header['s2top'])
    x_gap = float(header['s2hr']) + float(header['s2hl'])
    # intrument/aperture
    aperture = nx.NXslit(x_gap=nx.NXfield(x_gap, attrs={'units': 'mm'}),
                         y_gap=nx.NXfield(y_gap, attrs={'units': 'mm'}))
    # TO DO : find good distance
    collimator = nx.NXcollimator(length=nx.NXfield(1200, units='mm'),
                                 distance=nx.NXfield(20, units='mm'),
                                 s1bot=nx.NXfield(float(header['s1bot']), attrs={'units': 'mm'}),
                                 s1top=nx.NXfield(float(header['s1top']), attrs={'units': 'mm'}),
                                 s1hl=nx.NXfield(float(header['s1hl']), attrs={'units': 'mm'}),
                                 s1hr=nx.NXfield(float(header['s1hr']), attrs={'units': 'mm'}),
                                 s2bot=nx.NXfield(float(header['s1bot']), attrs={'units': 'mm'}),
                                 s2top=nx.NXfield(float(header['s1top']), attrs={'units': 'mm'}),
                                 s2hl=nx.NXfield(float(header['s2hl']), attrs={'units': 'mm'}),
                                 s2hr=nx.NXfield(float(header['s2hr']), attrs={'units': 'mm'}))
    # instrument/detector
    dist = float(header['SampleDistance']) * 1000
    detx = float(header['detx'])
    detz = float(header['detz'])
    pixSize1 = float(header['PSize_1']) * 1000
    pixSize2 = float(header['PSize_2']) * 1000
    x0 = float(header['Center_1'])
    y0 = float(header['Center_2'])
    detector = nx.NXdetector(data=fileObj.data,
                             distance=nx.NXfield(dist, attrs={'units': 'mm'}),
                             x_position=nx.NXfield(detx, attrs={'units': 'mm'}),
                             y_position=nx.NXfield(detz, attrs={'units': 'mm'}),
                             beam_center_x=nx.NXfield(x0, attrs={'units': 'pixel'}),
                             beam_center_y=nx.NXfield(y0, attrs={'units': 'pixel'}),
                             x_pixel_size=nx.NXfield(pixSize1, attrs={'units': 'mm'}),
                             y_pixel_size=nx.NXfield(pixSize2, attrs={'units': 'mm'}),
                             description='Pilatus 1M',
                             pixel_mask_applied=False,
                             pixel_mask=np.zeros_like(fileObj.data))
    # instrument/source
    wvl = float(header['WaveLength']) * 1e10
    sizeX = float(header['s1hr']) + float(header['s1hl'])
    sizeY = float(header['s1bot']) + float(header['s1top'])
    source = nx.NXsource(description='genix3D', radiation='x-ray',
                         incident_wavelength=nx.NXfield(wvl, attrs={'units': 'angstrom'}),
                         incident_wavelength_spread=0,
                         beam_size_x=nx.NXfield(sizeX, attrs={'units': 'mm'}),
                         beam_size_y=nx.NXfield(sizeY, attrs={'units': 'mm'}),
                         flux=nx.NXfield(float(header['pilai1']), attrs={'units': '1/s'}))
    entry.instrument = instrument
    entry.instrument.insert(detector)
    entry.instrument.insert(aperture)
    entry.instrument.insert(collimator)
    entry.instrument.insert(source)

    sample = nx.NXsample(sample_name=header['Comment'],
                         thickness=nx.NXfield(1.0, attrs={'units': 'mm'}),
                         transmission=float(header['pilroi1']),
                         x_position=nx.NXfield(float(header['x']), attrs={'units': 'mm'}),
                         y_position=nx.NXfield(float(header['z']), attrs={'units': 'mm'}),
                         om=nx.NXfield(float(header['om']), attrs={'units': 'deg'}),
                         phi=nx.NXfield(float(header['phi']), attrs={'units': 'deg'}),
                         rx=nx.NXfield(float(header['rx']), attrs={'units': 'deg'}),
                         ry=nx.NXfield(float(header['ry']), attrs={'units': 'deg'}),
                         temperature=nx.NXfield(float(header['Temperature']),
                                                attrs={'units': '°C'}),
                         count_time=nx.NXfield(float(header['count_time']),
                                               attrs={'units': 's'}),
                         description=header['Comment']
                         )
    entry.insert(sample)
    entry.data = nx.NXdata(attrs={'interpretation': b"image",
                                  'signal': "data"})
    root.entry0 = entry
    root.entry0.data.makelink(root.entry0.instrument.detector.data)
    new_name = filename.split('.')[0]
    new_name += '.nxs'
    try:
        root.save(new_name, mode='w')
        root.close()
    except NeXusError:
        print('error')
        # if os.path.exists(new_name):
        #     print('already here')
        #     os.remove(new_name)
        #     root.save(new_name, mode='w')
        #     root.unlock()
        # else:
        #     print('something else')

    return root

def build_rxnexus_from_edf(fileList, directbeam, outputFile):
    fileObj = fabio.open(fileList[0])
    header = fileObj.header
    root = nx.NXroot()
    root.attrs['default'] = 'entry0'
    root.attrs['NX_class'] = b'NXroot'
    entry = nx.NXentry()
    entry.title = nx.NXfield(value=header['title'])
    entry.run = nx.NXfield(value=header['run'])
    entry.definition = nx.NXfield(definition='NXsas')
    entry.attrs['default'] = 'data'
    entry.attrs['version'] = '1.0'
    entry.title = nx.NXfield(value=header['title'])
    entry.run = nx.NXfield(value=header['run'])
    entry.definition = nx.NXfield(definition='RX_Xeuss')
    # building instrument
    instrument = nx.NXinstrument(description='Xeuss')
    y_gap = float(header['s2bot']) + float(header['s2top'])
    x_gap = float(header['s2hr']) + float(header['s2hl'])
    # intrument/aperture
    aperture = nx.NXslit(x_gap=nx.NXfield(x_gap, attrs={'units': 'mm'}),
                         y_gap=nx.NXfield(y_gap, attrs={'units': 'mm'}))
    # TO DO : find good distance
    collimator = nx.NXcollimator(length=nx.NXfield(1200, units='mm'),
                                 distance=nx.NXfield(20, units='mm'),
                                 s1bot=nx.NXfield(float(header['s1bot']), attrs={'units': 'mm'}),
                                 s1top=nx.NXfield(float(header['s1top']), attrs={'units': 'mm'}),
                                 s1hl=nx.NXfield(float(header['s1hl']), attrs={'units': 'mm'}),
                                 s1hr=nx.NXfield(float(header['s1hr']), attrs={'units': 'mm'}),
                                 s2bot=nx.NXfield(float(header['s1bot']), attrs={'units': 'mm'}),
                                 s2top=nx.NXfield(float(header['s1top']), attrs={'units': 'mm'}),
                                 s2hl=nx.NXfield(float(header['s2hl']), attrs={'units': 'mm'}),
                                 s2hr=nx.NXfield(float(header['s2hr']), attrs={'units': 'mm'}))
    # instrument/detector
    dist = float(header['SampleDistance']) * 1000
    detx = float(header['detx'])
    detz = float(header['detz'])
    pixSize1 = float(header['PSize_1']) * 1000
    pixSize2 = float(header['PSize_2']) * 1000
    x0 = float(header['Center_1'])
    y0 = float(header['Center_2'])
    detector = nx.NXdetector(data=fileObj.data,
                             distance=nx.NXfield(dist, attrs={'units': 'mm'}),
                             x_position=nx.NXfield(detx, attrs={'units': 'mm'}),
                             y_position=nx.NXfield(detz, attrs={'units': 'mm'}),
                             beam_center_x=nx.NXfield(x0, attrs={'units': 'pixel'}),
                             beam_center_y=nx.NXfield(y0, attrs={'units': 'pixel'}),
                             x_pixel_size=nx.NXfield(pixSize1, attrs={'units': 'mm'}),
                             y_pixel_size=nx.NXfield(pixSize2, attrs={'units': 'mm'}),
                             description='Pilatus 1M',
                             pixel_mask_applied=False,
                             pixel_mask=np.zeros_like(fileObj.data))
    # instrument/source
    wvl = float(header['WaveLength']) * 1e10
    sizeX = float(header['s1hr']) + float(header['s1hl'])
    sizeY = float(header['s1bot']) + float(header['s1top'])
    source = nx.NXsource(description='genix3D', radiation='x-ray',
                         incident_wavelength=nx.NXfield(wvl, attrs={'units': 'angstrom'}),
                         incident_wavelength_spread=0,
                         beam_size_x=nx.NXfield(sizeX, attrs={'units': 'mm'}),
                         beam_size_y=nx.NXfield(sizeY, attrs={'units': 'mm'}),
                         flux=nx.NXfield(float(header['pilai1']), attrs={'units': '1/s'}))
    entry.instrument = instrument
    entry.instrument.insert(detector)
    entry.instrument.insert(aperture)
    entry.instrument.insert(collimator)
    entry.instrument.insert(source)

    sample = nx.NXsample(sample_name=header['Comment'],
                         thickness=nx.NXfield(1.0, attrs={'units': 'mm'}),
                         transmission=float(header['pilroi1']),
                         x_position=nx.NXfield(float(header['x']), attrs={'units': 'mm'}),
                         y_position=nx.NXfield(float(header['z']), attrs={'units': 'mm'}),
                         # om=nx.NXfield(float(header['om']), attrs={'units': 'deg'}),
                         phi=nx.NXfield(float(header['phi']), attrs={'units': 'deg'}),
                         rx=nx.NXfield(float(header['rx']), attrs={'units': 'deg'}),
                         ry=nx.NXfield(float(header['ry']), attrs={'units': 'deg'}),
                         temperature=nx.NXfield(float(header['Temperature']),
                                                attrs={'units': '°C'}),
                         count_time=nx.NXfield(float(header['count_time']),
                                               attrs={'units': 's'}),
                         description=header['Comment']
                         )
    entry.insert(sample)
    entry.data = nx.NXdata(attrs={'interpretation': b"image",
                                  'signal': "data"})
    root.entry0 = entry
    root.entry0.data.makelink(root.entry0.instrument.detector.data)
    new_name = filename.split('.')[0]
    new_name += '.nxs'
    try:
        root.save(outputFile, mode='w')
        root.close()
    except NeXusError:
        print('error')
        # if os.path.exists(new_name):
        #     print('already here')
        #     os.remove(new_name)
        #     root.save(new_name, mode='w')
        #     root.unlock()
        # else:
        #     print('something else')

    return root


def get_last_entry_key(root):
    """
    Returns the last NXentry key assuming that the entries are formated using the format : "entry0", "entry1", ....

    Args:
        root (NXroot): NXroot object
    Returns:
        str : the last entry key

    """
    last_key = 'entry0'
    i = 0
    for key in root.keys():
        if int(key[-1]) > i:
            last_key = key
            i = int(key[-1])
    return last_key


def create_new_entry(root):
    """
    Create a new NXentry which key is formated using the format : "entry0", "entry1", ....
    The new NXentry is a hard copy of the last NXentry

    Args:
        root (NXroot): NXroot object
    Returns:
        str : The key of the newly created entry

    """
    last_key = get_last_entry_key(root)
    new_entry = root[last_key].copy()
    new_key = 'entry' + str(int(last_key[-1]) + 1)
    root[new_key] = new_entry
    return new_key


def delete_last_entry(root):
    """
     Delete the last NXentry which key assuming that the entries are formated using the format : "entry0", "entry1", ....

     Args:
         root (NXroot): NXroot object
     Returns:
         None : None

     """
    last_key = get_last_entry_key(root)
    if last_key != 'entry0':
        del (root[last_key])
    last_key = get_last_entry_key(root)
    root.attrs['default'] = last_key
    return


def delete_all_entry(root):
    """
     Delete the all NXentry  except "entry0"

     Args:
         root (NXroot): NXroot object
     Returns:
         None : None

     """
    # root = nx.nxload(file, mode='rw')
    last_key = get_last_entry_key(root)
    while last_key != 'entry0':
        del (root[last_key])
        delete_last_entry(root)
        last_key = get_last_entry_key(root)
    last_key = get_last_entry_key(root)
    root.attrs['default'] = last_key
    return


def get_last_signal_key(root):
    """
    Returns the path of the signal associated to the default NXdata of a NXroot object

    Args:
        root (NXroot): NXroot object
    Returns:
        str : default signal path

    """
    key = get_last_entry_key(root)
    key += '/' + root[key].default
    key += '/' + root[key].signal
    return key


def function_performed(root, function_name):
    """
    Check if a treatment function has already been performed

    Args:
        root (NXroot): NXroot object
        function_name (str): treatment function name

    Returns:
        bool: True if performed, False if not

    """
    performed = False
    for entry_key in root:
        if 'process' in root[entry_key]:
            if root[entry_key + '/process/program'] == function_name:
                performed = True
                return performed
    return performed


def treatment_function(func):
    def wrapper(*args, **kwargs):
        file = args[0]
        root = loadfile(file, mode='rw')
        with root.nxfile:
            # check if function is already performed
            # TODO : solve problem when azimutal over the three detector
            # if nr.function_performed(root, func.__name__):
            #     print('function ' + func.__name__ + ' already performed.')
            #     return
            if 'new_entry' in kwargs:
                if kwargs['new_entry']:
                    last_key = create_new_entry(root)
                else:
                    last_key = get_last_entry_key(root)
            else:
                last_key = create_new_entry(root)
            root.attrs['default'] = last_key

            # entry = root[last_key]
            if 'new_entry' in kwargs:
                kwargs.pop('new_entry')
            args = list(args)
            args.pop(0)
            func(root, *args, **kwargs)
            last_key = get_last_entry_key(root)
            entry = root[last_key]
            if 'process' not in entry:
                entry.process = nx.NXprocess(program=func.__name__,
                                             arguments=str(kwargs),
                                             date=str(datetime.datetime.today()),
                                             version='nxread-'+NXREAD_VERSION
                                             )
            else:
                del root[last_key + '/' + 'process']
                entry.process = nx.NXprocess(program=func.__name__,
                                             arguments=str(kwargs),
                                             date=str(datetime.datetime.today()),
                                             version='nxread-' + NXREAD_VERSION
                                             )
            return

    # wrapper._original = func
    # wrapper.__name__ = func.__name__
    return wrapper


def loadfile(file, mode='rw'):
    name, extension = os.path.splitext(file)
    if extension == '.edf':
        build_nexus_from_edf(file)
        extension = '.nxs'
    root = nx.nxload(name + extension, mode=mode, recursive=False)
    return root

def save_as_txt(filename):
    root = nx.nxload(filename, mode='r')
    last_key = get_last_entry_key(root)
    data = root[last_key + '/data']
    signal_key = data.signal
    if isinstance(data[signal_key], nx.tree.NXlinkfield):
        signal_key = data.nxsignal.nxtarget
    else:
        signal_key = last_key + '/data/' + data.signal
    signal = root[signal_key].nxdata
    if len(signal.shape) > 1:
        print("Cannot save 2D data as .txt file")
        return

    # if 'axes' in root.entry1.data._attrs:

    # if 'errors' in root[last_key+'/data']:
    #     root[last_key+'/data/errors'] /= flux
