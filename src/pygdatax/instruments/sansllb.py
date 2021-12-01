#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:37:16 2019

@author: achennev
"""


import fabio
import nexusformat.nexus as nx
import numpy as np
import time
from pygdatax import flib
from pygdatax import nxlib
import os
NXREAD_VERSION = '0.0'


@nxlib.treatment_function(new_entry=False)
def set_beam_center(root, detector=0, x0=64, y0=64):
    # root = loadfile(file, mode='rw')
    # with root.nxfile:
    entry = root[nxlib.get_last_entry_key(root)]
    # if x0 is None or y0 is None:
    #     # read direct beam and find center around the center entered in entry
    #     x0 = entry.instrument.detector0.beam_center_x.nxdata
    #     y0 = entry.instrument.detector0.beam_center_y.nxdata
    # if direct_beam_file is not None:
    #     # if direct file is given, it finds the driect beam using barycenter in a
    #     # a 40x40 pixel box around (x0,y0)
    #     root_db = nxlib.loadfile(direct_beam_file, mode='r')
    #     m = root_db.entry0.data0.data.nxdata
    #     i1 = int(x0 - 40)
    #     i2 = int(x0 + 40)
    #     j1 = int(y0 - 40)
    #     j2 = int(y0 + 40)
    #     crop_m = m[j1:j2, i1:i2]
    #     x0, y0 = flib.find_direct_beam(m, corners=[i1, i2, j1, j2])
    entry['instrument/detector'+str(detector)+'/beam_center_x'] = x0
    entry['instrument/detector' + str(detector) + '/beam_center_y'] = y0
    return


# @nxlib.treatment_function(new_entry=False)
# TODO: handle inf and nan pixels
def azimutal_integration(root, detector=0, mask_file=None, x0=None, y0=None, bins=90,
                         x_pixel_size=None, y_pixel_size=None):
    last_key = nxlib.get_last_entry_key(root)
    entry = root[last_key]
    if x0 is None:
        x0 = entry['instrument/detector'+str(detector)+'/beam_center_x'].nxdata
    else:
        entry['instrument/detector'+str(detector)+'/beam_center_x'].nxdata = x0

    if y0 is None:
        y0 = entry['instrument/detector'+str(detector)+'/beam_center_y'].nxdata
    else:
        entry['instrument/detector'+str(detector)+'/beam_center_y'].nxdata = y0
    if x_pixel_size is None:
        x_pixel_size = entry['instrument/detector'+str(detector)+'/x_pixel_size']
    else:
        entry['instrument/detector'+str(detector)+'/x_pixel_size'] = x_pixel_size

    if y_pixel_size is None:
        y_pixel_size = entry['instrument/detector'+str(detector)+'/y_pixel_size']
    else:
        entry['instrument/detector'+str(detector)+'/y_pixel_size'] = y_pixel_size

    m = entry['data'+str(detector)].nxsignal.nxdata
    mask_data = np.zeros_like(m)

    if mask_file is not None:
        if mask_file.endswith('.nxs'):
            maskroot = nxlib.loadfile(mask_file, mode='r')
            mask_data = maskroot.data.nxdata
        elif mask_file.endswith('.edf'):
            mask_data = fabio.open(mask_file).data
    if 'data_errors' in entry['data'+str(detector)]:
        error = entry['data'+str(detector)].nxerrors.nxdata
    else:
        error = None
    del entry['data' + str(detector)]
    entry['instrument/detector' + str(detector) + '/pixel_mask_applied'] = True
    entry['instrument/detector' + str(detector) + '/pixel_mask'] = mask_data
    r, i, sigma, dr = flib.regiso(m, mask_data, x0, y0, x_pixel_size, y_pixel_size,
                                  bins, error=error)

    # new_entry.data = nx.NXdata(attrs={'interpretation': 'spectrum',
    #                                   'signal': 'I', 'axes': ['Q']})
    entry['data'+str(detector)] = nx.NXdata()
    # new_entry.data.I = nx.NXfield(i, units='counts')  # uncertainties='I_errors'
    entry['data'+str(detector)].nxsignal = nx.NXfield(i, name='counts', attrs={'interpretation': 'spectrum'})
    entry['data'+str(detector)].nxerrors = sigma
    entry['data'+str(detector)].nxaxes = nx.NXfield(r, name='r', attrs={'units': 'mm', 'uncertainties': 'r_errors'})
    entry['data'+str(detector)].r_errors = nx.NXfield(dr, attrs={'units': 'mm'})
    # new_entry.data.nxerrors = sigma
    # new_entry.process = nx.NXprocess(program='azimutal integration',
    #                                  sequence_index=1,
    #                                  date=str(datetime.datetime.today()))
    return


@nxlib.treatment_function
def azimutal_integration_multidetector(root: nx.NXroot, mask_file0: str = None, x0: float = None, y0: float = None, bins0: int = 90,
                                       x_pixel_size0: float = None, y_pixel_size0: float = None,
                                       mask_file1: str = None, x1: float = None, y1: float = None, bins1: int = 90,
                                       x_pixel_size1: float = None, y_pixel_size1: float = None,
                                       mask_file2: str = None, x2: float = None, y2: float = None, bins2: int = 90,
                                       x_pixel_size2: float = None, y_pixel_size2: float = None
                                       ) -> None:
    azimutal_integration(root, detector=0, mask_file=mask_file0, x0=x0, y0=y0, bins=bins0,
                         x_pixel_size=x_pixel_size0, y_pixel_size=y_pixel_size0)
    azimutal_integration(root, detector=1, mask_file=mask_file1, x0=x1, y0=y1, bins=bins1,
                         x_pixel_size=x_pixel_size1, y_pixel_size=y_pixel_size1)
    azimutal_integration(root, detector=2, mask_file=mask_file2, x0=x2, y0=y2, bins=bins2,
                         x_pixel_size=x_pixel_size2, y_pixel_size=y_pixel_size2)


@nxlib.treatment_function(new_entry=True)
def reduction2D(root: nx.NXroot, sub_file: str = None, norm_file: str = None,
                thickness: float = None, transmission: float = None, distance: list = [None, None, None]) -> None:
    """
    Subtract and normalize the sans 2D spectra according to Brûlet, A., Lairez, D., Lapp, A., & Cotton, J. P. (2007). Improvement of data treatment in small-angle neutron scattering. Journal of Applied Crystallography, 40(1), 165-177.
    The resulting spectra is a 2D spectra
    Args:
        root:
        sub_file (str): substraction package file
        norm_file (str): normalization package file
        thickness (float): sample thickness
        transmission (float): sample transmisssion
        distance (list): list of detector distances

    Returns:

    """
    entry = root[nxlib.get_last_entry_key(root)]

    def delta(u, a):
        if u == 1:
            v = 1
        else:
            v = (1 - x ** a) / (-a * np.log(x))
        return v

    def tr_theta(t, angle):
        if t == 1:
            t_th = 1
        else:
            t_th = t * (t ** (1 - 1 / np.cos(angle)) - 1) / (np.log(t) * (1 - 1 / np.cos(angle)))
        return t_th

    if transmission is None:
        transmission = entry.sample.transmission.nxdata
        # if transmission == 0:
        #     entry.sample.transmission = 1.0
    else:
        entry.sample.transmission = transmission
    if thickness is None:
        thickness = entry.sample.thickness.nxdata
    else:
        entry.sample.thickness = thickness

    for i in range(3):
        distance_key = 'instrument/detector' + str(i) + '/distance'
        if distance[i] is None:
            distance[i] = entry[distance_key].nxdata
        else:
            entry[distance_key].nxdata = distance[i]

        i_sample = root['entry0/data'+str(i)]
        i_sample.nxerrors = nx.NXfield(np.sqrt(np.abs(i_sample.data.nxdata)))
        monitor_sample = root['entry0/monitor3/integral'].nxdata
        time_sample = root['entry0/instrument/detector'+str(i)+'/count_time'].nxdata
        shape = i_sample.nxsignal.nxdata.shape

        # uncack the substraction package file
        default_params = get_default_reduction_parameters(root)
        if sub_file is not None:
            sub_root = nx.nxload(sub_file, mode='rw')
            if 'dark' in sub_root:
                i_dark = sub_root['dark/data'+str(i)]
                # i_dark.nxerrors = np.zeros(shape)
                i_dark.nxerrors = nx.NXfield(np.sqrt(np.abs(i_dark.data.nxdata)))
                # monitor_dark = sub_root['dark/monitor3/integral'].nxdata
                time_dark = sub_root['dark/instrument/detector'+str(i)+'/count_time'].nxdata
            else:
                i_dark = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_dark'))
                i_dark.nxerrors = np.zeros(shape)
                # monitor_dark = 1
                time_dark = 1
            if 'empty cell' in sub_root:
                i_ec = sub_root['empty cell/data'+str(i)]
                i_ec.nxerrors = nx.NXfield(np.sqrt(np.abs(i_ec.nxsignal.nxdata)))
                monitor_ec = sub_root['empty cell/monitor3/integral'].nxdata
                time_ec = sub_root['empty cell/instrument/detector'+str(i)+'/count_time'].nxdata
                trans_ec = sub_root['empty cell/sample/transmission'].nxdata
            else:
                i_ec = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_ec'))
                i_ec.nxerrors = np.zeros(shape)
                monitor_ec = 1
                time_ec = 1
                trans_ec = 1
            if 'empty beam' in sub_root:
                i_eb = sub_root['empty beam/data'+str(i)]
                i_eb.nxerrors = nx.NXfield(np.sqrt(np.abs(i_eb.data.nxdata)))
                monitor_eb = sub_root['empty beam/monitor3/integral'].nxdata
                time_eb = sub_root['empty beam/instrument/detector'+str(i)+'/count_time'].nxdata
            else:
                i_eb = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_eb'))
                i_eb.nxerrors = np.zeros(shape)
                monitor_eb = 1
                time_eb = 1

            # centers for regroupement
            if 'beam_center_x' in sub_root['parameters/detector'+str(i)]:
                x0 = sub_root['parameters/detector'+str(i)+'/beam_center_x'].nxdata
                entry['instrument/detector'+str(i)+'/beam_center_x'] = x0
            else:
                x0 = default_params['x0']
            if 'beam_center_y' in sub_root['parameters/detector'+str(i)]:
                y0 = sub_root['parameters/detector'+str(i)+'/beam_center_y'].nxdata
                entry['instrument/detector' + str(i) + '/beam_center_y'] = y0
            else:
                y0 = default_params['y0']
        else:
            i_dark = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_dark'))
            i_dark.nxerrors = np.zeros(shape)
            # monitor_dark = 1
            time_dark = 1
            i_ec = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_ec'))
            i_ec.nxerrors = np.zeros(shape)
            monitor_ec = 1
            time_ec = 1
            trans_ec = 1
            i_eb = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_eb'))
            i_eb.nxerrors = np.zeros(shape)
            monitor_eb = 1
            time_eb = 1
            x0 = default_params['x0']
            y0 = default_params['y0']

        x_pixel_size = entry['instrument/detector' + str(i) + '/x_pixel_size'].nxdata
        y_pixel_size = entry['instrument/detector' + str(i) + '/y_pixel_size'].nxdata
        solid_angle = x_pixel_size * y_pixel_size / (distance[i] ** 2)
        y, x = np.indices(shape, dtype='float')
        y = (y - y0)*y_pixel_size
        x = (x - x0)*x_pixel_size
        r = np.sqrt(x ** 2 + y ** 2)
        theta = np.arctan(r / distance[i])

        if monitor_eb == 1:  # no empty beam
            fb = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_eb'))
            fb.nxerrors = np.zeros(shape)
        else:
            fb = i_eb - i_dark*time_eb/time_dark
            fb /= monitor_eb
        trans_ec = trans_ec ** 0.5
        if monitor_ec == 1 and time_ec == 1 and trans_ec == 1:  # no empty cell given
            z_fec = nx.NXdata(nx.NXfield(np.zeros_like(i_sample.nxsignal), name='z_fec'))
            z_fec.nxerrors = np.zeros_like(i_sample.nxsignal)
        else:
            z_fec = (i_ec - i_dark*time_ec/time_dark) / monitor_ec - fb * np.power(trans_ec, 2 / np.cos(theta))
            z_fec /= trans_ec ** (1 / np.cos(theta)) * (trans_ec + tr_theta(trans_ec, theta))

        # substarct the contributions
        fs = (i_sample - i_dark*time_sample/time_dark) / monitor_sample
        fs -= fb * trans_ec ** (2 / np.cos(theta)) * transmission ** (1 / np.cos(theta))
        fs -= z_fec * tr_theta(trans_ec, theta) * ((transmission * trans_ec) ** (1 / np.cos(theta)) +
                                                   trans_ec * transmission)
        fs /= tr_theta(transmission, theta) * trans_ec * trans_ec ** (1 / np.cos(theta)) * 0.1 * thickness
        fs /= solid_angle
        fs /= np.cos(theta) ** 3
        # normalization by water
        if norm_file:
            if os.path.exists(norm_file):
                norm_root = nx.nxload(norm_file, mode='r')
                if 'water_substracted' in norm_root.keys():
                    i_water = norm_root['water_substracted/data'+str(i)]
                else:
                    norm_root.close()
                    treat_normalization_package(norm_file)
                    norm_root = nx.nxload(norm_file, mode='r')
                    i_water = norm_root['water_substracted/data' + str(i)]
                fs /= i_water
            else:
                print('normalization file not found')
        del entry['data'+str(i)]
        entry['data'+str(i)] = fs
        # q_scale(root.file_name, distance=distance, new_entry=False)
    return


def treat_normalization_package(norm_file):

    def delta(u, a):
        if u == 1:
            v = 1
        else:
            v = (1 - x ** a) / (-a * np.log(x))
        return v

    def tr_theta(t, angle):
        if t == 1:
            t_th = 1
        else:
            t_th = t * (t ** (1 - 1 / np.cos(angle)) - 1) / (np.log(t) * (1 - 1 / np.cos(angle)))
        return t_th
    if not os.path.exists(norm_file):
        print('normalization file not found')
        return
    norm_root = nx.nxload(norm_file, mode='rw')
    if 'water' not in norm_root.keys():
        print('no water file in the normalisation package')
        return
    if 'water_substracted' not in norm_root.keys():
        new_entry = norm_root['water'].copy()
        norm_root['water_substracted'] = new_entry

    water_entry = norm_root['water_substracted']
    transmission = water_entry.sample.transmission.nxdata

    for i in range(3):
        distance_key = 'instrument/detector' + str(i) + '/distance'
        distance = water_entry[distance_key].nxdata
        i_water = water_entry['data'+str(i)]
        i_water.nxerrors = nx.NXfield(np.sqrt(np.abs(i_water.data.nxdata)))
        monitor_sample = water_entry['monitor3/integral'].nxdata
        time_sample = water_entry['instrument/detector'+str(i)+'/count_time'].nxdata
        shape = i_water.nxsignal.nxdata.shape

        # uncack the substraction package file
        if 'dark' in norm_root:
            i_dark = norm_root['dark/data'+str(i)]
            # i_dark.nxerrors = np.zeros(shape)
            i_dark.nxerrors = nx.NXfield(np.sqrt(np.abs(i_dark.data.nxdata)))
            # monitor_dark = sub_root['dark/monitor3/integral'].nxdata
            time_dark = norm_root['dark/instrument/detector'+str(i)+'/count_time'].nxdata
        else:
            i_dark = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_dark'))
            i_dark.nxerrors = np.zeros(shape)
            # monitor_dark = 1
            time_dark = 1
        if 'empty cell' in norm_root:
            i_ec = norm_root['empty cell/data'+str(i)]
            i_ec.nxerrors = nx.NXfield(np.sqrt(np.abs(i_ec.nxsignal.nxdata)))
            monitor_ec = norm_root['empty cell/monitor3/integral'].nxdata
            time_ec = norm_root['empty cell/instrument/detector'+str(i)+'/count_time'].nxdata
            trans_ec = norm_root['empty cell/sample/transmission'].nxdata
        else:
            i_ec = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_ec'))
            i_ec.nxerrors = np.zeros(shape)
            monitor_ec = 1
            time_ec = 1
            trans_ec = 1
        if 'empty beam' in norm_root:
            i_eb = norm_root['empty beam/data'+str(i)]
            i_eb.nxerrors = nx.NXfield(np.sqrt(np.abs(i_eb.data.nxdata)))
            monitor_eb = norm_root['empty beam/monitor3/integral'].nxdata
            time_eb = norm_root['empty beam/instrument/detector'+str(i)+'/count_time'].nxdata
        else:
            i_eb = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_eb'))
            i_eb.nxerrors = np.zeros(shape)
            monitor_eb = 1
            time_eb = 1

        # centers for regroupement
        if 'beam_center_x' in norm_root['parameters/detector'+str(i)]:
            x0 = norm_root['parameters/detector'+str(i)+'/beam_center_x'].nxdata
            water_entry['instrument/detector'+str(i)+'/beam_center_x'] = x0
        else:
            x0 = water_entry['instrument/detector'+str(i)+'/beam_center_x'].nxdata
        if 'beam_center_y' in norm_root['parameters/detector'+str(i)]:
            y0 = norm_root['parameters/detector'+str(i)+'/beam_center_y'].nxdata
            water_entry['instrument/detector' + str(i) + '/beam_center_y'] = y0
        else:
            y0 = water_entry['instrument/detector' + str(i) + '/beam_center_y']
        x_pixel_size = water_entry['instrument/detector' + str(i) + '/x_pixel_size'].nxdata
        y_pixel_size = water_entry['instrument/detector' + str(i) + '/y_pixel_size'].nxdata
        y, x = np.indices(shape, dtype='float')
        y = (y - y0)*y_pixel_size
        x = (x - x0)*x_pixel_size
        r = np.sqrt(x ** 2 + y ** 2)
        theta = np.arctan(r / distance)
        solid_angle = x_pixel_size * y_pixel_size / (distance ** 2)
        if monitor_eb == 1:  # no empty beam
            fb = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_eb'))
            fb.nxerrors = np.zeros(shape)
        else:
            fb = i_eb - i_dark*time_eb/time_dark
            fb /= monitor_eb
        trans_ec = trans_ec ** 0.5
        if monitor_ec == 1 and time_ec == 1 and trans_ec == 1:  # no empty cell given
            z_fec = nx.NXdata(nx.NXfield(np.zeros_like(i_water.nxsignal), name='z_fec'))
            z_fec.nxerrors = np.zeros_like(i_water.nxsignal)
        else:
            z_fec = (i_ec - i_dark*time_ec/time_dark) / monitor_ec - fb * np.power(trans_ec, 2 / np.cos(theta))
            z_fec /= trans_ec ** (1 / np.cos(theta)) * (trans_ec + tr_theta(trans_ec, theta))

        # substarct the contributions
        fs = (i_water - i_dark*time_sample/time_dark) / monitor_sample
        fs -= fb * trans_ec ** (2 / np.cos(theta)) * transmission ** (1 / np.cos(theta))
        fs -= z_fec * tr_theta(trans_ec, theta) * ((transmission * trans_ec) ** (1 / np.cos(theta)) +
                                                   trans_ec * transmission)
        fs /= tr_theta(transmission, theta) * trans_ec * trans_ec ** (1 / np.cos(theta)) * 0.1
        fs /= solid_angle
        fs /= np.cos(theta) ** 3
        del water_entry['data'+str(i)]
        water_entry['data'+str(i)] = fs
        # q_scale(root.file_name, distance=distance, new_entry=False)
        norm_root.close()
    return


def get_default_reduction_parameters(root):
    """
    Return the default parameters for azimutal integration stored within the NXroot object for the three detectors.
    These parameters are the beam center position and the default number of bins for each detector.

    The default bin nubre if defined as the diagonal length for the central detector and the number of tubes for the
    wide angle detectors.
    Args:
        root: NXroot

    Returns:
        parameters dictionnary e.g. {'x0':..., 'y0':..., 'bins0':...,
                                    'x1':..., 'y1':..., 'bins1':...,
                                    'x2':..., 'y2':..., 'bins2':...,}
    """
    entry = root[nxlib.get_last_entry_key(root)]
    x0 = entry.instrument.detector0.beam_center_x.nxdata
    y0 = entry.instrument.detector0.beam_center_y.nxdata
    bins0 = 90
    x1 = entry.instrument.detector1.beam_center_x.nxdata
    y1 = entry.instrument.detector1.beam_center_y.nxdata
    bins1 = 16
    x2 = entry.instrument.detector2.beam_center_x.nxdata
    y2 = root.entry0.instrument.detector2.beam_center_y.nxdata
    bins2 = 16
    dic = {'x0': x0, 'y0': y0, 'bins0': bins0,
           'x1': x1, 'y1': y1, 'bins1': bins1,
           'x2': x2, 'y2': y2, 'bins2': bins2
           }
    return dic


@nxlib.treatment_function(output_file=True)
def make_reduction_package(output_file,
                           dark_file=None, empty_cell_file=None, direct_beam_file=None,
                           water_file=None,
                           mask_file0=None, x0=None, y0=None, bins0=90,
                           mask_file1=None, x1=None, y1=None, bins1=16,
                           mask_file2=None, x2=None, y2=None, bins2=16,
                           ):
    root = nx.NXroot()
    root.save(output_file, mode='w')
    # nxfile = nx.NXFile(output_file, 'w')
    # root = nxfile.readfile()
    if dark_file is not None:
        dark_root = nx.nxload(dark_file, mode='r')
        nxlib.copy_entry(root, dark_root, 'dark', 'entry0')
        dark_root.close()
        def_params = get_default_reduction_parameters(dark_root)

    if empty_cell_file is not None:
        ec_root = nx.nxload(empty_cell_file, mode='r')
        nxlib.copy_entry(root, ec_root, 'empty cell', 'entry0')
        ec_root.close()
        def_params = get_default_reduction_parameters(ec_root)

    if direct_beam_file is not None:
        db_root = nx.nxload(direct_beam_file, mode='r')
        nxlib.copy_entry(root, dark_root, 'empty beam', 'entry0')
        db_root.close()
        def_params = get_default_reduction_parameters(db_root)

    if water_file is not None:
        water_root = nx.nxload(water_file, mode='r')
        nxlib.copy_entry(root, water_root, 'water', 'entry0')
        water_root.close()
        def_params = get_default_reduction_parameters(water_root)

    params_entry = nx.NXgroup()
    root.insert(params_entry, name='parameters')
    for i, maskfile, x, y, bin in zip([0, 1, 2],
                                      [mask_file0, mask_file1, mask_file2],
                                      [x0, x1, x2],
                                      [y0, y1, y2],
                                      [bins0, bins1, bins2]):
        detector = nx.NXgroup()
        params_entry.insert(detector, name='detector'+str(i))
        if maskfile is not None:
            mask_data = fabio.open(maskfile).data
        else:
            mask_data = None
        detector.insert(nx.NXfield(mask_data), name='mask')
        if x is None:
            detector.insert(nx.NXfield(def_params['x'+str(i)]), name='beam_center_x')
        else:
            detector.insert(nx.NXfield(x), name='beam_center_x')
        if y is None:
            detector.insert(nx.NXfield(def_params['y' + str(i)]), name='beam_center_y')
        else:
            detector.insert(nx.NXfield(y), name='beam_center_y')
        if bin is None:
            detector.insert(nx.NXfield(def_params['bins' + str(i)]), name='bins')
        else:
            detector.insert(nx.NXfield(bin), name='bins')
    root.close()


@nxlib.treatment_function(new_entry=False)
def q_scale(root: nx.NXroot, distance: float = None, detector: int = 0) -> None:
    last_key = nxlib.get_last_entry_key(root)
    entry = root[last_key]
    if 'r' not in root[last_key + '/data'+str(detector)]:
        print('the last entry data' + str(detector)+'does not present r fields')
        return
    dataField = 'data'+str(detector)
    r = entry[dataField].r.nxdata
    drField = entry[dataField].r.attrs['uncertainties']
    dr = entry[dataField+'/'+drField].nxdata
    wavelength = entry.instrument.velocity_selector.wavelength
    if distance is None:
        distance = root[last_key + '/instrument/detector'+str(detector)+'/distance'].nxdata
    else:
        root[last_key + '/instrument/detector'+str(detector)+'/distance'] = distance
    signal_key = root[last_key + '/data'+str(detector)].signal
    i = entry[dataField].nxsignal
    i_errors = entry[dataField][i.attrs['uncertainties']].nxdata
    data = nx.NXdata()
    data.nxsignal = i
    data.nxerrors = i_errors
    x1 = entry.instrument.collimator.source_aperture_x
    y1 = entry.instrument.collimator.source_aperture_y
    # convert the rectangular slit to the equivalent circular slit according to
    #  Optimization of the experimental resolution for small-angle scattering
    #     D.F.R. Mildner, J.M. Carpenter
    #     J. Appl. Cryst., 1984, 17, 249-256
    r1 = np.sqrt((x1**2+y1**2)/6)
    r2 = entry.instrument.collimator.snout_diameter.nxdata/2
    l1 = entry.instrument.collimator.length.nxdata
    l2 = distance
    dlsurl = root[last_key].instrument.velocity_selector.wavelength_spread.nxdata
    q, dq = flib.qResolPinhole(r, dr, wavelength, dlsurl, r1, r2, l1, l2)
    data.nxaxes = nx.NXfield(q, name='Q', attrs={'units': r'$\AA^{-1}$', 'uncertainties': 'Q_errors'})
    data.Q_errors = nx.NXfield(dq, name='Q_errors', units=r'$\AA^{-1}$')
    del root[last_key + '/' + dataField]
    entry[dataField] = data
    return


# @nxlib.treatment_function
# def q_scale2D(root, distance=[None, None, None]):
#     entry = root[nxlib.get_last_entry_key(root)]
#     for i, d in enumerate(distance):
#         data = entry['data'+str(i)]
#         shape = data.nxsignal.nxdata.shape
#         if len(shape) == 2:
#             x_pixel_size = entry['instrument/detector' + str(i) + '/x_pixel_size'].nxdata
#             y_pixel_size = entry['instrument/detector' + str(i) + '/y_pixel_size'].nxdata
#             x0 = entry['instrument/detector' + str(i) + '/beam_center_x'].nxdata
#             y0 = entry['instrument/detector' + str(i) + '/beam_center_y'].nxdata
#             solid_angle = x_pixel_size * y_pixel_size / (d ** 2)
#             y, x = np.indices(shape, dtype='float')
#             y = (y - y0)*y_pixel_size
#             x = (x - x0)*x_pixel_size
#             r = np.sqrt(x ** 2 + y ** 2)
#             thetaf = np.arctan(x / distance)/2
#             alphaf = np.arctan(y / distance)/
#             qx = nx.NXfield(np.cos(alphaf)*np.sin(2*thetaf),name='q_x', attrs={'units'='$$'})

def save_as_txt(filename):
    root = nx.nxload(filename, mode='r')
    last_key = nxlib.get_last_entry_key(root)
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

@nxlib.treatment_function
def divide_spectra(root, denominator_file=None):
    entry = root[nxlib.get_last_entry_key(root)]
    if denominator_file is not None:
        if os.path.exists(denominator_file):
            denom_root = nx.nxload(denominator_file, mode='r')
            denom_entry = denom_root[nxlib.get_last_entry_key(denom_root)]
            for i in range(3):
                entry['data'+str(i)].nxsignal /= denom_entry['data'+str(i)].nxsignal
                signal_key = entry['data'+str(i)].signal
                if signal_key+'_errors' in entry['data'+str(i)]:
                    entry['data'+str(i)].nxerrors /= denom_entry['data'+str(i)].nxsignal

# TODO : finish this
def compute_collimation(root):
    stateList = {}
    slitList = {}
    for key in root['entry0/instrument']:
        if isinstance(root['entry0/instrument/'+key], nx.NXguide):
            stateList[key.split('_')[-1]] = root['entry0/instrument/'+key].state


# TODO : estimation of side dectectors centers using the central detector and the geometrical parameters
# TODO : estimation of the side detectros centers using calibrant (silver behenate)

@nxlib.treatment_function(new_entry=False)
def set_transmission(root, trans_file=None, direct_beam_file=None, roi=[None, None, None, None]):
    """
    compute sample trasnmission and store it in the sample.transmission field. The computation is done over the roi.
    By default this roi is 10x10 pixel centered over the center given by the scattering file in detector0 field
    Args:
        root: NXroot of the scattering file
        trans_file (str): transmited filepath
        direct_beam_file (str): direct beaml filepath
        roi (list): roi extent [X1,Y1, X2,Y2] on which is computed the transmission

    Returns:

    """
    entry = root[nxlib.get_last_entry_key(root)]
    if roi == 4*[None]:
        x0 = entry['instrument/detector0/beam_center_x'].nxdata
        y0 = entry['instrument/detector0/beam_center_y'].nxdata
        roi = [int(x0-5), int(y0-5), int(x0+5), int(y0+5)]

    if trans_file is not None and direct_beam_file is not None:
        if os.path.exists(trans_file) and os.path.exists(direct_beam_file):
            direct_root = nx.nxload(direct_beam_file, mode='r')
            trans_root = nx.nxload(trans_file, mode='r')
            crop_direct = direct_root['entry0/instrument/detector0/data'].nxdata[roi[1]:roi[3], roi[0]:roi[2]]
            crop_direct = crop_direct.sum()
            crop_trans = trans_root['entry0/instrument/detector0/data'].nxdata[roi[1]:roi[3], roi[0]:roi[2]]
            crop_trans = crop_trans.sum()
            monitor_trans = trans_root['entry0/monitor3/integral'].nxdata
            monitor_direct = direct_root['entry0/monitor3/integral'].nxdata
            t = crop_trans/crop_direct*monitor_direct/monitor_trans

            entry['sample/transmission'] = t
            print(t)
        else:
            return
    else:
        print('No transmission file or direct beam file povided')
        return


def center_of_mass_central_detector(root, roi=None):
    """
    Find the center of mass of the central detector within the roi
    Args:
        root (NXroot): NXroot
        roi (list): roi extent [X1,Y1, X2,Y2] on which the calculatio is perfomed

    Returns:

    """




if __name__ == '__main__':
    import os
    folder = '/home/achennev/python/pa20_psi/rawdatafile'
    sub_file = '/home/achennev/python/pa20_psi/sub.nxs'
    norm_file = '/home/achennev/python/pa20_psi/norm.nxs'
    samples = {
        'dark': '/home/achennev/python/pa20_psi/b4c_gq.nxs',
        'ec': '/home/achennev/python/pa20_psi/ec_gq.nxs',
        'ec_tr': '/home/achennev/python/pa20_psi/ec_tr_gq.nxs',
        'water': '/home/achennev/python/pa20_psi/h2o_gq.nxs',
        'water_tr': '/home/achennev/python/pa20_psi/h2o_tr_gq.nxs',
        'direct_beam': '/home/achennev/python/pa20_psi/direct_beam_tr_gq.nxs',
        'njc74': '/home/achennev/python/pa20_psi/njc64_73p2_gq.nxs',
        'njc74_tr': '/home/achennev/python/pa20_psi/njc64_73p2_tr_gq.nxs'
    }
    mask0 = '/home/achennev/python/pa20_psi/mask_gq.edf'

    for key in samples:
        root = nx.nxload(samples[key], mode='rw')
        nxlib.delete_all_entry(root)
        root.close()

    set_transmission(samples['water'], direct_beam_file=samples['direct_beam'], trans_file=samples['water_tr'],roi=[57, 57, 74, 74])
    set_transmission(samples['ec'], direct_beam_file=samples['direct_beam'], trans_file=samples['ec_tr'],roi=[57, 57, 74, 74])
    set_transmission(samples['njc74'], direct_beam_file=samples['direct_beam'], trans_file=samples['njc74_tr'],
                     roi=[57, 57, 74, 74])

    # reduction package
    make_reduction_package(sub_file, dark_file=samples['dark'], empty_cell_file=samples['ec'], direct_beam_file=None,
                           mask_file0=mask0, mask_file1=None,
                           x0=62.542, y0=62.824, x1=30, y1=30, bins0=100
                           )
    make_reduction_package(norm_file, dark_file=samples['dark'], empty_cell_file=samples['ec'], direct_beam_file=None,
                           water_file=samples['water'],
                           mask_file0=mask0, mask_file1=None,
                           x0=62.542, y0=62.542, x1=30, y1=30, bins0=100
                           )
    # treat_normalization_package(norm_file)
    reduction2D(samples['njc74'], sub_file=sub_file, norm_file=norm_file, thickness=0.084)
    # reduction2D(samples['water'], sub_file=sub_file, norm_file=None)
    azimutal_integration_multidetector(samples['njc74'], mask_file0=mask0)
    # azimutal_integration(samples['njc74'], mask_file=mask0, detector=0, x0=64, y0=64)
    # azimutal_integration(samples['njc74'], detector=1)
    # azimutal_integration(samples['njc74'], detector=2)
    # azimutal_integration_multidetector(njc74)
    q_scale(samples['njc74'])
    # file1 = '/home/achennev/python/pa20_psi/rawdatafile/test_nexus_AC_v2.nxs'
    # root = nxlib.loadfile(file1, mode='rw')
    # with root.nxfile:
    #     nxlib.delete_all_entry(root)
    # azimutal_integration_multidetector(file1)
    # q_scale(file1, detector=0, new_entry=True)
    # q_scale(file1, detector=1, new_entry=False)
    # q_scale(file1, detector=2, new_entry=False)

    # azimutal_integration(file1, bins=300, detector=0, mask_file=None)
    # # build_nexus_from_txt('/home/achennev/Documents/xeuss/2018-10-25-TP_SAXS/ludox1.txt')
# we have the following monitors:
# 1: monitor of proton beam
# 2: monitor before the beam shutter. This monitor serves at the moment as a standard monitor to probe the neutron flux for several other instruments, for which we do not have a monitor available. The company where we bought them in the past do not offer them anymore in the required sized. It might be removed at later a later stage.
# 3. monitor after the shutters and before the selector (white beam). This one  will stay.
# 4. monitor at the exit of the neutron guide in front of the collimation section (monochromatic beam). Neutron first pass the monitor before the go through attenuator, apertures, collimator guides etc.