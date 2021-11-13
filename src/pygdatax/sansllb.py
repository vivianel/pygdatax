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
import flib
import nxlib

NXREAD_VERSION = '0.0'


@nxlib.treatment_function
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


@nxlib.treatment_function
def azimutal_integration(root, detector=0, mask_file=None, x0=None, y0=None, bins=90,
                         pixel_size=None):
    last_key = nxlib.get_last_entry_key(root)
    entry = root[last_key]
    if x0 is None:
        x0 = entry['instrument/detector'+str(detector)+'/beam_center_x']
    else:
        entry['instrument/detector'+str(detector)+'/beam_center_x'] = x0

    if y0 is None:
        y0 = entry['instrument/detector'+str(detector)+'/beam_center_y']
    else:
        entry['instrument/detector'+str(detector)+'/beam_center_y'] = y0
    if pixel_size is None:
        pixel_size = entry['instrument/detector'+str(detector)+'/x_pixel_size']
    else:
        entry['instrument/detector'+str(detector)+'/x_pixel_size'] = pixel_size

    m = entry['data'+str(detector)+'/data'].nxdata
    mask_data = np.zeros_like(m)

    if mask_file is not None:
        if mask_file.endswith('.nxs'):
            maskroot = nxlib.loadfile(mask_file, mode='r')
            mask_data = maskroot.data.nxdata
        elif mask_file.endswith('.edf'):
            mask_data = fabio.open(mask_file).data
    del entry['data' + str(detector)]
    entry['instrument/detector' + str(detector) + '/pixel_mask_applied'] = True
    entry['instrument/detector' + str(detector) + '/pixel_mask'] = mask_data
    r, i, sigma, dr = flib.regiso(m, mask_data, x0, y0, pixel_size, bins)

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
def azimutal_integration_multidetector(root, mask_file0=None, x0=None, y0=None, bins0=90, pixel_size0=None,
                                       mask_file1=None, x1=None, y1=None, bins1=90, pixel_size1=None,
                                       mask_file2=None, x2=None, y2=None, bins2=90, pixel_size2=None
                                       ):
    azimutal_integration(root.file_name, detector=0, mask_file=mask_file0, x0=x0, y0=y0, bins=bins0, pixel_size=pixel_size0,
                         new_entry=False)
    azimutal_integration(root.file_name, detector=1, mask_file=mask_file1, x0=x1, y0=y1, bins=bins1, pixel_size=pixel_size1,
                         new_entry=False)
    azimutal_integration(root.file_name, detector=2, mask_file=mask_file2, x0=x2, y0=y2, bins=bins2, pixel_size=pixel_size2,
                         new_entry=False)


# TODO: finish this function with monitor normalization
@nxlib.treatment_function
def reduction2D(root: nx.NXroot, sub_file=None, norm_file=None,
                thickness: float = None, transmission: float = None, distance: float = None) -> None:
    """

    Args:
        root:
        sub_file:
        norm_file:
        thickness:
        transmission:
        distance:

    Returns:

    """
    entry = root[nxlib.get_last_entry_key(root)]
    if distance is None:
        distance = entry.instrument.detector.distance.nxdata
    else:
        entry.instrument.detector.distance.nxdata = distance
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

    i_sample = root['entry0'].data
    i_sample.nxerrors = nx.NXfield(np.sqrt(np.abs(i_sample.data.nxdata)))
    monitor_sample = root['entry0/monitor3/integral'].nxdata
    time_sample = root['entry0/instrument/detector0/count_time'].nxdata
    shape = i_sample.data.nxdata.shape

    # uncack the substraction package file
    default_params = get_default_reduction_parameters(root)
    if sub_file is not None:
        sub_root = nx.nxload(sub_file, mode='r')
        if 'dark' in sub_root:
            i_dark = sub_root['dark/data0']
            i_dark.nxerrors = nx.NXfield(np.sqrt(np.abs(i_dark.data.nxdata)))
            monitor_dark = sub_root['dark/monitor3/integral'].nxdata
            time_dark = sub_root['dark/instrument/detector0/count_time'].nxdata
        else:
            i_dark = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_dark'))
            i_dark.nxerrors = np.zeros(shape)
            monitor_dark = 1
            time_dark = 1
        if 'empty cell' in sub_root:
            i_ec = sub_root['empty cell/data0']
            i_ec.nxerrors = nx.NXfield(np.sqrt(np.abs(i_ec.data.nxdata)))
            monitor_ec = sub_root['empty cell/monitor3/integral'].nxdata
            time_ec = sub_root['empty cell/instrument/detector0/count_time'].nxdata
            trans_ec = sub_root['empty cell/sample/transmission'].nxdata
        else:
            i_ec = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_ec'))
            i_ec.nxerrors = np.zeros(shape)
            monitor_ec = 1
            time_ec = 1
            trans_ec = 1
        if 'empty beam' in sub_root:
            i_eb = sub_root['empty beam/data0']
            i_eb.nxerrors = nx.NXfield(np.sqrt(np.abs(i_eb.data.nxdata)))
            monitor_eb = sub_root['empty beam/monitor3/integral'].nxdata
            time_eb = sub_root['empty beam/instrument/detector0/count_time'].nxdata
        else:
            i_eb = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_eb'))
            i_eb.nxerrors = np.zeros(shape)
            monitor_eb = 1
            time_eb = 1

        # centers for regroupement
        if 'beam_center_x' in sub_root['parameters/detector0']:
            x0 = sub_root['parameters/detector0/beam_center_x']
        else:
            x0 = default_params['x0']
        if 'beam_center_y' in sub_root['parameters/detector0']:
            y0 = sub_root['parameters/detector0/beam_center_y']
        else:
            y0 = default_params['y0']

    else:
        i_dark = nx.NXdata(nx.NXfield(np.zeros(shape), name='i_dark'))
        i_dark.nxerrors = np.zeros(shape)
        monitor_dark = 1
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

    y, x = np.indices(shape, dtype='float')
    y = y - y0
    x = x - x0
    r = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan(r / distance)

    x_pixel_size = entry.instrument.detector.x_pixel_size
    y_pixel_size = entry.instrument.detector.y_pixel_size
    distance = entry.instrument.detector.distance
    solid_angle = x_pixel_size * y_pixel_size / (distance ** 2)
    flux_sample = entry.instrument.source.flux.nxdata

    fb = i_eb - i_dark*time_eb/time_dark
    fb /= monitor_eb
    trans_ec = trans_ec ** 0.5
    if monitor_ec == 1 and time_ec == 1 and trans_ec == 1: # no empty cell given

    else:
        z_fec = (i_ec - i_dark*time_ec/time_dark) / monitor_ec - fb * np.power(trans_ec, 2 / np.cos(theta))
        z_fec /= trans_ec ** (1 / np.cos(theta)) * (trans_ec + tr_theta(trans_ec, theta))

    #
    # # load emty beam
    # if eb_file is not None:
    #     root_eb = nxlib.loadfile(eb_file, mode='rw')
    #     i_eb = root_eb['entry0/data']
    #     i_eb.nxerrors = nx.NXfield(np.sqrt(np.abs(i_eb.data.nxdata)))
    #     last_key = nxlib.get_last_entry_key(root_eb)
    #     i_eb /= root_eb[last_key + '/sample/count_time'].nxdata
    #     fb = i_eb - i_dark
    #     flux_eb = root_eb[last_key + '/instrument/source/flux'].nxdata
    #     fb /= flux_eb
    # else:
    #     root_eb = None
    #     fb = nx.NXdata(nx.NXfield(np.zeros_like(i_sample.nxsignal), name='fb'))
    #     fb.nxerrors = np.zeros_like(i_sample.nxsignal)
    #
    # # load emty cell file
    # if ec_file is not None:
    #     root_ec = nxlib.loadfile(ec_file, mode='rw')
    #     # normalize_by_time(root_ec.file_name, new_entry=False)
    #     i_ec = root_ec['entry0/data']
    #     i_ec.nxerrors = nx.NXfield(np.sqrt(np.abs(i_ec.data.nxdata)))
    #     last_key = nxlib.get_last_entry_key(root_ec)
    #     i_ec /= root_ec[last_key + '/sample/count_time']
    #     t_ec = root_ec[nxlib.get_last_entry_key(root_ec)].sample.transmission.nxdata
    #     t_ec = t_ec ** 0.5
    #     flux_ec = root_ec[last_key + '/instrument/source/flux'].nxdata
    #     z_fec = (i_ec - i_dark) / flux_ec - fb * np.power(t_ec, 2 / np.cos(theta))
    #     z_fec /= t_ec ** (1 / np.cos(theta)) * (t_ec + tr_theta(t_ec, theta))
    # else:
    #     root_ec = None
    #     t_ec = 1
    #     z_fec = nx.NXdata(nx.NXfield(np.zeros_like(i_sample.nxsignal), name='z_fec'))
    #     z_fec.nxerrors = np.zeros_like(i_sample.nxsignal)
    #     # z_fec = nx.NXfield(np.zeros_like(i_sample.signal), name='z_fec')
    #
    # # substarct the contributions
    # fs = (i_sample - i_dark) / flux_sample
    # fs -= fb * t_ec ** (2 / np.cos(theta)) * transmission ** (1 / np.cos(theta))
    # fs -= z_fec * tr_theta(t_ec, theta) * ((transmission * t_ec) ** (1 / np.cos(theta)) + t_ec * transmission)
    # fs /= tr_theta(transmission, theta) * t_ec * t_ec ** (1 / np.cos(theta)) * 0.1 * thickness
    # # fs = (i_sample - i_dark) / flux_sample - transmission ** (1 / np.cos(theta)) * fb*
    # # fs -= transmission * delta(t_ec, aT / 2) * (1 + (transmission / t_ec ** 0.5) ** aT) * z_fec
    # # fs /= 0.1 * thickness * transmission * t_ec ** (aT / 2) * delta(transmission / t_ec, aT)
    # fs /= solid_angle
    # fs /= np.cos(theta) ** 3
    #
    # data = nx.NXdata()
    # data.nxsignal = nx.NXfield(fs.nxsignal.nxdata, name='i', attrs={'units': r'cm$^{-1}$}'})
    # # data.nxerrors = fs[fs.nxsignal.attrs['uncertainties']]
    # # data.nxaxes = fs.nxaxes
    # # data.r_errors = fs.r_errors
    # del entry['data']
    # entry['data'] = fs
    # q_scale(root.file_name, distance=distance, new_entry=False)
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
    y2 = root.entry.instrument.detector2.beam_center_y.nxdata
    bins2 = 16
    dic = {'x0': x0, 'y0': y0, 'bins0': bins0,
           'x1': x1, 'y1': y1, 'bins1': bins1,
           'x2': x2, 'y2': y2, 'bins2': bins2
           }
    return dic


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
        dark_entry = dark_root['entry0'].copy()
        root.insert(dark_entry, name='dark')
        dark_root.close()
        def_params = get_default_reduction_parameters(dark_root)

    if empty_cell_file is not None:
        ec_root = nx.nxload(empty_cell_file, mode='r')
        ec_entry = ec_root['entry0'].copy()
        root.insert(ec_entry, name='empty cell')
        ec_root.close()
        def_params = get_default_reduction_parameters(ec_root)

    if direct_beam_file is not None:
        db_root = nx.nxload(direct_beam_file, mode='r')
        db_entry = db_root['entry0'].copy()
        root.insert(db_entry, name='empty beam')
        db_root.close()
        def_params = get_default_reduction_parameters(db_root)

    if water_file is not None:
        water_root = nx.nxload(water_file, mode='r')
        water_entry = water_root['entry0'].copy()
        root.insert(water_entry, name='water')
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


@nxlib.treatment_function
def q_scale(root, distance=None, detector=0):
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


if __name__ == '__main__':
    import os
    folder = '/home/achennev/Bureau/PIL pour tiago/PIL NP/2021-04-07-TC'
    file1 = os.path.join(folder, '2021-04-07-TC_0_87156.nxs')
    file2 = os.path.join(folder, '2021-04-07-TC_0_87156.nxs')
    maskfile1 = os.path.join(folder, 'mask_WAXS.edf')
    maskfile2 = os.path.join(folder, 'mask_MAXS.edf')
    output_file = os.path.join(folder, 'first_package.nxs')
    file = '/home/achennev/python/pa20_psi/rawdatafile/test_nexus_AC_v2.nxs'
    output_file = '/home/achennev/python/pa20_psi/rawdatafile/test_sub.nxs'
    make_reduction_package(output_file, dark_file=file, empty_cell_file=file,
                           mask_file0=maskfile1, mask_file1=maskfile2,
                           x0=1, y0=2, x1=3, y1=3, bins0=100
                           )
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