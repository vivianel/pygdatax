import fabio
import nexusformat.nexus as nx
import numpy as np
from pygdatax import nxlib, flib


@nxlib.treatment_function
def set_beam_center(root, x0=None, y0=None, direct_beam_file=None):
    # root = loadfile(file, mode='rw')
    # with root.nxfile:
    entry = root[nxlib.get_last_entry_key(root)]
    if x0 is None or y0 is None:
        # read direct beam and find center around the center entered in entry
        x0 = entry.instrument.detector.beam_center_x.nxdata
        y0 = entry.instrument.detector.beam_center_y.nxdata
    if direct_beam_file is not None:
        # if direct file is given, it finds the driect beam using barycenter in a
        # a 40x40 pixel box around (x0,y0)
        root_db = nxlib.loadfile(direct_beam_file, mode='r')
        m = root_db.entry0.data.data.nxdata
        i1 = int(x0 - 40)
        i2 = int(x0 + 40)
        j1 = int(y0 - 40)
        j2 = int(y0 + 40)
        crop_m = m[j1:j2, i1:i2]
        x0, y0 = flib.find_direct_beam(m, corners=[i1, i2, j1, j2])
    entry.instrument.detector.beam_center_x = x0
    entry.instrument.detector.beam_center_y = y0
    return


@nxlib.treatment_function
def azimutal_integration(root, mask=None, x0=None, y0=None, bins=900,
                         pixel_size=None):

    last_key = nxlib.get_last_entry_key(root)
    entry = root[last_key]
    if x0 is None:
        x0 = entry.instrument.detector.beam_center_x
    else:
        entry.instrument.detector.beam_center_x = x0

    if y0 is None:
        y0 = entry.instrument.detector.beam_center_y
    else:
        entry.instrument.detector.beam_center_y = y0
    if pixel_size is None:
        pixel_size = entry.instrument.detector.x_pixel_size
    else:
        entry.instrument.detector.x_pixel_size = pixel_size

    m = entry.data.nxsignal.nxdata

    if mask is not None:
        if mask.endswith('.nxs'):
            maskroot = nxlib.loadfile(mask, mode='r')
            mask_data = maskroot.data.nxdata
        elif mask.endswith('.edf'):
            mask_data = fabio.open(mask).data
        entry.instrument.detector.pixel_mask_applied = True
        entry.instrument.detector.pixel_mask = mask_data
    else:
        mask_data = np.zeros_like(m)
        entry.instrument.detector.pixel_mask = mask_data
    if 'data_errors' in entry.data:
        print('error')
        error = entry.data.nxerrors.nxdata
    else:
        error = None
    del entry['data']
    r, i, sigma, dr = flib.regiso(m, mask_data, x0, y0, pixel_size, bins, error=error)

    # new_entry.data = nx.NXdata(attrs={'interpretation': 'spectrum',
    #                                   'signal': 'I', 'axes': ['Q']})
    entry.data = nx.NXdata()
    # new_entry.data.I = nx.NXfield(i, units='counts')  # uncertainties='I_errors'
    entry.data.nxsignal = nx.NXfield(i, name='counts', attrs={'interpretation': 'spectrum',
                                                              'uncertainties': 'counts_errors'})
    entry.data.counts_errors = sigma
    entry.data.nxaxes = nx.NXfield(r, name='r', attrs={'units': 'mm', 'uncertainties': 'r_errors'})
    entry.data.r_errors = nx.NXfield(dr, attrs={'units': 'mm'})
    # new_entry.data.nxerrors = sigma
    # new_entry.process = nx.NXprocess(program='azimutal integration',
    #                                  sequence_index=1,
    #                                  date=str(datetime.datetime.today()))
    return


@nxlib.treatment_function
def resu(root, dark_file=None, ec_file=None, eb_file=None,
         thickness=None, transmission=None, distance=None):
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

    def delta(x, a):
        if x == 1:
            y = 1
        else:
            y = (1 - x ** a) / (-a * np.log(x))
        return y

    # normalize_by_time(root.file_name, new_entry=False)  # .file_name, new_entry=False)
    i_sample = root[nxlib.get_last_entry_key(root)].data
    last_key = nxlib.get_last_entry_key(root)
    i_sample /= entry['sample/count_time']
    r = entry.data.r.nxdata

    theta = np.arctan(r / distance)

    # aT = 1 / np.cos(theta) - 1
    # aT = 1-1/np.cos(theta)

    def tr_theta(t, angle):
        if t == 1:
            t_th = 1
        else:
            t_th = t * (t ** (1 - 1 / np.cos(angle)) - 1) / (np.log(t) * (1 - 1 / np.cos(angle)))
        return t_th

    x_pixel_size = entry.instrument.detector.x_pixel_size
    y_pixel_size = entry.instrument.detector.y_pixel_size
    distance = entry.instrument.detector.distance
    solid_angle = x_pixel_size * y_pixel_size / (distance ** 2)
    flux_sample = entry.instrument.source.flux.nxdata

    # load dark file
    if dark_file is not None:
        root_dark = nxlib.loadfile(dark_file, mode='rw')
        # normalize_by_time(root_dark.file_name, new_entry=False)
        i_dark = root_dark[nxlib.get_last_entry_key(root_dark) + '/data']
        last_key = nxlib.get_last_entry_key(root_dark)
        i_dark /= root_dark[last_key + '/sample/count_time']
    else:
        root_dark = None
        # i_dark = np.zeros_like(i_sample.signal)
        i_dark = nx.NXdata(nx.NXfield(np.zeros_like(i_sample.nxsignal), name='i_dark'))
        i_dark.nxerrors = np.zeros_like(i_sample.nxsignal)

    # load emty beam
    if eb_file is not None:
        root_eb = nxlib.loadfile(eb_file, mode='rw')
        i_eb = root_eb[nxlib.get_last_entry_key(root_eb) + '/data']
        last_key = nxlib.get_last_entry_key(root_eb)
        i_eb /= root_eb[last_key + '/sample/count_time'].nxdata
        fb = i_eb - i_dark
        flux_eb = root_eb[last_key + '/instrument/source/flux'].nxdata
        fb /= flux_eb
    else:
        root_eb = None
        fb = nx.NXdata(nx.NXfield(np.zeros_like(i_sample.nxsignal), name='fb'))
        fb.nxerrors = np.zeros_like(i_sample.nxsignal)

    # load emty cell file
    if ec_file is not None:
        root_ec = nxlib.loadfile(ec_file, mode='rw')
        # normalize_by_time(root_ec.file_name, new_entry=False)
        i_ec = root_ec[nxlib.get_last_entry_key(root_ec) + '/data']
        last_key = nxlib.get_last_entry_key(root_ec)
        i_ec /= root_ec[last_key + '/sample/count_time']
        t_ec = root_ec[nxlib.get_last_entry_key(root_ec)].sample.transmission.nxdata
        t_ec = t_ec ** 0.5
        flux_ec = root_ec[last_key + '/instrument/source/flux'].nxdata
        truc = i_ec - i_dark

        z_fec = (i_ec - i_dark) / flux_ec - fb * np.power(t_ec, 2 / np.cos(theta))
        z_fec /= t_ec ** (1 / np.cos(theta)) * (t_ec + tr_theta(t_ec, theta))
    else:
        root_ec = None
        t_ec = 1
        z_fec = nx.NXdata(nx.NXfield(np.zeros_like(i_sample.nxsignal), name='z_fec'))
        z_fec.nxerrors = np.zeros_like(i_sample.nxsignal)
        # z_fec = nx.NXfield(np.zeros_like(i_sample.signal), name='z_fec')

    # substarct the contributions
    fs = (i_sample - i_dark) / flux_sample
    fs -= fb*t_ec ** (2 / np.cos(theta)) * transmission ** (1 / np.cos(theta))
    fs -= z_fec * tr_theta(t_ec, theta) * ((transmission * t_ec) ** (1 / np.cos(theta)) + t_ec * transmission)
    fs /= tr_theta(transmission, theta) * t_ec * t_ec ** (1 / np.cos(theta)) * 0.1 * thickness
    # fs = (i_sample - i_dark) / flux_sample - transmission ** (1 / np.cos(theta)) * fb*
    # fs -= transmission * delta(t_ec, aT / 2) * (1 + (transmission / t_ec ** 0.5) ** aT) * z_fec
    # fs /= 0.1 * thickness * transmission * t_ec ** (aT / 2) * delta(transmission / t_ec, aT)
    fs /= solid_angle
    fs /= np.cos(theta) ** 3

    data = nx.NXdata()
    data.nxsignal = nx.NXfield(fs.nxsignal.nxdata, name='i', attrs={'units': r'cm$^{-1}$}'})
    data.nxerrors = fs[fs.nxsignal.attrs['uncertainties']]
    data.nxaxes = fs.nxaxes
    data.r_errors = fs.r_errors
    del entry['data']
    entry['data'] = data
    q_scale(root.file_name, distance=distance, new_entry=False)
    return

# TODO : handle error propagation
@nxlib.treatment_function
def resu2D(root, dark_file=None, ec_file=None, eb_file=None,
           thickness=None, transmission=None, distance=None):
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

    # normalize_by_time(root.file_name, new_entry=False)  # .file_name, new_entry=False)
    i_sample = root['entry0'].data
    i_sample.nxerrors = nx.NXfield(np.sqrt(np.abs(i_sample.data.nxdata)))
    last_key = nxlib.get_last_entry_key(root)
    i_sample /= entry['sample/count_time']

    shape = i_sample.data.nxdata.shape
    y, x = np.indices(shape, dtype='float')
    x0 = root[nxlib.get_last_entry_key(root)].instrument.detector.beam_center_x.nxdata
    y0 = root[nxlib.get_last_entry_key(root)].instrument.detector.beam_center_y.nxdata
    y = y - y0
    x = x - x0
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan(r / distance)

    # aT = 1 / np.cos(theta) - 1
    # aT = 1-1/np.cos(theta)

    def tr_theta(t, angle):
        if t == 1:
            t_th = 1
        else:
            t_th = t * (t ** (1 - 1 / np.cos(angle)) - 1) / (np.log(t) * (1 - 1 / np.cos(angle)))
        return t_th

    x_pixel_size = entry.instrument.detector.x_pixel_size
    y_pixel_size = entry.instrument.detector.y_pixel_size
    distance = entry.instrument.detector.distance
    solid_angle = x_pixel_size * y_pixel_size / (distance ** 2)
    flux_sample = entry.instrument.source.flux.nxdata

    # load dark file
    if dark_file is not None:
        root_dark = nxlib.loadfile(dark_file, mode='rw')
        # normalize_by_time(root_dark.file_name, new_entry=False)
        i_dark = root_dark['entry0/data']
        i_dark.nxerrors = nx.NXfield(np.sqrt(np.abs(i_dark.data.nxdata)))
        last_key = nxlib.get_last_entry_key(root_dark)
        i_dark /= root_dark[last_key + '/sample/count_time']
    else:
        root_dark = None
        # i_dark = np.zeros_like(i_sample.signal)
        i_dark = nx.NXdata(nx.NXfield(np.zeros_like(i_sample.nxsignal), name='i_dark'))
        i_dark.nxerrors = np.zeros_like(i_sample.nxsignal)

    # load emty beam
    if eb_file is not None:
        root_eb = nxlib.loadfile(eb_file, mode='rw')
        i_eb = root_eb['entry0/data']
        i_eb.nxerrors = nx.NXfield(np.sqrt(np.abs(i_eb.data.nxdata)))
        last_key = nxlib.get_last_entry_key(root_eb)
        i_eb /= root_eb[last_key + '/sample/count_time'].nxdata
        fb = i_eb - i_dark
        flux_eb = root_eb[last_key + '/instrument/source/flux'].nxdata
        fb /= flux_eb
    else:
        root_eb = None
        fb = nx.NXdata(nx.NXfield(np.zeros_like(i_sample.nxsignal), name='fb'))
        fb.nxerrors = np.zeros_like(i_sample.nxsignal)

    # load emty cell file
    if ec_file is not None:
        root_ec = nxlib.loadfile(ec_file, mode='rw')
        # normalize_by_time(root_ec.file_name, new_entry=False)
        i_ec = root_ec['entry0/data']
        i_ec.nxerrors = nx.NXfield(np.sqrt(np.abs(i_ec.data.nxdata)))
        last_key = nxlib.get_last_entry_key(root_ec)
        i_ec /= root_ec[last_key + '/sample/count_time']
        t_ec = root_ec[nxlib.get_last_entry_key(root_ec)].sample.transmission.nxdata
        t_ec = t_ec ** 0.5
        flux_ec = root_ec[last_key + '/instrument/source/flux'].nxdata
        z_fec = (i_ec - i_dark) / flux_ec - fb * np.power(t_ec, 2 / np.cos(theta))
        z_fec /= t_ec ** (1 / np.cos(theta)) * (t_ec + tr_theta(t_ec, theta))
    else:
        root_ec = None
        t_ec = 1
        z_fec = nx.NXdata(nx.NXfield(np.zeros_like(i_sample.nxsignal), name='z_fec'))
        z_fec.nxerrors = np.zeros_like(i_sample.nxsignal)
        # z_fec = nx.NXfield(np.zeros_like(i_sample.signal), name='z_fec')

    # substarct the contributions
    fs = (i_sample - i_dark) / flux_sample
    fs -= fb*t_ec ** (2 / np.cos(theta)) * transmission ** (1 / np.cos(theta))
    fs -= z_fec * tr_theta(t_ec, theta) * ((transmission * t_ec) ** (1 / np.cos(theta)) + t_ec * transmission)
    fs /= tr_theta(transmission, theta) * t_ec * t_ec ** (1 / np.cos(theta)) * 0.1 * thickness
    # fs = (i_sample - i_dark) / flux_sample - transmission ** (1 / np.cos(theta)) * fb*
    # fs -= transmission * delta(t_ec, aT / 2) * (1 + (transmission / t_ec ** 0.5) ** aT) * z_fec
    # fs /= 0.1 * thickness * transmission * t_ec ** (aT / 2) * delta(transmission / t_ec, aT)
    fs /= solid_angle
    fs /= np.cos(theta) ** 3

    data = nx.NXdata()
    data.nxsignal = nx.NXfield(fs.nxsignal.nxdata, name='i', attrs={'units': r'cm$^{-1}$}'})
    # data.nxerrors = fs[fs.nxsignal.attrs['uncertainties']]
    # data.nxaxes = fs.nxaxes
    # data.r_errors = fs.r_errors
    del entry['data']
    entry['data'] = fs
    # q_scale(root.file_name, distance=distance, new_entry=False)
    return

# @treatment_function
# def divide_by_flux(root, flux=None):
#     # root = loadfile(file, mode='rw')
#     last_key = get_last_entry_key(root)
#     entry = root[last_key]
#     if flux is not None:
#         entry.instrument.source.flux.nxdata = flux
#     else:
#         flux = entry.instrument.source.flux.nxdata
#     if isinstance(entry.data.nxsignal, nx.tree.NXlinkfield):
#         signal_key = entry.data.nxsignal.nxtarget
#         root[signal_key] /= flux
#     else:
#         signal_key = last_key + '/data/' + entry.data.signal
#         root[signal_key] /= flux
#         if 'errors' in root[last_key + '/data']:
#             root[last_key + '/data/errors'] /= flux
#     return


# @treatment_function
# def normalize_by_time(root, count_time=None):
#     # root = loadfile(file, mode='rw')
#     last_key = get_last_entry_key(root)
#     entry = root[last_key]
#     if count_time is None:
#         count_time = entry.sample.count_time.nxdata
#     else:
#         entry.sample.count_time.nxdata = count_time
#     if isinstance(entry.data.nxsignal, nx.tree.NXlinkfield):
#         signal_key = entry.data.nxsignal.nxtarget
#         root[signal_key] /= count_time
#     else:
#         signal_key = last_key + '/data/' + entry.data.signal
#         root[signal_key] /= count_time
#         if 'errors' in root[last_key + '/data']:
#             root[last_key + '/data/errors'] /= count_time
#     return
#
#
# @treatment_function
# def divide_by_transmission(root, transmission=None, new_entry=True):
#     last_key = get_last_entry_key(root)
#     entry = root[last_key]
#     if transmission is not None:
#         entry.sample.transmission.nxdata = transmission
#     else:
#         transmission = entry.sample.transmission.nxdata
#     if isinstance(entry.data.nxsignal, nx.tree.NXlinkfield):
#         signal_key = entry.data.nxsignal.nxtarget
#     else:
#         signal_key = last_key + '/data/' + entry.data.signal
#         if 'errors' in root[last_key + '/data']:
#             root[last_key + '/data/errors'] /= transmission
#     root[signal_key] /= transmission
#     return
#
#
# @treatment_function
# def divide_by_thickness(root, thickness=None):
#     # divide by thickness concerted in cm
#     # root = loadfile(file, mode='rw')
#     last_key = get_last_entry_key(root)
#     entry = root[last_key]
#     if thickness is not None:
#         entry.sample.thickness.nxdata = thickness
#     else:
#         thickness = entry.sample.thickness.nxdata
#     if isinstance(entry.data.nxsignal, nx.tree.NXlinkfield):
#         signal_key = entry.data.nxsignal.nxtarget
#     else:
#         signal_key = last_key + '/data/' + entry.data.signal
#         if 'errors' in root[last_key + '/data']:
#             root[last_key + '/data/errors'] /= (thickness * 10)
#     root[signal_key] /= (thickness * 10)
#     return


@nxlib.treatment_function
def q_scale(root, distance=None):
    last_key = nxlib.get_last_entry_key(root)
    entry = root[last_key]
    if 'r' not in root[last_key + '/data']:
        print('the last entry data does not present r fields')
        return
    r = entry.data.r.nxdata
    dr = entry.data[entry.data.r.attrs['uncertainties']].nxdata
    wavelength = entry.instrument.source.incident_wavelength
    if distance is None:
        distance = entry.instrument.detector.distance.nxdata
    else:
        entry.instrument.detector.distance = nx.NXfield(distance, attrs={'units' : 'mm'})
    signal_key = root[last_key + '/data'].signal
    i = entry.data.nxsignal
    i_errors = entry.data[i.attrs['uncertainties']].nxdata
    data = nx.NXdata()
    data.nxsignal = i
    data.nxerrors = i_errors
    x1 = entry.instrument.source.beam_size_x.nxdata/2
    y1 = entry.instrument.source.beam_size_y.nxdata / 2
    x2 = entry.instrument.slit.x_gap.nxdata / 2
    y2 = entry.instrument.slit.y_gap.nxdata / 2
    l1 = entry.instrument.collimator.length.nxdata
    l2 = distance
    l = root[last_key].instrument.collimator.distance.nxdata
    dlsurl = root[last_key].instrument.source.incident_wavelength_spread.nxdata
    q, dq = flib.qResolSlits(r, dr, wavelength, dlsurl, x1, y1, x2, y2, l1, l2, l)
    data.nxaxes = nx.NXfield(q, name='Q', attrs={'units': r'$\AA^{-1}$', 'uncertainties': 'Q_errors'})
    data.Q_errors = nx.NXfield(dq, name='Q_errors', attrs={'units': r'$\AA^{-1}$'})
    del root[last_key + '/data']
    entry.data = data
    return


@nxlib.treatment_function
def azimutal_integration2D(root, mask=None, x0=None, y0=None, distance=None,
                           r_bins=900, chi_bins=360, pixel_size=None):
    last_key = nxlib.get_last_entry_key(root)
    entry = root[last_key]
    if x0 is None:
        x0 = entry.instrument.detector.beam_center_x.nxdata
    else:
        entry.instrument.detector.beam_center_x = x0
    if y0 is None:
        y0 = entry.instrument.detector.beam_center_y.nxdata
    else:
        entry.instrument.detector.beam_center_x = y0
    if distance is None:
        distance = entry.instrument.detector.distance
    else:
        entry.instrument.detector.distance = distance
    if pixel_size is None:
        pixel_size = entry.instrument.detector.x_pixel_size
    else:
        entry.instrument.detector.x_pixel_size = pixel_size
    m = entry.data.data.nxdata
    mask_data = np.zeros_like(m)
    if mask is not None:
        if mask.endswith('.nxs'):
            maskroot = nxlib.loadfile(mask, mode='r')
            mask_data = maskroot.data.nxdata
        elif mask.endswith('.edf'):
            mask_data = fabio.open(mask).data
        entry.instrument.detector.pixel_mask_applied = True
    entry.instrument.detector.pixel_mask = mask_data
    del entry['data']
    r_grid, chi_grid, masked_data = flib.xy2polar(m, mask_data, x0, y0, r_bins, chi_bins)
    wavelength = entry.instrument.source.incident_wavelength
    q = 4*np.pi*np.sin(np.arctan(r_grid*pixel_size/distance)/2)/wavelength

    # new_entry.data = nx.NXdata(attrs={'interpretation': 'spectrum',
    #                                   'signal': 'I', 'axes': ['Q']})
    y = nx.NXfield(q, name='Q', attrs={'units': r'$\AA^{-1}$'})
    x = nx.NXfield(chi_grid, name='chi', attrs={'units': 'deg'})
    entry.data = nx.NXdata(nx.NXfield(masked_data, name='counts'), [y, x])
    entry.data.attrs['interpretation'] = 'image'
    return


@nxlib.treatment_function
def polar_cut(root, q=None, pixel_width=1):
    last_key = nxlib.get_last_entry_key(root)
    entry = root[last_key]
    x = entry.data.Q.nxdata
    if q is None:
        q = x[0]
    chi = entry.data.chi.nxdata
    counts = entry.data.counts.nxdata
    closest = np.min(np.abs(q-x))
    k = np.argwhere(np.abs(q-x) == closest)[0][0]
    k -= int(pixel_width/2)
    index = []
    for p in range(pixel_width):
        index.append(k)
        k = k+1
    i = counts[index, :]
    i = np.ma.masked_less(i, 0)
    # i_errors = np.sqrt(np.sum(i,axis=0)).compressed()
    i = np.mean(i, axis=0)
    chi = np.ma.masked_array(chi, mask=i.mask).compressed()
    i = i.compressed()
    del entry['data']
    x = nx.NXfield(chi, name='chi', attrs={'units': 'deg'})
    entry.data = nx.NXdata(nx.NXfield(i, name='counts'), [x])
    # entry.data.errors=nx.NXfield(i_errors, name='errors')
    entry.data.q_cut = nx.NXfield(q, name='q_cut', attrs={'units': r'$\AA^{-1}$'})
    entry.data.line_number = nx.NXfield(pixel_width, name='line_number')
    return


@nxlib.treatment_function
def normalization_factor(root, factor=None):
    last_key = nxlib.get_last_entry_key(root)
    entry = root[last_key]
    entry.data.nxsignal *= factor
    entry.data.nxerrors *= factor
    return


@nxlib.treatment_function
def bkg_substraction(root, bkg=0):
    """
    Substract conctant the the nxsignal
    Args:
        root (NXRoot): NXroot
        bkg (float): Constant to be substracted

    Returns:

    """
    last_key = nxlib.get_last_entry_key(root)
    entry = root[last_key]
    if np.ndim(entry.data.nxsignal.nxdata) == 1:
        entry.data.nxsignal -= bkg
        return


@nxlib.treatment_function
def ref_substraction(root, ref_file=None, prefactor=1):
    """
    Substract reference file (i.e. solvant) from the data wiegthd by a prefector
    I = I_sample - I_ref * prefactor
    Args:
        root (NXroot): NXroot
        ref_file (str): file to substract
        prefactor (float): prefactor

    Returns:

    """
    last_key = nxlib.get_last_entry_key(root)
    entry = root[last_key]
    if ref_file:
        root_ref = nxlib.loadfile(ref_file, mode='r')
        date_ref = root_ref[nxlib.get_last_entry_key(root_ref)].data
        if np.ndim(entry.data.nxsignal.nxdata) == 1:
            entry.data.nxsignal -= date_ref.nxsignal*prefactor
    return


@nxlib.treatment_function
def cut(root, xmin=None, xmax=None):
    last_key = nxlib.get_last_entry_key(root)
    entry = root[last_key]
    if np.ndim(entry.data.nxsignal.nxdata) == 1:
        x = entry.data.nxaxes[0].nxdata
        y = entry.data
        if xmin is None:
            xmin = np.min(x)
        if xmax is None:
            xmax = np.max(x)
        index1 = np.argwhere(x >= xmin)
        index2 = np.argwhere(x <= xmax)
        index = np.intersect1d(index1, index2)
        data = nx.NXdata()
        i = entry.data.nxsignal[index]
        i_errors = entry.data[i.attrs['uncertainties']].nxdata[index]
        data.nxsignal = i
        data.nxerrors = i_errors
        data.nxaxes = entry.data.nxaxes[0][index]
        x_errors_key = entry.data.nxaxes[0].attrs['uncertainties']
        data.insert(entry.data[x_errors_key][index], name=x_errors_key)
        del entry['data']
        entry.data = data
        return


@nxlib.treatment_function
def concat(root, file=None):
    last_key = nxlib.get_last_entry_key(root)
    entry = root[last_key]
    new_root = nxlib.loadfile(file, mode='r')
    new_entry = new_root[nxlib.get_last_entry_key(new_root)]
    # check if nxdata fields are identical
    keys = set(entry.data.keys())
    new_keys = set(new_entry.data.keys())
    if keys == new_keys and entry.data.nxsignal.nxdata.ndim == 1:
        i1 = entry.data.nxsignal.nxdata
        if 'uncertainties' in entry.data.nxsignal.attrs:
            i1_errors = entry.data[entry.data.nxsignal.attrs['uncertainties']].nxdata
        else:
            i1_errors = np.zeros_like(i1)
        if 'units' in entry.data.nxsignal.attrs:
            i_units = entry.data.nxsignal.attrs['units']
        else:
            i_units = None
        x1_key = entry.data.axes[0]
        x1 = entry.data[x1_key].nxdata
        if 'uncertainties' in entry.data[x1_key].attrs:
            x1_errors = entry.data[entry.data[x1_key].attrs['uncertainties']].nxdata
        else:
            x1_errors = np.zeros_like(x1)
        if 'units' in entry.data[x1_key].attrs:
            x_units = entry.data[x1_key].attrs['units']
        else:
            x_units = None

        i2 = new_entry.data.nxsignal.nxdata
        if 'uncertainties' in new_entry.data.nxsignal.attrs:
            i2_errors = new_entry.data[new_entry.data.nxsignal.attrs['uncertainties']].nxdata
        else:
            i2_errors = np.zeros_like(i2)
        x2_key = new_entry.data.axes[0]
        x2 = new_entry.data[x2_key].nxdata
        if 'uncertainties' in new_entry.data[x2_key].attrs:
            x2_errors = new_entry.data[new_entry.data[x2_key].attrs['uncertainties']].nxdata
        else:
            x2_errors = np.zeros_like(x2)

        matini = np.array([np.append(x1, x2),
                           np.append(i1, i2),
                           np.append(i1_errors, i2_errors),
                           np.append(x1_errors, x2_errors)])
        mat = matini[:, np.argsort(matini[0, :])]
        data = nx.NXdata()
        data.nxsignal = nx.NXfield(mat[1, :],
                                   name=entry.data.attrs['signal'],
                                   units=i_units)
        data.nxerrors = nx.NXfield(mat[2, :],
                                   name=entry.data.attrs['signal'] + '_errors')
        data.nxaxes = nx.NXfield(mat[0, :],
                                 name=entry.data.attrs['axes'][0],
                                 units=x_units)
        data.insert(mat[3, :],
                    name=entry.data.attrs['axes'][0] + '_errors')
        # data.nxsignal = nx.NXfield(np.append(i1, i2),
        #                            name=entry.data.attrs['signal'],
        #                            units=i_units)
        # data.nxerrors = nx.NXfield(np.append(i1_errors, i2_errors),
        #                            name=entry.data.attrs['signal']+'_errors')
        # data.nxaxes = nx.NXfield(np.append(x1, x2),
        #                          name=entry.data.attrs['axes'][0],
        #                          units=x_units)
        # data.insert(np.append(x1_errors, x2_errors),
        #             name=entry.data.attrs['axes'][0]+'_errors')
        del entry['data']
        entry.data = data
        return


def save_as_txt(filename):
    if not filename.endswith('.nxs'):
        print('cannot convert non nexus file')
        return
    root = nx.nxload(filename, mode='r')
    last_key = nxlib.get_last_entry_key(root)
    try:
        sample_description = root[last_key].sample.sample_name.nxdata.decode()
    except AttributeError:
        sample_description = root[last_key].sample.sample_name.nxdata
    data = root[last_key + '/data']
    signal_key = data.signal
    if isinstance(data[signal_key], nx.tree.NXlinkfield):
        signal_key = data.nxsignal.nxtarget
    else:
        signal_key = last_key + '/data/' + data.signal
    y = root[signal_key].nxdata
    dy = data.nxerrors
    if len(y.shape) > 1:
        print("Cannot save 2D data as .txt file")
        return

    x = data.nxaxes[0].nxdata
    dx = None
    if data.axes+'_errors' in data:
        dx = data[data.axes+'_errors'].nxdata
    if dx is not None and dy is not None:
        mat = np.stack((x, y, dy, dx))
    elif dx is None and dy is not None:
        mat = np.stack((x, y, dy))
    else:
        mat = np.stack((x, y))
    newfile = filename.split('.')[0]+'_'+sample_description+'.txt'
    np.savetxt(newfile, mat.transpose(), delimiter='\t')


if __name__ == '__main__':
    # build_nexus_from_edf('AgBe.edf')
    # set_beam_center('AgBe.nxs', x0=200, y0=300, new_entry=False)
    # file = 'AgBe.edf'
    # empty_cell_GQ = '/home/achennev/Documents/xeuss/fc aout 2020/kapton_gq.edf'
    # dark_GQ = '/home/achennev/Documents/xeuss/fc aout 2020/dark_GQ.edf'
    # file1_GQ = '/home/achennev/Documents/xeuss/fc aout 2020/147_1.edf'
    # file2_GQ = '/home/achennev/Documents/xeuss/fc aout 2020/147_2.edf'
    # mask = '/home/achennev/Documents/xeuss/fc aout 2020/mask_silx_GQ.nxs'
    # nxlib.build_nexus_from_txt('/home/achennev/Documents/xeuss/fc aout 2020/147_1_pasi.txt')
    # nxlib.build_nexus_from_txt('/home/achennev/Documents/xeuss/fc aout 2020/147_2_pasi.txt')
    # nxlib.build_nexus_from_edf(file1_GQ)
    # nxlib.build_nexus_from_edf(file2_GQ)
    #
    # file1_GQ = file1_GQ.split('.')[0] + '.nxs'
    # file2_GQ = file2_GQ.split('.')[0] + '.nxs'
    # # set_beam_center(empty_cell, x0=x0, y0=y0, new_entry=False)# direct_beam_file=directbeam, new_entry=False)
    # azimutal_integration(empty_cell_GQ, bins=900, mask=mask)
    # # # set_beam_center(dark, x0=x0, y0=y0, new_entry=False)#direct_beam_file=directbeam, new_entry=False)
    # azimutal_integration(dark_GQ, bins=900, mask=mask)
    # dark_GQ = dark_GQ.split('.')[0] + '.nxs'
    # empty_cell_GQ = empty_cell_GQ.split('.')[0] + '.nxs'

    # azimutal_integration(file1_GQ, bins=900, mask=mask)
    # azimutal_integration(file2_GQ, bins=900, mask=mask)
    #
    # resu(file1_GQ, dark_file=dark_GQ, ec_file=empty_cell_GQ, thickness=10)
    # resu(file2_GQ, dark_file=dark_GQ, ec_file=empty_cell_GQ, thickness=10)

    # long distance
    empty_cell_PQ = '/home/achennev/Documents/xeuss/fc aout 2020/kapton_PQ.edf'
    dark_PQ = '/home/achennev/Documents/xeuss/fc aout 2020/dark_GQ.edf'
    file1_PQ = '/home/achennev/Documents/xeuss/fc aout 2020/147_1_PQ.edf'
    file2_PQ = '/home/achennev/Documents/xeuss/fc aout 2020/147_2_PQ.edf'
    mask_PQ = '/home/achennev/Documents/xeuss/fc aout 2020/mask_silx_PQ.nxs'
    nxlib.build_nexus_from_txt('/home/achennev/Documents/xeuss/fc aout 2020/147_1_PQ_pasi.txt')
    nxlib.build_nexus_from_txt('/home/achennev/Documents/xeuss/fc aout 2020/147_2_PQ_pasi.txt')
    nxlib.build_nexus_from_edf(file1_PQ)
    nxlib.build_nexus_from_edf(file2_PQ)
    nxlib.build_nexus_from_edf(dark_PQ)
    nxlib.build_nexus_from_edf(empty_cell_PQ)

    file1_PQ = file1_PQ.split('.')[0] + '.nxs'
    file2_PQ = file2_PQ.split('.')[0] + '.nxs'
    dark_PQ = dark_PQ.split('.')[0] + '.nxs'
    file1_PQ_2D = '/home/achennev/Documents/xeuss/fc aout 2020/147_1_PQ _2D.nxs'
    root2D = nxlib.loadfile(file1_PQ_2D)
    nxlib.delete_all_entry(root2D)
    root2D.close()
    empty_cell_PQ = empty_cell_PQ.split('.')[0] + '.nxs'

    print(2)
    azimutal_integration(empty_cell_PQ, bins=900, mask=mask_PQ)
    print(3)
    azimutal_integration(dark_PQ, bins=900, mask=mask_PQ)
    print(4)
    azimutal_integration(file1_PQ, bins=900, mask=mask_PQ)
    print(5)
    azimutal_integration(file2_PQ, bins=900, mask=mask_PQ)

    resu(file1_PQ, dark_file=dark_PQ, ec_file=empty_cell_PQ, thickness=1)
    resu2D(file1_PQ_2D, dark_file=dark_PQ, ec_file=empty_cell_PQ)
    print(1)
    azimutal_integration(file1_PQ_2D, bins=900, mask=mask_PQ)
    q_scale(file1_PQ_2D)
    # resu(file2_PQ, dark_file=dark_PQ, ec_file=empty_cell_PQ, thickness=10)
    # # normalization_factor(file1_PQ, factor=10)
    # # cut(file1_PQ, xmin=None, xmax=0.1)
    # # concat(file1_PQ, file=file1_GQ)
    # # # build_nexus_from_txt('/home/achennev/Documents/xeuss/2018-10-25-TP_SAXS/ludox1.txt')
    #
    # # test 2D resu

