import numpy as np
import nexusformat.nexus as nx
import time

"""
Python Script creating nexus file architecture for SANS_LLB using nexusformat 
python API.
"""

# nexus root definition
root = nx.NXroot()
root.attrs['default'] = 'entry0'
root.attrs['NX_class'] = 'NXroot'
entry0 = nx.NXentry(attrs={'default': 'data0'})
root.entry0 = entry0
# start and end time of the acquisition
entry0.start_time = nx.NXfield(time.asctime())
entry0.end_time = nx.NXfield(time.asctime())
# Defnition of the instrument field
instrument = nx.NXinstrument()
instrument.name = "SANS_LLB"
entry0.instrument = instrument

# velocitity selector
# entry0/instrument/velocity_selector
vs = nx.NXvelocity_selector()
vs.wavelength = nx.NXfield(5.5, attrs={'units': 'Angstrom'})
vs.wavelength_set = nx.NXfield(5.4, attrs={'units': 'Angstrom'})
vs.wavelength_spread = nx.NXfield(0.11)
vs.type = nx.NXfield('NVS model/type.')
vs.rotation_speed = nx.NXfield(3600.0, attrs={'units': 'rpm'})
vs.rotation_speed_set = nx.NXfield(3602.3, attrs={'units': 'rpm'})
vs.twist = nx.NXfield(2, attrs={'units': 'deg'})
vs.twist_set = nx.NXfield(2.0, attrs={'units': 'deg'})
vs.twist_axis = nx.NXfield("horizontal")
instrument.velocity_selector = vs

# attenuator
# entry0/instrument/attenuator
# value is the position number on the wheel
attenuator = nx.NXattenuator(type=nx.NXfield("PMMA"),
                             position=nx.NXfield(1235, attrs={'units': 'steps'}),
                             position_set=nx.NXfield(1233, attrs={'units': 'steps'}),
                             value=1
                             )
instrument.attenuator = attenuator

# polarizer
# entry0/instrument/polarizer
instrument.polarizer = nx.NXpolarizer(position=nx.NXfield(23356, attrs={'units': 'deg'}),
                                      position_set=nx.NXfield(2353, attrs={'units': 'deg'}),
                                      state=nx.NXfield('off')
                                      )
# flipper
# entry0/instrument/flipper
instrument.flipper = nx.NXflipper(current=nx.NXfield(2.5, attrs={'units': 'A'}),
                                  current_set=nx.NXfield(2.5, attrs={'units': 'A'}),
                                  voltage=nx.NXfield(25, attrs={'units': 'V'}),
                                  voltage_set=nx.NXfield(25, attrs={'units': 'V'}),
                                  state=nx.NXfield('off')
                                  )

# entry0/instrument/beam_stop
instrument.beam_stop = nx.NXbeam_stop(x=nx.NXfield(2.32, attrs={'units': 'mm'}),
                                      x_set=nx.NXfield(2.32, attrs={'units': 'mm'}),
                                      y=nx.NXfield(2.32, attrs={'units': 'mm'}),
                                      y_set=nx.NXfield(2.32, attrs={'units': 'mm'}),
                                      size=nx.NXfield(25, attrs={'units': 'mm'}),
                                      distance_to_detector=nx.NXfield(35, attrs={'units': 'mm'})
                                      )
# entry0/instrument/collimator
collimator = nx.NXcollimator()

# slits positions read from encoder (4 blades per slits)
# x_gap and y_gap are computed thanks to the blades positions
# the same for the beam center
# previous calibration of the center have to be performed during hot commisionning

collimator.slit_0 = nx.NXslit(left_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              right_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              bottom_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              upper_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              x_gap=nx.NXfield(6.8, attrs={'units': 'mm'}),
                              y_gap=nx.NXfield(6.8, attrs={'units': 'mm'}),
                              x_center=nx.NXfield(0, attrs={'units': 'mm'}),
                              y_center=nx.NXfield(0, attrs={'units': 'mm'})
                              )

collimator.slit_1 = nx.NXslit(left_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              right_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              bottom_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              upper_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              x_gap=nx.NXfield(6.8, attrs={'units': 'mm'}),
                              y_gap=nx.NXfield(6.8, attrs={'units': 'mm'}),
                              x_center=nx.NXfield(0, attrs={'units': 'mm'}),
                              y_center=nx.NXfield(0, attrs={'units': 'mm'})
                              )
collimator.slit_2 = nx.NXslit(left_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              right_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              bottom_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              upper_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              x_gap=nx.NXfield(6.8, attrs={'units': 'mm'}),
                              y_gap=nx.NXfield(6.8, attrs={'units': 'mm'}),
                              x_center=nx.NXfield(0, attrs={'units': 'mm'}),
                              y_center=nx.NXfield(0, attrs={'units': 'mm'})
                              )
collimator.slit_3 = nx.NXslit(left_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              right_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              bottom_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              upper_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              x_gap=nx.NXfield(6.8, attrs={'units': 'mm'}),
                              y_gap=nx.NXfield(6.8, attrs={'units': 'mm'}),
                              x_center=nx.NXfield(0, attrs={'units': 'mm'}),
                              y_center=nx.NXfield(0, attrs={'units': 'mm'})
                              )
collimator.slit_4 = nx.NXslit(left_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              right_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              bottom_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              upper_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              x_gap=nx.NXfield(6.8, attrs={'units': 'mm'}),
                              y_gap=nx.NXfield(6.8, attrs={'units': 'mm'}),
                              x_center=nx.NXfield(0, attrs={'units': 'mm'}),
                              y_center=nx.NXfield(0, attrs={'units': 'mm'})
                              )
collimator.slit_5 = nx.NXslit(left_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              right_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              bottom_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              upper_position=nx.NXfield(3.4, attrs={'units': 'mm'}),
                              x_gap=nx.NXfield(6.8, attrs={'units': 'mm'}),
                              y_gap=nx.NXfield(6.8, attrs={'units': 'mm'}),
                              x_center=nx.NXfield(0, attrs={'units': 'mm'}),
                              y_center=nx.NXfield(0, attrs={'units': 'mm'})
                              )
# guides position read from encoder. According to this value, the guide state is either
# 'in', 'out' ,' 'undefined'. The last case appears for miss positioning of the guide and allow to
# check if there is a problem
collimator.guide_0 = nx.NXguide(position=nx.NXfield(2356, units='steps'),
                                position_set=nx.NXfield(2356, attrs={'units': 'steps'}),
                                m_value=nx.NXfield(3),
                                state='in'
                                )
collimator.guide_1 = nx.NXguide(position=nx.NXfield(2356, attrs={'units': 'steps'}),
                                position_set=nx.NXfield(2356, attrs={'units': 'steps'}),
                                m_value=nx.NXfield(3),
                                state='in'
                                )
collimator.guide_2 = nx.NXguide(position=nx.NXfield(2356, attrs={'units': 'steps'}),
                                position_set=nx.NXfield(2356, attrs={'units': 'steps'}),
                                m_value=nx.NXfield(3),
                                state='in'
                                )
collimator.guide_3 = nx.NXguide(position=nx.NXfield(2356, attrs={'units': 'steps'}),
                                position_set=nx.NXfield(2356, attrs={'units': 'steps'}),
                                m_value=nx.NXfield(3),
                                state='in'
                                )
collimator.guide_4 = nx.NXguide(position=nx.NXfield(2356, attrs={'units': 'steps'}),
                                position_set=nx.NXfield(2356, attrs={'units': 'steps'}),
                                m_value=nx.NXfield(3),
                                state='in'
                                )

# Sourcde aperture computed from the slits and guides positions
# This should be the aperture of the first slit after the last guide that is in 'in' position
collimator.source_aperture_x = nx.NXfield(6.6, attrs={'units': 'mm'})
collimator.source_aperture_y = nx.NXfield(6.6, attrs={'units': 'mm'})
# Diameter of the snout diaphragm to be set be the instrument scientist.
collimator.snout_diameter = nx.NXfield(7.6, attrs={'units': 'mm'})
# pressure inside the collimator (if available)
collimator.pressure = nx.NXfield(0.001, attrs={'units': 'mbar'})
# "Base" collimation length (used to calculate the actual collimation length).
# Retrieved by control software from guide settings.
collimator.length = nx.NXfield(8000.0, attrs={'units': 'mm'})
# collimator/length_offset : Distance between last slit and sample (used to calculate the actual collimation length).
collimator.length_offset = nx.NXfield(50.0, attrs={'units': 'mm'})
collimator.collimation_length = collimator.length_offset + collimator.length

instrument.collimator = collimator
# detector0
instrument.detector0 = nx.NXdetector(count_time=nx.NXfield(2500, attrs={'units': 's'}),
                                     detector_vessel_pressure=nx.NXfield(0.001, attrs={'units': 'mbar'}),
                                     raw_data=nx.NXfield(100*np.random.rand(128, 128)),
                                     data=nx.NXfield(100*np.random.rand(128, 128)),
                                     distance=nx.NXfield(2554, attrs={'units': 'mm'}),
                                     x_pixel_size=nx.NXfield(5.0, attrs={'units': 'mm'}),
                                     y_pixel_size=nx.NXfield(5.0, attrs={'units': 'mm'}),
                                     type=nx.NXfield('monoblock'),
                                     deadtime=nx.NXfield(3.5e-6, attrs={'units': 's'}),
                                     z_position=nx.NXfield(2000.0, attrs={'units': 'mm'}),
                                     z_position_set=nx.NXfield(2000.0, attrs={'units': 'mm'}),
                                     z_position_offset=nx.NXfield(2000.0, attrs={'units': 'mm'}),
                                     x_position=nx.NXfield(10.0, attrs={'units': 'mm'}),
                                     x_position_set=nx.NXfield(10.0, attrs={'units': 'mm'}),
                                     x_position_offset=nx.NXfield(10.0, attrs={'units': 'mm'}),
                                     beam_center_x=nx.NXfield(64.2, attrs={'units': 'pixels'}),
                                     beam_center_y=nx.NXfield(64.2, attrs={'units': 'pixels'}),
                                     pixel_mask=nx.NXfield(np.zeros((128, 128))),
                                     pixel_mask_applied=False
                                     )
# detector 1 (left) :
# raw data 256x16 pixels
# tube size : 12.7mm
# raw_pixel size = 2.5 x 12.7 mm
# targeted filtered pixel size = 10 x 12.7 mm
# filtered data shape = 64x16 pixels
instrument.detector1 = nx.NXdetector(count_time=nx.NXfield(2500, attrs={'units': 's'}),
                                     detector_vessel_pressure=nx.NXfield(0.001, attrs={'units': 'mbar'}),
                                     raw_data=nx.NXfield(np.zeros((256, 16))),
                                     data=nx.NXfield(100*np.random.rand(64, 16)),
                                     distance=nx.NXfield(2554, attrs={'units': 'mm'}),
                                     x_pixel_size=nx.NXfield(12.7, attrs={'units': 'mm'}),
                                     y_pixel_size=nx.NXfield(10.0, attrs={'units': 'mm'}),
                                     type=nx.NXfield('monoblock'),
                                     deadtime=nx.NXfield(3.5e-6, attrs={'units': 's'}),
                                     z_position=nx.NXfield(2000.0, attrs={'units': 'mm'}),
                                     z_position_set=nx.NXfield(2000.0, attrs={'units': 'mm'}),
                                     z_position_offset=nx.NXfield(2000.0, attrs={'units': 'mm'}),
                                     x_position=nx.NXfield(10.0, attrs={'units': 'mm'}),
                                     x_position_set=nx.NXfield(10.0, attrs={'units': 'mm'}),
                                     x_position_offset=nx.NXfield(10.0, attrs={'units': 'mm'}),
                                     beam_center_x=nx.NXfield(64.2, attrs={'units': 'pixels'}),
                                     beam_center_y=nx.NXfield(64.2, attrs={'units': 'pixels'}),
                                     pixel_mask=nx.NXfield(np.zeros((64, 16))),
                                     pixel_mask_applied=False,
                                     )
# detector 2 (bottom)
instrument.detector2 = nx.NXdetector(count_time=nx.NXfield(2500, attrs={'units': 's'}),
                                     detector_vessel_pressure=nx.NXfield(0.001, attrs={'units': 'mbar'}),
                                     raw_data=nx.NXfield(100*np.random.rand(16, 256)),
                                     data=nx.NXfield(100*np.random.rand(16, 64)),
                                     distance=nx.NXfield(2554, attrs={'units': 'mm'}),
                                     x_pixel_size=nx.NXfield(10, attrs={'units': 'mm'}),
                                     y_pixel_size=nx.NXfield(12.7, attrs={'units': 'mm'}),
                                     type=nx.NXfield('monoblock'),
                                     deadtime=nx.NXfield(3.5e-6, attrs={'units': 's'}),
                                     z_position=nx.NXfield(2000.0, attrs={'units': 'mm'}),
                                     z_position_set=nx.NXfield(2000.0, attrs={'units': 'mm'}),
                                     z_position_offset=nx.NXfield(2000.0, attrs={'units': 'mm'}),
                                     x_position=nx.NXfield(10.0, attrs={'units': 'mm'}),
                                     x_position_set=nx.NXfield(10.0, attrs={'units': 'mm'}),
                                     x_position_offset=nx.NXfield(10.0, attrs={'units': 'mm'}),
                                     beam_center_x=nx.NXfield(64.2, attrs={'units': 'pixels'}),
                                     beam_center_y=nx.NXfield(64.2, attrs={'units': 'pixels'}),
                                     pixel_mask=nx.NXfield(np.zeros((16, 64))),
                                     pixel_mask_applied=False,
                                     )

# Distance between sample and detector's vessel window (used to calculate the actual sample to detector distance).
# Set by instrument scientist at the beginning of the experiment.
instrument.sample2window_distance = nx.NXfield(250.0, attrs={'units': 'mm'})

# Definition of the monitors here we take into account 4 monitors according to the description given
# by Joachim
# entry0/monitor0
entry0.monitor0 = nx.NXmonitor(mode='timer',
                               preset=nx.NXfield(2500, attrs={'units': 's'}),
                               integral=nx.NXfield(25000),
                               description='monitor of proton beam'
                               )
# entry0/monitor1
entry0.monitor1 = nx.NXmonitor(mode='timer',
                               preset=nx.NXfield(2500, attrs={'units': 's'}),
                               integral=nx.NXfield(25000),
                               description='monitor before the beam shutter'
                               )
# entry0/monitor2
entry0.monitor2 = nx.NXmonitor(mode='timer',
                               preset=nx.NXfield(2500, attrs={'units': 's'}),
                               integral=nx.NXfield(25000),
                               description=' monitor after the shutters and before the selector '
                               )
# entry0/monitor3
entry0.monitor3 = nx.NXmonitor(mode='timer',
                               preset=nx.NXfield(2500, attrs={'units': 's'}),
                               integral=nx.NXfield(25000),
                               description='monitor at the exit of the neutron guide'
                               )

# user definition
entry0.user = nx.NXuser(user_name='Bob',
                        local_contact='Hary',
                        proposal_id='2020-9-235',
                        proposal_title='Nice experiment',
                        user_email='user@gmail.com'
                        )
# sample definition
# the sample2window distance is
sample = nx.NXsample(sample_name='Good sample',
                     desciption='it is indeed a good sample',
                     transmission=1.0,
                     thickness=nx.NXfield(1.0, attrs={'units': 'mm'})
                     )

# sample environement
sample_environment = nx.NXenvironment(name='sample changer',
                                      type='type of sample'
                                      )
# sample table
sample_environment.table = nx.NXcollection(z=nx.NXfield(35.3, attrs={'units': 'mm'}),
                                           z_set=nx.NXfield(35.3, attrs={'units': 'mm'}),
                                           x=nx.NXfield(35.3, attrs={'units': 'mm'}),
                                           x_set=nx.NXfield(35.3, attrs={'units': 'mm'}),
                                           y=nx.NXfield(35.3, attrs={'units': 'mm'}),
                                           y_set=nx.NXfield(35.3, attrs={'units': 'mm'}),
                                           rotation=nx.NXfield(0.0, attrs={'units': 'deg'}),
                                           goniometer=nx.NXfield(1.0, attrs={'units': 'deg'})
                                           )
# here we define sub environement entries taht will be store in the NXenvironement field according to the presnce or the abscence of the sample environnement

# sample.sample_changer = nx.NXcollection(x_position=nx.NXfield(1.0, attrs={'units': 'mm'}),
#                                         y_position=nx.NXfield(1.0, attrs={'units': 'mm'}),
#                                         temperature=nx.NXfield(1.0, attrs={'units': '°C'}),
#                                         temperature_set=nx.NXfield(1.0, attrs={'units': '°C'})
#                                         )
# sample.magnet = nx.NXcollection(type='cryomagnet',
#                                 current=nx.NXfield(2.5, attrs={'units': 'A'}),
#                                 current_set=nx.NXfield(2.5, attrs={'units': 'A'}),
#                                 voltage=nx.NXfield(2.5, attrs={'units': 'V'}),
#                                 voltage_set=nx.NXfield(2.5, attrs={'units': 'A'}),
#                                 magnetic_field=nx.NXfield(0.5, attrs={'units': 'T'}),
#                                 magnetic_set=nx.NXfield(0.5, attrs={'units': 'T'})
#                                 )

entry0.insert(sample, name='sample')
entry0.insert(sample_environment, name='environement')

# data : link to detector data
entry0.data0 = nx.NXdata(attrs={'interpretation': b"image",
                                'signal': "data"})
entry0.data0.makelink(entry0.instrument.detector0.data)

entry0.data1 = nx.NXdata(attrs={'interpretation': b"image",
                                'signal': "data"})
entry0.data1.makelink(entry0.instrument.detector1.data)

entry0.data2 = nx.NXdata(attrs={'interpretation': b"image",
                                'signal': "data"})
entry0.data2.makelink(entry0.instrument.detector2.data)

# save the NXroot to the file
try:
    filepath = '/home/achennev/python/pygdatax/src/dataformat/test_nexus_AC_v2.nxs'
    root.save(filename=filepath, mode='w')
    # root.save(filename='test_nexus_AC_v2.nxs', mode='w')
    import os
    os.system('silx view '+filepath)
except nx.NeXusError:
    print('pb')
