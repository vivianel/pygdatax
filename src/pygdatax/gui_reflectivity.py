import sys
import os
import warnings
from matplotlib import mplDeprecation
import numpy as np
import fabio
import silx.gui.hdf5
from silx.gui import qt, colors
from silx.gui.plot import PlotWindow, Profile
from silx .gui.plot.tools.roi import RegionOfInterestManager
from silx.gui.plot.items.roi import RectangleROI, CrossROI
import silx.io as sio
from silx.io.utils import is_group, is_dataset, is_file
from silx.io.nxdata import is_NXroot_with_default_NXdata, get_default
import nexusformat.nexus as nx
from numpy.random import randint
from scipy.ndimage.measurements import center_of_mass
from pygdatax.icons import getQIcon
from pygdatax import nxlib
from pathlib import Path
import yaml


class DataView(PlotWindow):

    def __init__(self):
        super().__init__(backend=None, resetzoom=True,
                         autoScale=True, logScale=True,
                         grid=False, curveStyle=True, colormap=True,
                         aspectRatio=True, yInverted=True,
                         copy=True, save=True, print_=True,
                         control=True, position= [('X', lambda x, y: x),
                                                  ('Y', lambda x, y: y),
                                                  ('Data', self._zValue)],
                         roi=False, mask=True, fit=True)
        # """Widget displaying information"""
        # posInfo = [('X', lambda x, y: x),
        #            ('Y', lambda x, y: y),
        #            ('Data', self._zValue)]
        # plot widget
        """Widget displaying information"""
        posInfo = [('X', lambda x, y: x),
                   ('Y', lambda x, y: y),
                   ('Data', self._zValue)]
        self.plotWindow = PlotWindow(backend=None, resetzoom=True,
                                     autoScale=True, logScale=True,
                                     grid=False, curveStyle=True, colormap=True,
                                     aspectRatio=True, yInverted=True,
                                     copy=True, save=True, print_=True,
                                     control=True, position=posInfo,
                                     roi=False, mask=True, fit=True)
        self.roiManager = RegionOfInterestManager(self)
        self.roiManager.setColor('pink')
        self.roiManager.sigRoiAdded.connect(self.updateAddedRegionOfInterest)
        self.addToolBarBreak()
        action = self.roiManager.getInteractionModeAction(RectangleROI)
        toolbar = qt.QToolBar('')
        toolbar.addAction(action)
        findCenterAction = qt.QAction(self)
        findCenterAction.setIcon(getQIcon('target.ico'))
        findCenterAction.setToolTip('find beam on the current ROI')
        findCenterAction.triggered.connect(self.findCenter)
        toolbar.addAction(findCenterAction)
        self.addToolBar(toolbar)
        legendWidget = self.getLegendsDockWidget()
        # self.plotWindow.setDockOptions()
        legendWidget.setAllowedAreas(qt.Qt.TopDockWidgetArea)
        self.profileTools = Profile.ProfileToolBar(parent=self,
                                                   plot=self)
        self.addToolBar(self.profileTools)
        self.setDefaultColormap(colors.Colormap(name='jet', normalization='log',
                                                vmin=None, vmax=None, autoscaleMode='stddev3')
                                )

    def updateAddedRegionOfInterest(self, roi):
        """Called for each added region of interest: set the name"""
        roisList = self.roiManager.getRois()
        if len(roisList)>1:
            self.roiManager.removeRoi(roisList[0])
        roisList = self.roiManager.getRois()
        if roi.getName() == '':
            roi.setName('ROI %d' % len(self.roiManager.getRois()))
        # if isinstance(roi, LineMixIn):
        #     roi.setLineWidth(1)
        #     roi.setLineStyle('--')
        # if isinstance(roi, SymbolMixIn):
        #     roi.setSymbolSize(5)
        if isinstance(roi, RectangleROI):
            roi.setSelectable(True)
            roi.setEditable(True)

    def _zValue(self, x, y):
        value = '-'
        valueZ = - float('inf')
        for image in self.getAllImages():
            data = image.getData(copy=False)
            if image.getZValue() >= valueZ:  # This image is over the previous one
                ox, oy = image.getOrigin()
                sx, sy = image.getScale()
                row, col = (y - oy) / sy, (x - ox) / sx
                if row >= 0 and col >= 0:
                    # Test positive before cast otherwise issue with int(-0.5) = 0
                    row, col = int(row), int(col)
                    if row < data.shape[0] and col < data.shape[1]:
                        value = data[row, col]
                        valueZ = image.getZValue()
        return value

    def findCenter(self):
        roisList = self.roiManager.getRois()
        if roisList:
            if isinstance(roisList[0], RectangleROI):
                point = roisList[0].getOrigin()
                size = roisList[0].getSize()
                for image in self.getAllImages():
                    data = image.getData(copy=False)
                    center = center_of_mass(data[int(point[1]):int(point[1])+int(size[1]),int(point[0]):int(point[0])+int(size[0])])
                    roiCenter = CrossROI()
                    absoluteCenter = [center[1]+int(point[0]), center[0]+int(point[1])]
                    roiCenter.setPosition(absoluteCenter)
                    label = 'beam center\n x0 : %.3f y0: %.3f' % (absoluteCenter[0], absoluteCenter[1])
                    roiCenter.setName(label)
                    self.roiManager.addRoi(roiCenter)

            #
            # if image.getZValue() >= valueZ:  # This image is over the previous one
            #     ox, oy = image.getOrigin()
            #     sx, sy = image.getScale()
            #     row, col = (y - oy) / sy, (x - ox) / sx
            #     if row >= 0 and col >= 0:
            #         # Test positive before cast otherwise issue with int(-0.5) = 0
            #         row, col = int(row), int(col)
            #         if row < data.shape[0] and col < data.shape[1]:
            #             value = data[row, col]
            #             valueZ = image.getZValue()