import sys
import os
import warnings
from matplotlib import mplDeprecation
import numpy as np
import fabio
import silx.gui.hdf5
from silx.gui import qt, colors
from silx.gui.plot import PlotWindow, Profile, PlotWidget
from silx .gui.plot.tools.roi import RegionOfInterestManager
from silx.gui.plot.items.roi import RectangleROI, CrossROI
import silx.io as sio
from silx.io.utils import is_group, is_dataset, is_file
from silx.io.nxdata import is_NXroot_with_default_NXdata, get_default
import nexusformat.nexus as nx
from numpy.random import randint
from scipy.ndimage.measurements import center_of_mass
from pathlib import Path
import yaml
from pygdatax.icons import getQIcon
from pygdatax import nxlib, moduledescription
from pygdatax.instruments import xeuss


def get_edf_description(edf_filepath):
    des = []
    if os.path.isfile(edf_filepath):
        if edf_filepath.endswith('.edf'):
            dataObj = fabio.open(edf_filepath)
            try:
                des.append(os.path.basename(edf_filepath))
                des.append(dataObj.header['Comment'])
                des.append(str(1.5))
                des.append(dataObj.header['pilroi1'])
                distance = float(dataObj.header['SampleDistance'])
                # convert to mm
                distance *= 1000.0
                des.append(str(distance))
                des.append(dataObj.header['count_time'])
                des.append(dataObj.header['Date'])
            except KeyError:
                des.append(os.path.split(edf_filepath)[1])
                des += 5 * ['']
    return des


def get_nxs_description(nxs_filepath):
    if os.path.isfile(nxs_filepath):
        if nxs_filepath.endswith('.nxs'):
            try:
                des = []
                root = nx.nxload(nxs_filepath, mode='r')
                entry = root[nxlib.get_last_entry_key(root)]
                des.append(os.path.basename(nxs_filepath))
                des.append(str(entry.sample.sample_name.nxdata.decode()))
                des.append(str(str(entry.sample.thickness.nxdata)))
                des.append(str(entry.sample.transmission.nxdata))
                des.append(str(entry.instrument.detector.distance.nxdata))
                des.append(str(entry.sample.count_time.nxdata))
                des.append(root.attrs['file_time'])
                root.close()
            except (KeyError, nx.NeXusError, ValueError):
                des = []
                des.append(os.path.split(nxs_filepath)[1])
                des += 5 * ['']
            # compatibility with windows
            except AttributeError:
                des = []
                root = nx.nxload(nxs_filepath, mode='r')
                entry = root[nxlib.get_last_entry_key(root)]
                des.append(os.path.basename(nxs_filepath))
                des.append(str(entry.sample.sample_name.nxdata))
                des.append(str(str(entry.sample.thickness.nxdata)))
                des.append(str(entry.sample.transmission.nxdata))
                des.append(str(entry.instrument.detector.distance.nxdata))
                des.append(str(entry.sample.count_time.nxdata))
                des.append(root.attrs['file_time'])
                root.close()
    return des


class EdfFileTable(qt.QTableWidget):
    directory = ''
    file_extension = '.edf'
    fileSelectedChanged = qt.pyqtSignal(str)
    emptyCellFile = None
    darkFile = None
    emptyBeamFile = None
    maskFile = None
    trashFiles = []
    treatedFiles = []

    def __init__(self):
        super(EdfFileTable, self).__init__()
        self.setColumnCount(6)
        self.setRowCount(3)
        self.setRowCount(4)
        self.setSelectionBehavior(qt.QAbstractItemView.SelectRows)
        self.setHorizontalHeaderLabels(['File', 'comment', 'e (mm)','tr', 'distance',  'counting time', ' date'])
        self.currentItemChanged.connect(self.on_selectionChanged)
        self.setContextMenuPolicy(qt.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.generateMenu)

    def setDirectory(self, directory):
        folderPath = Path(directory)
        if folderPath.is_dir():
            self.directory = folderPath
            self.refresh()

    def refresh(self):
        self.currentItemChanged.disconnect()
        if os.path.isdir(self.directory):
            l = os.listdir(self.directory)
            # l.sort()
            fileList = []
            for item in l:
                if os.path.splitext(item)[1] == self.file_extension:
                    fileList.append(item)
            # self.clearContents()
            self.setRowCount(len(fileList))
            for i, file in enumerate(fileList):
                description = get_edf_description(os.path.join(self.directory, file))
                for j, des in enumerate(description):
                    item = qt.QTableWidgetItem(des)
                    if j == 2:
                        item.setFlags(qt.Qt.ItemIsEditable | qt.Qt.ItemIsEnabled | qt.Qt.ItemIsSelectable)
                    else:
                        item.setFlags(qt.Qt.ItemIsSelectable | qt.Qt.ItemIsEnabled)
                    self.setItem(i, j, item)
                # check of the file was current set to parametric file
                filepath = os.path.join(self.directory, file)
                if filepath == self.darkFile:
                    self.set_row_bkg(i, qt.QColor("blue"))
                    self.item(i, 0).setIcon(getQIcon('dark.ico'))
                elif filepath == self.emptyCellFile:
                    self.set_row_bkg(i, qt.QColor("cyan"))
                    self.item(i, 0).setIcon(getQIcon('empty_cell.ico'))
                elif filepath == self.emptyBeamFile:
                    self.set_row_bkg(i, qt.QColor("red"))
                    self.item(i, 0).setIcon(getQIcon('beam.ico'))
                # elif filepath in self.trashFiles:
                #     self.set_row_bkg(i, qt.QColor("grey"))
                #     self.item(i, 0).setIcon(getQIcon('cross.ico'))
                elif filepath == self.maskFile:
                    self.set_row_bkg(i, qt.QColor("white"))
                    self.item(i, 0).setIcon(getQIcon('mask.ico'))
                elif os.path.exists(os.path.splitext(filepath)[0]+'.nxs'):
                    self.item(i, 0).setIcon(getQIcon('check.ico'))

        self.sortItems(0, qt.Qt.AscendingOrder)
        self.currentItemChanged.connect(self.on_selectionChanged)

    def on_selectionChanged(self):
        row = self.currentRow()
        file = self.item(row, 0).text()
        self.fileSelectedChanged.emit(os.path.join(self.directory, file))

    def generateMenu(self, event):
        current_item = self.itemAt(event)
        # current_item = self.selectedItems()
        menu = qt.QMenu()
        # emptyCellAction = qt.QAction(qt.QIcon(qt.QPixmap('../ressources/empty_cell.ico')), 'empty cell')
        emptyCellAction = qt.QAction(getQIcon('empty_cell.ico'), 'empty cell')
        emptyCellAction.triggered.connect(self._set_empty_cell)
        darkAction = qt.QAction(getQIcon('dark.ico'), 'dark')
        darkAction.triggered.connect(self._set_dark)
        emptyBeamAction = qt.QAction(getQIcon('beam.ico'), 'empty beam')
        emptyBeamAction.triggered.connect(self._set_empty_beam)
        # trashAction = qt.QAction(getQIcon('cross.ico'), 'trash')
        # trashAction.triggered.connect(self._set_trash)
        sampleAction = qt.QAction('sample')
        sampleAction.triggered.connect(self._set_sample)
        maskAction = qt.QAction(getQIcon('mask.ico'), 'mask')
        maskAction.triggered.connect(self._set_mask)
        # build menu
        menu.addAction(darkAction)
        menu.addAction(emptyCellAction)
        menu.addAction(emptyBeamAction)
        # menu.addAction(trashAction)
        menu.addAction(sampleAction)
        menu.addAction(maskAction)
        action = menu.exec_(self.mapToGlobal(event))
        # print('###################################### \n')
        # print(('dark : %s \n '
        #        'empty cell: %s \n'
        #        'empty beam %s') %
        #       (self.darkFile, self.emptyCellFile, self.emptyBeamFile))

    def set_row_bkg(self, row, color):
        for i in range(self.columnCount()):
            item = self.item(row, i)
            item.setBackground(color)

    def _set_empty_cell(self):
        current_item = self.currentItem()
        if current_item is not None:
            current_ec_item = self.findItems(os.path.basename(str(self.emptyCellFile)), qt.Qt.MatchExactly)
            # remove the previous empty cell icons
            if current_ec_item:
                self.set_row_bkg(current_ec_item[0].row(), qt.QColor("white"))
                filepath = os.path.join(self.directory, current_ec_item[0].text())
                if os.path.exists(os.path.splitext(filepath)[0] + '.nxs'):
                    current_ec_item[0].setIcon(getQIcon('check.ico'))
                else:
                    current_ec_item[0].setIcon(qt.QIcon())
            row = current_item.row()
            ncol = self.columnCount()
            first_col_item = self.item(row, 0)
            file = first_col_item.text()
            first_col_item.setIcon(getQIcon('empty_cell.ico'))
            self.set_row_bkg(row, qt.QColor("cyan"))
            fullfile = os.path.join(self.directory, file)
            self.emptyCellFile = fullfile
            # remove double reference
            if fullfile == self.darkFile:
                self.darkFile = None
            elif fullfile == self.emptyBeamFile:
                self.emptyBeamFile = None
            elif fullfile == self.maskFile:
                self.maskFile = None
            else:
                pass

    def _set_sample(self):
        for current_item in self.selectedItems():
            if current_item is not None:
                row = current_item.row()
                ncol = self.columnCount()
                first_col_item = self.item(row, 0)
                file = first_col_item.text()
                self.set_row_bkg(row, qt.QColor("white"))
                first_col_item.setIcon(qt.QIcon())
                fullfile = os.path.join(self.directory, file)
                if os.path.exists(os.path.splitext(fullfile)[0] + '.nxs'):
                    first_col_item.setIcon(getQIcon('check.ico'))
                else:
                    first_col_item.setIcon(qt.QIcon())
                # remove double reference
                if fullfile == self.emptyCellFile:
                    self.emptyCellFile = None
                elif fullfile == self.darkFile:
                    self.darkFile = None
                elif fullfile == self.emptyBeamFile:
                    self.emptyBeamFile = None
                elif fullfile == self.maskFile:
                    self.maskFile = None
                # elif fullfile in self.trashFiles:
                #     self.trashFiles.remove(fullfile)

    def _set_dark(self, event):
        current_item = self.currentItem()
        if current_item is not None:
            current_dark_item = self.findItems(os.path.basename(str(self.darkFile)), qt.Qt.MatchExactly)
            # remove the previous empty cell icons
            if current_dark_item:
                self.set_row_bkg(current_dark_item[0].row(), qt.QColor("white"))
                filepath = os.path.join(self.directory, current_dark_item[0].text())
                if os.path.exists(os.path.splitext(filepath)[0] + '.nxs'):
                    current_dark_item[0].setIcon(getQIcon('check.ico'))
                else:
                    current_dark_item[0].setIcon(qt.QIcon())
            row = current_item.row()
            ncol = self.columnCount()
            first_col_item = self.item(row, 0)
            file = first_col_item.text()
            self.set_row_bkg(row, qt.QColor("blue"))
            first_col_item.setIcon(getQIcon('dark.ico'))
            # remove double reference
            fullfile = os.path.join(self.directory, file)
            self.darkFile = fullfile
            if fullfile == self.emptyCellFile:
                self.emptyCellFile = None
            elif fullfile == self.emptyBeamFile:
                self.emptyBeamFile = None
            elif fullfile == self.maskFile:
                self.maskFile = None
            else:
                pass

    def _set_empty_beam(self):
        current_item = self.currentItem()
        if current_item is not None:
            current_eb_item = self.findItems(os.path.basename(str(self.emptyBeamFile)), qt.Qt.MatchExactly)
            # remove the previous empty cell icons
            if current_eb_item:
                self.set_row_bkg(current_eb_item[0].row(), qt.QColor("white"))
                filepath = os.path.join(self.directory, current_eb_item[0].text())
                if os.path.exists(os.path.splitext(filepath)[0] + '.nxs'):
                    current_eb_item[0].setIcon(getQIcon('check.ico'))
                else:
                    current_eb_item[0].setIcon(qt.QIcon())
            row = current_item.row()
            ncol = self.columnCount()
            first_col_item = self.item(row, 0)
            file = first_col_item.text()
            self.set_row_bkg(row, qt.QColor("red"))
            first_col_item.setIcon(getQIcon('beam.ico'))
            # remove double reference
            fullfile = os.path.join(self.directory, file)
            self.emptyBeamFile = fullfile
            if fullfile == self.emptyCellFile:
                self.emptyCellFile = None
            elif fullfile == self.darkFile:
                self.darkFile = None
            elif fullfile == self.maskFile:
                self.maskFile = None
            else:
                pass

    # def _set_trash(self):
    #     # can be applied to the overall selection
    #     for current_item in self.selectedItems():
    #         if current_item is not None:
    #             row = current_item.row()
    #             ncol = self.columnCount()
    #             first_col_item = self.item(row, 0)
    #             file = first_col_item.text()
    #             self.set_row_bkg(row, qt.QColor("grey"))
    #             first_col_item.setIcon(getQIcon('cross.ico'))
    #             # remove double reference
    #             fullfile = os.path.join(self.directory, file)
    #             self.trashFiles.append(fullfile)
    #             if fullfile == self.emptyCellFile:
    #                 self.emptyCellFile = None
    #             elif fullfile == self.darkFile:
    #                 self.darkFile = None
    #             elif fullfile == self.emptyBeamFile:
    #                 self.emptyBeamFile = None
    #             elif fullfile == self.maskFile:
    #                 self.maskFile = None
    #             else:
    #                 pass

    def _set_mask(self):
        current_item = self.currentItem()
        if current_item is not None:
            current_mask_item = self.findItems(os.path.basename(str(self.maskFile)), qt.Qt.MatchExactly)
            # remove the previous empty cell icons
            if current_mask_item:
                self.set_row_bkg(current_mask_item[0].row(), qt.QColor("white"))
                filepath = os.path.join(self.directory, current_mask_item[0].text())
                if os.path.exists(os.path.splitext(filepath)[0] + '.nxs'):
                    current_mask_item[0].setIcon(getQIcon('check.ico'))
                else:
                    current_mask_item[0].setIcon(qt.QIcon())
            row = current_item.row()
            ncol = self.columnCount()
            first_col_item = self.item(row, 0)
            file = first_col_item.text()
            first_col_item.setIcon(getQIcon('mask.ico'))
            # self.set_row_bkg(row, qt.QColor("cyan"))

            fullfile = os.path.join(self.directory, file)
            self.maskFile = fullfile
            # remove double reference
            if fullfile == self.darkFile:
                self.darkFile = None
            elif fullfile == self.emptyBeamFile:
                self.emptyBeamFile = None
            elif fullfile == self.emptyCellFile:
                self.maskFile = None
            else:
                pass

    def get_sample_files(self):
        sampleList = []
        thicknessList = []
        for current_item in self.selectedItems():
            row = current_item.row()
            ncol = self.columnCount()
            first_col_item = self.item(row, 0)
            file = first_col_item.text()
            # remove double reference
            fullfile = os.path.join(self.directory, file)
            self.trashFiles.append(fullfile)
            if fullfile not in [self.emptyCellFile, self.darkFile, self.emptyBeamFile]:
                sampleList.append(fullfile)
                try:
                    thicknessList.append(float(self.item(row, 2).text()))
                except ValueError:
                    thicknessList.append(1.5)
        return sampleList, thicknessList


class EdfTreatmentWidget(qt.QWidget):
    edfSelectionChanged = qt.pyqtSignal(str)
    edfTreatmentClicked = qt.pyqtSignal()

    def __init__(self, parent=None):
        super(EdfTreatmentWidget, self).__init__()
        # self.directoryLineEdit = qt.QLineEdit(parent=self)
        # self.directoryPickerButton = qt.QPushButton()
        # self.directoryPickerButton.setIcon(qt.QIcon(qt.QPixmap('../ressources/directory.ico')))
        # self.refreshButton = qt.QPushButton()
        # self.refreshButton.setIcon(qt.QIcon(qt.QPixmap('../ressources/refresh.ico')))
        self.table = EdfFileTable()
        # beam center coordinates
        self.x0LineEdit = qt.QLineEdit()
        # self.x0LineEdit.setValidator(qt.QDoubleValidator())
        self.y0LineEdit = qt.QLineEdit()
        # self.y0LineEdit.setValidator(qt.QDoubleValidator())
        # sample to detector distance
        self.distanceLineEdit = qt.QLineEdit()
        # self.distanceLineEdit.setValidator(qt.QDoubleValidator())
        # define the number of bins for azimutal averaging
        self.binsLineEdit = qt.QLineEdit('900')
        self.binsLineEdit.setValidator(qt.QIntValidator())
        # load and save integration parameters
        self.saveConfigButton = qt.QPushButton('save')
        self.saveConfigButton.setToolTip('Save treatment parameters\n'
                                         'and the subtraction files')
        self.loadConfigButton = qt.QPushButton('load')
        self.loadConfigButton.setToolTip('Save treatment parameters\n'
                                         'and the subtraction files')

        # button to treat data
        self.treatButton = qt.QPushButton('treat selected')
        self.treatButton.setIcon(getQIcon('gear.ico'))
        self.treatButton.setToolTip('Perform azimutal integration and \n '
                                    'data substraction\n'
                                    'on selected files')
        # parameter form layout
        formLayout = qt.QFormLayout()
        formLayout.addRow('x0 (pixels):', self.x0LineEdit)
        formLayout.addRow('y0 (pixels):', self.y0LineEdit)
        formLayout.addRow('distance (mm):', self.distanceLineEdit)
        formLayout.addRow('bins :', self.binsLineEdit)
        # parameter total layout
        paramLayout = qt.QHBoxLayout()
        configLayout= qt.QVBoxLayout()
        configLayout.addWidget(self.loadConfigButton)
        configLayout.addWidget(self.saveConfigButton)
        paramLayout.addLayout(formLayout)
        paramLayout.addLayout(configLayout)
        # general layout
        hlayout = qt.QHBoxLayout()
        # hlayout.addWidget(qt.QLabel('directory :'))
        # hlayout.addWidget(self.directoryLineEdit)
        # hlayout.addWidget(self.directoryPickerButton)
        # hlayout.addWidget(self.refreshButton)
        vlayout = qt.QVBoxLayout()
        vlayout.addLayout(paramLayout)
        vlayout.addLayout(hlayout)
        vlayout.addWidget(self.table)
        vlayout.addWidget(self.treatButton)
        self.setLayout(vlayout)
        # connect signals
        # self.directoryLineEdit.textChanged.connect(self.set_directory)
        # self.directoryPickerButton.clicked.connect(self.choose_directory)
        # self.refreshButton.clicked.connect(self.table.refresh)
        # we unconnect the treatbutton here because of interactions with treeview
        # self.treatButton.clicked.connect(self.treat)
        self.table.fileSelectedChanged.connect(self.on_file_selected)
        self.saveConfigButton.clicked.connect(self.saveConfig_clicked)
        self.loadConfigButton.clicked.connect(self.loadConfig_clicked)

    def on_file_selected(self, file):
        self.edfSelectionChanged.emit(file)

    def treat(self):
        directory = self.table.directory
        try:
            x0 = float(self.x0LineEdit.text())
        except ValueError:
            x0 = None
        try:
            y0 = float(self.y0LineEdit.text())
        except ValueError:
            y0 = None
        try:
            distance = float(self.distanceLineEdit.text())
        except ValueError:
            distance = None
        nbins = int(self.binsLineEdit.text())
        mask_file = self.table.maskFile
        # treat the reference files
        # dark file
        dark_file = self.table.darkFile
        if dark_file is not None:
            nxlib.build_nexus_from_edf(dark_file)
            dark_file = dark_file.split('.')[0] + '.nxs'
            xeuss.set_beam_center(dark_file, x0=x0, y0=y0)
            xeuss.azimutal_integration(dark_file, bins=nbins, mask=mask_file)
        # empty cell
        ec_file = self.table.emptyCellFile
        if ec_file is not None:
            nxlib.build_nexus_from_edf(ec_file)
            ec_file = ec_file.split('.')[0] + '.nxs'
            xeuss.set_beam_center(ec_file, x0=x0, y0=y0)
            xeuss.azimutal_integration(ec_file, bins=nbins, mask=mask_file)
        # empty beam
        eb_file = self.table.emptyBeamFile
        if eb_file is not None:
            nxlib.build_nexus_from_edf(eb_file)
            eb_file = eb_file.split('.')[0] + '.nxs'
            xeuss.set_beam_center(eb_file, x0=x0, y0=y0)
            xeuss.azimutal_integration(eb_file, bins=nbins, mask=mask_file)

        sampleList, thicknessList = self.table.get_sample_files()
        for sample, thickness in zip(sampleList, thicknessList):
            # nxfile = sample.split('.')[0]+'.nxs'
            try:
                nxlib.build_nexus_from_edf(sample)
                file = sample.split('.')[0]+'.nxs'
                xeuss.set_beam_center(file, x0=x0, y0=y0)  # direct_beam_file=directbeam, new_entry=False)
                xeuss.azimutal_integration(file, bins=nbins, mask=mask_file)
                xeuss.resu(file, dark_file=dark_file, ec_file=ec_file, eb_file=eb_file,
                           distance=distance, thickness=thickness)
            except (KeyError, ValueError):
                print(('%s was ignored during the treatment') % (sample))
        self.table.refresh()

    def saveConfig_clicked(self):
        dic = {}
        if self.table.directory:
            basedir = self.table.directory
        else:
            basedir = os.path.expanduser("~")
        fname, ext = qt.QFileDialog.getSaveFileName(self, 'Save treatment parameters', str(basedir),'YAML files (*.yaml);; all files (*.*)',
                                                      options=qt.QFileDialog.DontUseNativeDialog)
        if fname:
            try:
                dic['x0'] = float(self.x0LineEdit.text())
            except ValueError:
                dic['x0'] = None
            try:
                dic['y0'] = float(self.y0LineEdit.text())
            except ValueError:
                dic['y0'] = None
            try:
                dic['distance'] = float(self.distanceLineEdit.text())
            except ValueError:
                dic['distance'] = None
            dic['nbins'] = int(self.binsLineEdit.text())
            dic['mask_file'] = self.table.maskFile
            dic['dark_file'] = self.table.darkFile
            dic['ec_file'] = self.table.emptyCellFile
            dic['eb_file'] = self.table.emptyBeamFile
            p = Path(fname)
            if p.suffix != '.yaml':
                fname = str(p.with_suffix('.yaml'))
            with open(fname, 'w') as fid:
                yaml.dump(dic, fid)

    def loadConfig_clicked(self):
        if self.table.directory:
            basedir = self.table.directory
        else:
            basedir = os.path.expanduser("~")
        fname, ext = qt.QFileDialog.getOpenFileName(self, 'Load treatment parameters', str(basedir),
                                                    'YAML files (*.yaml);; all files (*.*)',
                                                    options=qt.QFileDialog.DontUseNativeDialog)
        if fname:
            with open(fname, 'r') as fid:
                params = yaml.safe_load(fid)
            if params['x0'] is not None:
                self.x0LineEdit.setText(str(params['x0']))
            else:
                self.x0LineEdit.setText('')
            if params['y0'] is not None:
                self.y0LineEdit.setText(str(params['y0']))
            else:
                self.y0LineEdit.setText('')
            if params['distance'] is not None:
                self.distanceLineEdit.setText(str(params['distance']))
            else:
                self.distanceLineEdit.setText('')
            if params['nbins'] is not None:
                self.binsLineEdit.setText(str(params['nbins']))
            else:
                self.binsLineEdit.setText('900')
            if self.table.directory:
                if os.path.dirname(params['ec_file']) == self.table.directory:
                    self.table.emptyCellFile = params['ec_file']
                if os.path.dirname(params['dark_file']) == self.table.directory:
                    self.table.darkFile = params['dark_file']
                if os.path.dirname(params['eb_file']) == self.table.directory:
                    self.table.emptyBeamFile = params['eb_file']
                if os.path.dirname(params['mask_file']) == self.table.directory:
                    self.table.maskFile = params['mask_file']
                self.table.refresh()


class FileSurvey(qt.QWidget):

    def __init__(self):
        super(FileSurvey, self).__init__()
        self.directoryLineEdit = qt.QLineEdit(parent=self)
        self.directoryPickerButton = qt.QPushButton()
        self.directoryPickerButton.setIcon(getQIcon('directory.ico'))
        self.directoryPickerButton.setToolTip('open data directory')
        self.refreshButton = qt.QPushButton()
        self.refreshButton.setIcon(getQIcon('refresh.ico'))
        self.refreshButton.setToolTip('refresh directory')

        self.tabWidget = qt.QTabWidget()
        self.edfTab = EdfTreatmentWidget()
        self.nxsTab = NexusTreatmentWidget()
        self.tabWidget.addTab(self.edfTab, 'edf data')
        self.tabWidget.addTab(self.nxsTab, 'treated data')

        # layout
        hlayout = qt.QHBoxLayout()
        hlayout.addWidget(qt.QLabel('directory :'))
        hlayout.addWidget(self.directoryLineEdit)
        hlayout.addWidget(self.directoryPickerButton)
        hlayout.addWidget(self.refreshButton)
        vlayout = qt.QVBoxLayout()
        vlayout.addLayout(hlayout)
        vlayout.addWidget(self.tabWidget)
        self.setLayout(vlayout)

        # connect signals
        self.directoryLineEdit.textChanged.connect(self.set_directory)
        self.directoryPickerButton.clicked.connect(self.choose_directory)
        self.refreshButton.clicked.connect(self.edfTab.table.refresh)
        self.refreshButton.clicked.connect(self.nxsTab.tableWidget.refresh)
        # self.edfTab.edfTreatmentClicked.connect(self.on_treatment)
        self.edfTab.treatButton.clicked.connect(self.on_treatment_clicked)

    def on_treatment_clicked(self):
        model = self.nxsTab.treeWidget.treeview.findHdf5TreeModel()
        model.clear()
        self.edfTab.treat()
        self.nxsTab.tableWidget.refresh()

    def set_directory(self):
        text = self.directoryLineEdit.text()
        self.edfTab.table.setDirectory(text)
        self.nxsTab.set_directory(text)

    def choose_directory(self):
        basedir = os.path.expanduser("~")
        fname = qt.QFileDialog.getExistingDirectory(self, 'Select data directory', basedir,
                                                    options=qt.QFileDialog.DontUseNativeDialog)
        if fname:
            self.directoryLineEdit.setText(fname)
            # self.edfTab.table.setDirectory(fname)
            # self.nxsTab.setDirectory(fname)


class FunctionWorker(qt.QObject):
    progress = qt.pyqtSignal()
    finished = qt.pyqtSignal()

    def __init__(self, cmdList, fileList):
        super(FunctionWorker, self).__init__()
        self.cmdList = cmdList
        self.fileList = fileList

    def run(self):
        self.progress.emit()
        for file in self.fileList:
            for script in self.cmdList:
                cmd = 'xeuss.' + script.replace('root', '\'' + file.replace('\\', '/') + '\'')
                print(cmd)
                eval(cmd)
        self.finished.emit()


class SaxsUtily(qt.QMainWindow):
    """
    This window show an example of use of a Hdf5TreeView.

    The tree is initialized with a list of filenames. A panel allow to play
    with internal property configuration of the widget, and a text screen
    allow to display events.
    """

    def __init__(self):
        """

        """

        qt.QMainWindow.__init__(self)
        self.setWindowTitle("pygdatax GUI")
        self.setWindowIcon(getQIcon('logo_llb.ico'))

        self.__asyncload = False
        central_wigdet = self.centralWidget()

        self.fileSurvey = FileSurvey()
        # self.nxsFileTable = NexusFileTable()

        # plot widget
        self.plotWindow = DataView()
        # treatment dock widget
        self.treatmentDock = qt.QDockWidget('treatment', self)
        # self.treatmentDock.setStyleSheet("border: 5px solid black")
        self.treatmentDock.setFeatures(qt.QDockWidget.DockWidgetFloatable |
                                       qt.QDockWidget.DockWidgetMovable)
        self.editor = CommandTreatmentWidget(self, module=xeuss)
        self.treatmentDock.setWidget(self.editor)
        self.treatmentDock.setFloating(False)
        # replace the addTabbedwidget metho of the plot window
        self.plotWindow._dockWidgets.append(self.treatmentDock)
        self.plotWindow.addDockWidget(qt.Qt.BottomDockWidgetArea,self.treatmentDock)
        # self.treatmentDock.setAllowedAreas(qt.Qt.BottomDockWidgetArea)
        # self.addDockWidget(qt.Qt.RightDockWidgetArea, self.treatmentDock)
        self.treatmentDock.show()
        # directory picker layout
        spliter = qt.QSplitter(qt.Qt.Horizontal)
        spliter.addWidget(self.fileSurvey)
        spliter.addWidget(self.plotWindow)
        spliter.setStretchFactor(1, 1)
        main_panel = qt.QWidget(self)
        layout = qt.QVBoxLayout()
        layout.addWidget(spliter)
        layout.setStretchFactor(spliter, 1)
        main_panel.setLayout(layout)
        self.setCentralWidget(main_panel)


        # connect signals
        # edf table dispplay
        self.fileSurvey.edfTab.edfSelectionChanged.connect(self.displayEdf)
        # nxs table dispplay
        self.fileSurvey.nxsTab.nxsSelectionChanged.connect(self.displayNxs)
        # treatement run
        self.editor.runClicked.connect(self.run_function)
        self.fileSurvey.nxsTab.treeWidget.selectedNodeChanged.connect(self.displayNxTree)

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
        roi.setSelectable(True)
        roi.setEditable(True)

    def run_functionTest(self, cmdList):
        selectedFiles = self.fileSurvey.nxsTab.tableWidget.get_selectedFiles()
        model = self.fileSurvey.nxsTab.treeWidget.treeview.findHdf5TreeModel()
        model.clear()
        self.thread = qt.QThread()
        self.worker = FunctionWorker(cmdList, selectedFiles)
        self.worker.moveToThread(self.thread)

        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(lambda: self.editor.run_btn.setEnabled(False))
        # when thread finishes
        self.worker.finished.connect(lambda: self.editor.run_btn.setEnabled(True))
        self.worker.finished.connect(self.worker.deleteLater)
        self.worker.finished.connect(self.thread.deleteLater)
        self.thread.finished.connect(self.fileSurvey.nxsTab.tableWidget.on_selectionChanged)
        self.thread.start()

    def run_function(self, cmdList):
        selectedFiles = self.fileSurvey.nxsTab.tableWidget.get_selectedFiles()
        model = self.fileSurvey.nxsTab.treeWidget.treeview.findHdf5TreeModel()
        model.clear()
        for file in selectedFiles:
            for script in cmdList:
                for line in script.splitlines():
                    cmd = 'xeuss.' + line.replace('root', '\'' + file.replace('\\', '/') + '\'')
                    try:
                        eval(cmd)
                        # print(cmd)
                    except:
                        print('command : '+cmd+'not performed on:' + file)
            # model.insertFile(file)
        self.fileSurvey.nxsTab.treeWidget.operationPerformed.emit()
        self.fileSurvey.nxsTab.tableWidget.on_selectionChanged()

    def displayEdf(self, file):
        """
        Called to update the dataviewer with the selected data.
        """
        self.plotWindow.clear()
        self.plotWindow.setKeepDataAspectRatio(True)
        self.plotWindow.setXAxisLogarithmic(False)
        self.plotWindow.setYAxisLogarithmic(False)
        if file.endswith('.edf'):
            dataObj = fabio.open(file)
            data = dataObj.data
            self.plotWindow.clear()
            if 'Comments' in dataObj.header:
                legend = dataObj.header['Comments']
            else:
                legend = dataObj.filename
            self.plotWindow.addImage(data, replace=True,
                                       legend=legend, xlabel='x [pixel]',
                                       ylabel='y [pixel]')

    def displayNxs(self, files):
        self.plotWindow.clear()
        c = ['blue', 'red', 'green', 'black', 'yellow', 'grey', 'magenta', 'cyan',
             'darkGreen', 'darkBrown', 'darkCyan', 'darkYellow', 'darkMagenta']

        for i, file in enumerate(files):
            h5file = sio.open(file)
            if is_NXroot_with_default_NXdata(h5file):
                nxd = get_default(h5file)
                legend = os.path.basename(file).split('.')[0] +'/'
                # legend += h5file['entry0/sample/sample_name'].asstr()[()]
                try:
                    legend += h5file['entry0/sample/sample_name'][()].decode()
                except AttributeError:
                    legend += h5file['entry0/sample/sample_name'][()]
                if nxd.is_curve:
                    xlabel = nxd.axes_names[0]
                    ylabel = nxd.signal_name
                    if 'units' in nxd.axes[0].attrs:
                        xlabel += ' [' + nxd.axes[0].attrs['units'] + ']'
                    if 'units' in nxd.signal.attrs:
                        ylabel += ' [' + nxd.signal.attrs['units'] + ']'
                    try:
                        col = colors.COLORDICT[c[i]]
                    except IndexError:
                        col = c[randint(len(c))]
                    self.plotWindow.addCurve(nxd.axes[0], nxd.signal,
                                               yerror=nxd.errors,
                                               legend=legend,
                                               replace=False,
                                               color=col,
                                               xlabel=xlabel,
                                               ylabel=ylabel,
                                               resetzoom=True)
                    self.plotWindow.setGraphXLabel(xlabel)
                    self.plotWindow.setGraphYLabel(ylabel)
                    self.plotWindow.setKeepDataAspectRatio(False)
                elif nxd.is_image:
                    if nxd.axes_names == [None, None]:
                        origin = (0., 0.)
                        scale = (1., 1.)
                        xlabel = 'x [pixel]'
                        ylabel = 'y [pixel]'
                        self.plotWindow.setKeepDataAspectRatio(True)
                        self.plotWindow.setXAxisLogarithmic(False)
                        self.plotWindow.setYAxisLogarithmic(False)
                    else:
                        aspect_button = self.plotWindow.getKeepDataAspectRatioButton()
                        self.plotWindow.setKeepDataAspectRatio(False)
                        origin = (nxd.axes[1][0], nxd.axes[0][1])
                        scale_x = np.abs(nxd.axes[1][0] - nxd.axes[1][-1]) / len(nxd.axes[1])
                        scale_y = np.abs(nxd.axes[0][0] - nxd.axes[0][-1]) / len(nxd.axes[0])
                        scale = (scale_x, scale_y)
                        xlabel = nxd.axes_names[1]
                        ylabel = nxd.axes_names[0]
                        if 'units' in nxd.axes[1].attrs:
                            xlabel += ' [' + nxd.axes[1].attrs['units'] + ']'
                        if 'units' in nxd.axes[0].attrs:
                            ylabel += ' [' + nxd.axes[0].attrs['units'] + ']'

                    self.plotWindow.addImage(nxd.signal, replace=True,
                                               legend=legend, xlabel='x',
                                               ylabel='y',
                                               origin=origin,
                                               scale=scale)
                    self.plotWindow.setGraphXLabel(xlabel)
                    self.plotWindow.setGraphYLabel(ylabel)
                else:
                    return

    def displayNxTree(self, nodes):
        self.plotWindow.clear()
        c = ['blue', 'red', 'green', 'black', 'yellow', 'grey', 'magenta', 'cyan',
             'darkGreen', 'darkBrown', 'darkCyan', 'darkYellow', 'darkMagenta']

        for i, s in enumerate(nodes):
            if is_group(s.h5py_object) or is_file(s.h5py_object):
                nxd = get_default(s.h5py_object)
                if nxd is None:
                    return
                elif not nxd.is_valid:
                    return
                legend = os.path.split(s.h5py_object.file.filename)[1] + s.h5py_object.name
                if nxd.is_curve:
                    xlabel = nxd.axes_names[0]
                    ylabel = nxd.signal_name
                    if 'units' in nxd.axes[0].attrs:
                        xlabel += ' [' + nxd.axes[0].attrs['units'] + ']'
                    if 'units' in nxd.signal.attrs:
                        ylabel += ' [' + nxd.signal.attrs['units'] + ']'
                    self.plotWindow.addCurve(nxd.axes[0], nxd.signal,
                                               yerror=nxd.errors,
                                               legend=legend,
                                               replace=False,
                                               color=colors.COLORDICT[c[i]],
                                               xlabel=xlabel,
                                               ylabel=ylabel,
                                               resetzoom=True)
                    self.plotWindow.setGraphXLabel(xlabel)
                    self.plotWindow.setGraphYLabel(ylabel)
                    self.plotWindow.setKeepDataAspectRatio(False)
                elif nxd.is_image:
                    if nxd.axes_names == [None, None]:
                        origin = (0., 0.)
                        scale = (1., 1.)
                        xlabel = 'x [pixel]'
                        ylabel = 'y [pixel]'
                        self.plotWindow.setKeepDataAspectRatio(True)
                        self.plotWindow.setXAxisLogarithmic(False)
                        self.plotWindow.setYAxisLogarithmic(False)
                    else:
                        aspect_button = self.plotWindow.getKeepDataAspectRatioButton()
                        self.plotWindow.setKeepDataAspectRatio(False)
                        origin = (nxd.axes[1][0], nxd.axes[0][1])
                        scale_x = np.abs(nxd.axes[1][0] - nxd.axes[1][-1]) / len(nxd.axes[1])
                        scale_y = np.abs(nxd.axes[0][0] - nxd.axes[0][-1]) / len(nxd.axes[0])
                        scale = (scale_x, scale_y)
                        xlabel = nxd.axes_names[1]
                        ylabel = nxd.axes_names[0]
                        if 'units' in nxd.axes[1].attrs:
                            xlabel += ' [' + nxd.axes[1].attrs['units'] + ']'
                        if 'units' in nxd.axes[0].attrs:
                            ylabel += ' [' + nxd.axes[0].attrs['units'] + ']'

                    self.plotWindow.addImage(nxd.signal, replace=True,
                                               legend=legend, xlabel='x',
                                               ylabel='y',
                                               origin=origin,
                                               scale=scale)
                    self.plotWindow.setGraphXLabel(xlabel)
                    self.plotWindow.setGraphYLabel(ylabel)
                else:
                    return
            elif is_dataset(s.h5py_object):
                legend = os.path.split(s.h5py_object.file.filename)[1] + s.h5py_object.name
                if len(s.h5py_object.shape) == 2:
                    origin = (0., 0.)
                    scale = (1., 1.)
                    xlabel = 'x [pixel]'
                    ylabel = 'y [pixel]'
                    self.plotWindow.setKeepDataAspectRatio(True)
                    self.plotWindow.addImage(s.h5py_object, replace=True,
                                               legend=legend, xlabel='x',
                                               ylabel='y',
                                               origin=origin,
                                               scale=scale)
                    self.plotWindow.setGraphXLabel(xlabel)
                    self.plotWindow.setGraphYLabel(ylabel)
                elif len(s.h5py_object.shape) == 1:
                    xlabel = 'index'
                    ylabel = s.h5py_object.name
                    if 'units' in s.h5py_object.attrs:
                        ylabel += ' [' + str(s.h5py_object.attrs['units']) + ']'
                    x = np.arange(len(s.h5py_object)) + 1
                    self.plotWindow.addCurve(x, s.h5py_object,
                                             legend=legend,
                                             replace=False,
                                             color=colors.COLORDICT[c[i]],
                                             xlabel=xlabel,
                                             ylabel=ylabel,
                                             resetzoom=True)
                    self.plotWindow.setGraphXLabel(xlabel)
                    self.plotWindow.setGraphYLabel(ylabel)
                    self.plotWindow.setKeepDataAspectRatio(False)
                else:
                    return



class NexusFileTable(qt.QTableWidget):
    directory = ''
    file_extension = '.nxs'
    fileSelectedChanged = qt.pyqtSignal(list)

    def __init__(self, parent=None):
        super(NexusFileTable, self).__init__(parent=None)
        self.setColumnCount(6)
        self.setRowCount(3)
        self.setRowCount(4)
        self.setSelectionBehavior(qt.QAbstractItemView.SelectRows | qt.QAbstractItemView.MultiSelection)
        self.setHorizontalHeaderLabels(['File', 'comment', 'e (mm)','tr', 'distance',  'counting time', ' date'])
        self.itemSelectionChanged.connect(self.on_selectionChanged)

    def setDirectory(self, directory):
        forlderPath = Path(directory)
        if forlderPath.is_dir():
            self.directory = forlderPath
            self.refresh()

    def refresh(self):
        self.itemSelectionChanged.disconnect()
        if os.path.isdir(self.directory):
            l = os.listdir(self.directory)
            # l.sort()
            fileList = []
            for item in l:
                if os.path.splitext(item)[1] == self.file_extension:
                    fileList.append(item)
            # self.clearContents()
            self.setRowCount(len(fileList))
            for i, file in enumerate(fileList):
                description = get_nxs_description(os.path.join(self.directory, file))
                for j, des in enumerate(description):
                    item = qt.QTableWidgetItem(des)
                    item.setFlags(qt.Qt.ItemIsSelectable | qt.Qt.ItemIsEnabled)
                    self.setItem(i, j, item)

            self.sortItems(0, qt.Qt.AscendingOrder)
            self.itemSelectionChanged.connect(self.on_selectionChanged)

    def on_selectionChanged(self):
        items = self.selectedItems()
        selectedFiles = []
        for item in items:
            row = item.row()
            file = self.item(row, 0).text()
            selectedFiles.append(os.path.join(self.directory, file))
        self.fileSelectedChanged.emit(selectedFiles)

    def get_selectedFiles(self):
        items = self.selectedItems()
        selectedFiles = []
        for item in items:
            row = item.row()
            file = self.item(row, 0).text()
            selectedFiles.append(os.path.join(self.directory, file))
        return selectedFiles


class NexusTreeWidget(qt.QWidget):
    operationPerformed = qt.pyqtSignal()
    selectedNodeChanged = qt.pyqtSignal(list)

    def __init__(self):
        super(NexusTreeWidget, self).__init__()
        self.clear_last_btn = qt.QPushButton('clear last')
        self.clear_last_btn.setToolTip('Clear last performed calculation')
        self.clear_all_btn = qt.QPushButton('clear all')
        self.clear_all_btn.setToolTip('clear all performed calculations')
        # self.sync_btn = qt.QPushButton('sync')
        # self.sync_btn.setToolTip('synchronize the .nxs files')
        """Silx HDF5 TreeView"""
        self.treeview = silx.gui.hdf5.Hdf5TreeView(self)
        treemodel = silx.gui.hdf5.Hdf5TreeModel(self.treeview,
                                                ownFiles=True
                                                )
        # treemodel.sigH5pyObjectLoaded.connect(self.__h5FileLoaded)
        # treemodel.sigH5pyObjectRemoved.connect(self.__h5FileRemoved)
        # treemodel.sigH5pyObjectSynchronized.connect(self.__h5FileSynchonized)
        treemodel.setDatasetDragEnabled(False)
        # self.treeview.setModel(treemodel)
        self.__treeModelSorted = silx.gui.hdf5.NexusSortFilterProxyModel(self.treeview)
        self.__treeModelSorted.setSourceModel(treemodel)
        self.__treeModelSorted.sort(0, qt.Qt.AscendingOrder)
        self.__treeModelSorted.setSortCaseSensitivity(qt.Qt.CaseInsensitive)
        self.treeview.setModel(self.__treeModelSorted)
        self.treeview.setSelectionMode(qt.QAbstractItemView.ExtendedSelection)
        # layout
        hlayout = qt.QHBoxLayout()
        hlayout.addWidget(self.clear_last_btn)
        hlayout.addWidget(self.clear_all_btn)
        # hlayout.addWidget(self.sync_btn)
        vlayout = qt.QVBoxLayout()
        vlayout.addLayout(hlayout)
        vlayout.addWidget(self.treeview)
        self.setLayout(vlayout)

        # connect signals
        # self.sync_btn.clicked.connect(self.sync_all)
        self.clear_last_btn.clicked.connect(self.clear_last)
        self.clear_all_btn.clicked.connect(self.clear_all)
        self.treeview.selectionModel().selectionChanged.connect(self.on_tree_selection)

    def load_files(self, files):
        model = self.treeview.findHdf5TreeModel()
        model.clear()
        for file in files:
            model.insertFile(file, row=-1)
        self.treeview.expandToDepth(0)

    def clear_last(self):
        model = self.treeview.findHdf5TreeModel()
        nrow = model.rowCount()
        files = []
        for n in range(nrow):
            index = model.index(n, 0, qt.QModelIndex())
            node = model.nodeFromIndex(index)
            filename = node.obj.filename
            model.removeH5pyObject(node.obj)
            # node.obj.close()
            root = nxlib.loadfile(filename, mode='rw')
            nxlib.delete_last_entry(root)
            root.close()
            model.insertFile(filename, row=n)
        self.treeview.expandToDepth(0)
        self.operationPerformed.emit()

    def clear_all(self):
        model = self.treeview.findHdf5TreeModel()
        nrow = model.rowCount()
        for n in range(nrow):
            index = model.index(n, 0, qt.QModelIndex())
            node = model.nodeFromIndex(index)
            filename = node.obj.filename
            model.removeH5pyObject(node.obj)
            root = nxlib.loadfile(filename, mode='rw')
            nxlib.delete_all_entry(root)
            root.close()
            model.insertFile(filename, row=n)
        self.treeview.expandToDepth(0)
        self.operationPerformed.emit()

    def sync_all(self):
        model = self.treeview.findHdf5TreeModel()
        nrow = model.rowCount()

        for n in range(nrow):
            index = model.index(n, 0, qt.QModelIndex())
            node = model.nodeFromIndex(index)
            filename = node.obj.filename
            model.removeH5pyObject(node.obj)
            root = nxlib.loadfile(filename, mode='rw')
            nxlib.delete_all_entry(root)
            root.close()
            model.insertFile(filename, row=n)
        self.treeview.expandToDepth(0)
        self.operationPerformed.emit()

    def on_tree_selection(self):
        selected = list(self.treeview.selectedH5Nodes())
        self.selectedNodeChanged.emit(selected)


class NexusTreatmentWidget(qt.QWidget):
    nxsSelectionChanged = qt.pyqtSignal(list)

    def __init__(self):
        super(NexusTreatmentWidget, self).__init__()
        self.tableWidget = NexusFileTable()
        self.treeWidget = NexusTreeWidget()
        # layout
        spliter = qt.QSplitter(qt.Qt.Vertical)
        spliter.addWidget(self.tableWidget)
        spliter.addWidget(self.treeWidget)
        # spliter.setStretchFactor(0, 3)
        layout = qt.QVBoxLayout()
        layout.addWidget(spliter)
        self.setLayout(layout)
        # connect
        self.tableWidget.fileSelectedChanged.connect(self.on_file_selected)
        self.treeWidget.operationPerformed.connect(self.on_tree_operation)
        # context menu for quick functions
        self.tableWidget.setContextMenuPolicy(qt.Qt.CustomContextMenu)
        self.tableWidget.customContextMenuRequested.connect(self.generateMenu)

    def set_directory(self, directory):
        self.tableWidget.setDirectory(directory)

    def on_file_selected(self, files):
        self.treeWidget.load_files(files)
        self.nxsSelectionChanged.emit(files)

    def on_tree_operation(self):
        self.tableWidget.on_selectionChanged()

    def generateMenu(self, event):
        menu = qt.QMenu()
        concatAction = qt.QAction(getQIcon('concat.ico'), 'concat')
        concatAction.triggered.connect(self._concat)
        convertAction = qt.QAction(getQIcon('nxs2text.ico'), 'convert to .txt')
        convertAction.triggered.connect(self._convert2txt)
        getPathAction = qt.QAction(getQIcon('clipboard.ico'), 'copy file path')
        getPathAction.triggered.connect(self._copyPath2clipboard)
        menu.addAction(concatAction)
        menu.addAction(convertAction)
        menu.addAction(getPathAction)
        action = menu.exec_(self.mapToGlobal(event))

    def _copyPath2clipboard(self):
        items = self.tableWidget.selectedItems()
        fileList = []
        for item in items:
            row = item.row()
            # fileList.append('"'+os.path.join(self.tableWidget.directory, self.tableWidget.item(row, 0).text())+'"')
            clip = '"' + Path(self.tableWidget.directory, self.tableWidget.item(row, 0).text()).__str__()+'"'
            fileList.append(clip)
        if fileList:
            qt.QApplication.clipboard().setText(fileList[0])

    def _concat(self, items):
        items = self.tableWidget.selectedItems()
        fileList = []
        for item in items:
            row = item.row()
            fileList.append(os.path.join(self.tableWidget.directory, self.tableWidget.item(row, 0).text()))
        treemodel = self.treeWidget.treeview.findHdf5TreeModel()
        # clear model before using the filess
        treemodel.clear()
        if len(fileList) > 1:
            firstFile = fileList.pop(0)
            for file in fileList:
                xeuss.concat(firstFile, file=file)
        treemodel.insertFile(firstFile)
        for file in fileList:
            treemodel.insertFile(file)

    def _convert2txt(self):
        items = self.tableWidget.selectedItems()
        fileList = []
        for item in items:
            row = item.row()
            fileList.append(os.path.join(self.tableWidget.directory, self.tableWidget.item(row, 0).text()))
        treemodel = self.treeWidget.treeview.findHdf5TreeModel()
        # clear model before using the filess
        treemodel.clear()
        for file in fileList:
            xeuss.save_as_txt(file)
        for file in fileList:
            treemodel.insertFile(file)


class CommandTreatmentWidget(qt.QWidget):
    runClicked = qt.pyqtSignal(list)

    def __init__(self, parent, module=None):
        super(CommandTreatmentWidget, self).__init__()
        self.run_btn = qt.QPushButton('run')
        self.runAll_btn = qt.QPushButton('run all')
        self.run_btn.clicked.connect(self.run)
        self.runAll_btn.clicked.connect(self.runAll)
        self.tabWidget = qt.QTabWidget()
        self.add_btn = qt.QPushButton('+')
        self.add_btn.setFixedSize(20,20)
        self.add_btn.clicked.connect(self.addTab)
        self.tabWidget.setCornerWidget(self.add_btn, corner=qt.Qt.TopLeftCorner)
        self.tabWidget.setTabsClosable(True)
        self.tabWidget.setMaximumHeight(500)
        self.tabWidget.tabCloseRequested.connect(self.closeTabs)
        if module is not None:
            self.commandList = moduledescription.get_commandList(module)
        else:
            self.commandList = []
        self.tabWidget.addTab(MultiLineCodeEditor(parent=self, completerList=self.commandList), 'cmd1')
        hlayout = qt.QHBoxLayout()
        layout = qt.QVBoxLayout(self)
        hlayout.addWidget(self.run_btn)
        hlayout.addWidget(self.runAll_btn)
        hlayout.addStretch()
        layout.addLayout(hlayout)
        layout.addWidget(self.tabWidget)
        layout.addStretch()
        # self.formLayout = qt.QFormLayout(self)
        # layout.addLayout(self.formLayout)
        self.setLayout(layout)

    def closeTabs(self, index):
        count = self.tabWidget.count()
        if count > 1:
            self.tabWidget.removeTab(index)

    def addTab(self):
        count = self.tabWidget.count()
        # widget = qt.QWidget()
        # layout = qt.QVBoxLayout()
        # layout.addWidget(CodeEditor(self))
        # layout.addStretch()
        # widget.setLayout(layout)
        self.tabWidget.addTab(MultiLineCodeEditor(parent=self, completerList=self.commandList), 'cmd'+str(count+1))
        # self.tabWidget.addTab(widget, 'cmd'+str(count+1))

    def run(self):
        widget = self.tabWidget.currentWidget()
        text = widget.toPlainText()
        self.runClicked.emit([text])

    def runAll(self):
        count = self.tabWidget.count()
        l = []
        for i in range(count):
            widget = self.tabWidget.widget(i)
            l.append(widget.toPlainText())
        self.runClicked.emit(l)

    def on_comboBox(self, text):
        widget = self.tabWidget.currentWidget()
        widget.setText(text)


class CodeEditor(qt.QLineEdit):
    """
    QLineEdit widget with treatment function autocompletion
    """
    def __init__(self, parent=None, completerList=[]):
        super().__init__()
        # self.setTabStopDistance(
            # qt.QFontMetricsF(self.font()).horizontalAdvance(' ') * 4)
        # self.highlighter = PythonHighlighter(self.document())

        completer = qt.QCompleter(completerList)
        self.setCompleter(completer)
        # self.setFixedHeight(30)
        self.setContextMenuPolicy(qt.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.generateMenu)
        self.completerList = completerList

    def generateMenu(self, event):
        menu = qt.QMenu()
        for fun in moduledescription.get_commandList(xeuss):
            menu.addAction(fun, self.actionClicked)
        menu.exec_(self.mapToGlobal(event))

    def actionClicked(self):
        action = self.sender()
        completer = self.completer()
        # print(self.completer.model())
        self.setText(action.text())


class MyCompleter(qt.QCompleter):
    insertText = qt.pyqtSignal(str)
    def __init__(self, parent=None, completerList=[]):
        super(MyCompleter, self).__init__(completerList, parent)
        self.setCompletionMode(qt.QCompleter.PopupCompletion)
        self.highlighted.connect(self.setHighlighted)

    def setHighlighted(self, text):
        self.lastSelected = text

    def getSelected(self):
        return self.lastSelected


class MultiLineCodeEditor(qt.QPlainTextEdit):
    def __init__(self, parent=None, completerList=[]):
        super(MultiLineCodeEditor, self).__init__(parent)
        self.completer = MyCompleter(completerList=completerList)
        self.completer.setWidget(self)
        self.completer.insertText.connect(self.insertCompletion)
        self.setContextMenuPolicy(qt.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.generateMenu)
        self.completerList = completerList
        self.setMinimumHeight(25)
        self.highlighter = qt.Qt

    def insertCompletion(self, completion):
        tc = self.textCursor()
        extra = (len(completion) - len(self.completer.completionPrefix()))
        tc.movePosition(qt.QTextCursor.Left)
        tc.movePosition(qt.QTextCursor.EndOfWord)
        tc.insertText(completion[-extra:])
        self.setTextCursor(tc)
        self.completer.popup().hide()

    def focusInEvent(self, event):
        if self.completer:
            self.completer.setWidget(self)
        qt.QPlainTextEdit.focusInEvent(self, event)

    def keyPressEvent(self, event):
        tc = self.textCursor()
        if event.key() == qt.Qt.Key_Tab and self.completer.popup().isVisible():
            self.completer.insertText.emit(self.completer.getSelected())
            self.completer.setCompletionMode(qt.QCompleter.PopupCompletion)
            return

        qt.QPlainTextEdit.keyPressEvent(self, event)
        tc.select(qt.QTextCursor.WordUnderCursor)
        cr = self.cursorRect()

        if len(tc.selectedText()) > 0:
            self.completer.setCompletionPrefix(tc.selectedText())
            popup = self.completer.popup()
            popup.setCurrentIndex(self.completer.completionModel().index(0, 0))

            cr.setWidth(self.completer.popup().sizeHintForColumn(0)
                        + self.completer.popup().verticalScrollBar().sizeHint().width())
            self.completer.complete(cr)
        else:
            self.completer.popup().hide()

    def generateMenu(self, event):
        menu = qt.QMenu()
        for fun in self.completerList:
            menu.addAction(fun, self.actionClicked)
        menu.exec_(self.mapToGlobal(event))

    def actionClicked(self):
        action = self.sender()
        text = action.text()
        tc = self.textCursor()
        # extra = len(text)
        # tc.movePosition(qt.QTextCursor.Left)
        # tc.movePosition(qt.QTextCursor.EndOfWord)
        # tc.insertText(completion[-extra:])
        tc.insertText(text)
        tc.movePosition(qt.QTextCursor.EndOfWord)
        self.setTextCursor(tc)


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
        # self.plotWindow = PlotWindow(backend=None, resetzoom=True,
        #                              autoScale=True, logScale=True,
        #                              grid=False, curveStyle=True, colormap=True,
        #                              aspectRatio=True, yInverted=True,
        #                              copy=True, save=True, print_=True,
        #                              control=True, position=posInfo,
        #                              roi=False, mask=True, fit=True)
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
        # roi table widget
        # self.roiDock = qt.QDockWidget('rois', self)
        # # self.treatmentDock.setStyleSheet("border: 5px solid black")
        # self.roiDock.setFeatures(qt.QDockWidget.DockWidgetFloatable |
        #                          qt.QDockWidget.DockWidgetMovable)
        # self.roiTableWidget = RegionOfInterestTableWidget(self)
        # self.roiTableWidget.setRegionOfInterestManager(self.roiManager)
        # self.roiDock.setWidget(self.roiTableWidget)
        # self.roiDock.setFloating(False)
        # replace the addTabbedwidget metho of the plot window
        # self._dockWidgets.append(self.roiDock)
        # self.addDockWidget(qt.Qt.BottomDockWidgetArea, self.roiDock)

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
# TODO : finish it
class DataView3Dectectors(PlotWindow):

    def __init__(self):
        super().__init__(backend=None, resetzoom=True,
                         autoScale=True, logScale=True,
                         grid=False, curveStyle=True, colormap=True,
                         aspectRatio=True, yInverted=True,
                         copy=True, save=True, print_=True,
                         control=True, position=[('X', lambda x, y: x),
                                                 ('Y', lambda x, y: y),
                                                 ('Data', self._zValue)],
                         roi=False, mask=True, fit=True)
        layout = qt.QGridLayout()
        layout.addWidget(self.getWidgetHandle(), 0, 1)
        self.left_detector = PlotWidget(parent=self)
        self.left_detector.setMaximumWidth(200)
        self.bot_detector = PlotWidget(parent=self)
        self.bot_detector.setMaximumHeight(200)
        layout.addWidget(self.left_detector, 0, 0, 2, 1)
        layout.addWidget(self.bot_detector, 1, 1)
        layout.setColumnStretch(0, 1)
        layout.setColumnStretch(1, 0)
        layout.setRowStretch(0, 1)
        layout.setRowStretch(1, 0)
        centralWidget = qt.QWidget(self)
        centralWidget.setLayout(layout)
        self.setCentralWidget(centralWidget)

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


class TreatmentWidget(qt.QWidget):
    stackCommandClicked = qt.pyqtSignal(str)
    runClicked = qt.pyqtSignal(str)

    def __init__(self, parent=None):
        super(TreatmentWidget, self).__init__(parent=parent)
        self.functionComboBox = qt.QComboBox(self)
        self.descriptionDict = {}
        self.functionComboBox.currentIndexChanged.connect(self.on_function_selected)

    def setModule(self, module):
        self.descriptionDict = moduledescription.get_descriptionDict(module, decorator='@nxlib.treatment_function')
        self.functionComboBox.clear()
        self.functionComboBox.addItems(self.descriptionDict.keys())

    def on_function_selected(self,i):
        self.findChildren()


class ParametersWidget(qt.QWidget):

    def __init__(self, function_description, parent=None):
        super(ParametersWidget, self).__init__(parent=parent)
        self.function_description = function_description
        qt.QFormLayout(self)
        fo




def main():
    # unlock hdf5 files for file access during plotting
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    warnings.filterwarnings("ignore", category=mplDeprecation)
    app = qt.QApplication([])
    # sys.excepthook = qt.exceptionHandler
    window = SaxsUtily()
    window.show()
    result = app.exec_()
    # remove ending warnings relative to QTimer
    app.deleteLater()
    sys.exit(result)



if __name__ == "__main__":
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    app = qt.QApplication([])
    # splash_pix = qt.QPixmap('/home/achennev/python/pygdatax/src/pygdatax/resources/empty_cell.png')
    #
    # splash = qt.QSplashScreen(splash_pix, qt.Qt.WindowStaysOnTopHint)
    # splash.setMask(splash_pix.mask())
    # splash.show()
    # app.processEvents()
    #
    warnings.filterwarnings("ignore", category=mplDeprecation)
    ##################################################################
    # window = SaxsUtily()
    # window.show()
    # folder = '/home/achennev/Bureau/PIL pour tiago/PIL NP/2021-07-21_TOC'
    # window.fileSurvey.directoryLineEdit.setText(folder)
    ####################################################################
    w = DataView3Dectectors()
    w.show()
    # splash.finish(window)
    result = app.exec_()
    app.deleteLater()
    sys.exit(result)

