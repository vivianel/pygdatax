import sys
import os
import fabio
import silx.gui.hdf5
from silx.gui import qt, colors
from silx.gui.plot import PlotWindow, Profile
import silx.io as sio
from silx.io.nxdata import is_NXroot_with_default_NXdata, get_default
from pygdatax import xeuss
from pygdatax import nxlib
import nexusformat.nexus as nx
from numpy.random import randint
from pygdatax.icons import getQIcon

COMPLETER_NAMES = ['azimutal_integration(root,x0=None,y0=None,mask=None,bins=900)',
                   'azimutal_integration2D(root, mask=None, x0=None, y0=None, distance=None,r_bins=900, chi_bins=360',
                   'bkg_substraction(root, bkg=0)',
                   'concat(root, file=None)',
                   'cut(root, xmin=None, xmax=None)',
                   'normalization_factor(root, factor=None)',
                   'polar_cut(root, q=None, pixel_width=1)',
                   'q_scale(root, distance=None)',
                   'ref_substraction(root, ref_file=None, prefactor=0)',
                   'resu(root, dark_file=None, ec_file=None, thickness=None)',
                   'azimutal_integration2D(root, mask=None, x0=None, y0=None, distance=None,r_bins=900, chi_bins=360)',
                   'save_as_txt(root)'
                   'set_beam_center(root,x0=None,y0=None)'
                   ]


def get_edf_description(edf_filepath):
    des = []
    if os.path.isfile(edf_filepath):
        if edf_filepath.endswith('.edf'):
            dataObj = fabio.open(edf_filepath)
            try:
                des.append(dataObj.header['title'])
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
    des = []
    if os.path.isfile(nxs_filepath):
        if nxs_filepath.endswith('.nxs'):
            root = nx.nxload(nxs_filepath, mode='r')
            entry = root[nxlib.get_last_entry_key(root)]
            try:
                des.append(os.path.basename(nxs_filepath))
                des.append(str(entry.sample.sample_name.nxdata.decode()))
                des.append(str(str(entry.sample.thickness.nxdata)))
                des.append(str(entry.sample.transmission.nxdata))
                des.append(str(entry.instrument.detector.distance.nxdata))
                des.append(str(entry.sample.count_time.nxdata))
                des.append(root.attrs['file_time'])
            except KeyError:
                des.append(os.path.split(nxs_filepath)[1])
                des += 5 * ['']
            # compatibility with windows
            except AttributeError:
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
        if os.path.isdir(directory):
            self.directory = directory
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
                    # self.item(i, 0).setIcon(qt.QIcon(qt.QPixmap('../ressources/dark.ico')))
                    self.item(i, 0).setIcon(getQIcon('dark.ico'))
                elif filepath == self.emptyCellFile:
                    self.set_row_bkg(i, qt.QColor("cyan"))
                    self.item(i, 0).setIcon(getQIcon('empty_cell.ico'))
                    self.item(i, 0).setIcon(getQIcon('empty_cell.ico'))
                elif filepath == self.emptyBeamFile:
                    self.set_row_bkg(i, qt.QColor("red"))
                    self.item(i, 0).setIcon(qt.QIcon(qt.QPixmap('../ressources/beam.ico')))
                elif filepath in self.trashFiles:
                    self.set_row_bkg(i, qt.QColor("grey"))
                    self.item(i, 0).setIcon(qt.QIcon(qt.QPixmap('../ressources/cross.ico')))
                elif filepath == self.maskFile:
                    self.set_row_bkg(i, qt.QColor("white"))
                    self.item(i, 0).setIcon(qt.QIcon(qt.QPixmap('../ressources/mask.ico')))
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
        trashAction = qt.QAction(getQIcon('cross.ico'), 'trash')
        trashAction.triggered.connect(self._set_trash)
        sampleAction = qt.QAction('sample')
        sampleAction.triggered.connect(self._set_sample)
        maskAction = qt.QAction(getQIcon('mask.ico'), 'mask')
        maskAction.triggered.connect(self._set_mask)
        # build menu
        menu.addAction(darkAction)
        menu.addAction(emptyCellAction)
        menu.addAction(emptyBeamAction)
        menu.addAction(trashAction)
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
            row = current_item.row()
            ncol = self.columnCount()
            first_col_item = self.item(row, 0)
            file = first_col_item.text()
            first_col_item.setIcon(getQIcon('empty_cell.ico'))
            self.set_row_bkg(row, qt.QColor("cyan"))

            current_ec_item = self.findItems(os.path.basename(str(self.emptyCellFile)), qt.Qt.MatchExactly)
            # remove the previous empty cell icons
            if current_ec_item:
                self.set_row_bkg(current_ec_item[0].row(), qt.QColor("white"))
                current_ec_item[0].setIcon(qt.QIcon())
            fullfile = os.path.join(self.directory, file)
            self.emptyCellFile = fullfile
            # remove double reference
            if fullfile == self.darkFile:
                self.darkFile = None
            elif fullfile == self.emptyBeamFile:
                self.emptyBeamFile = None
            elif fullfile == self.maskFile:
                self.maskFile = None
            elif fullfile in self.trashFiles:
                self.trashFiles.remove(fullfile)

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
                # remove double reference
                if fullfile == self.emptyCellFile:
                    self.emptyCellFile = None
                elif fullfile == self.darkFile:
                    self.darkFile = None
                elif fullfile == self.emptyBeamFile:
                    self.emptyBeamFile = None
                elif fullfile == self.maskFile:
                    self.maskFile = None
                elif fullfile in self.trashFiles:
                    self.trashFiles.remove(fullfile)

    def _set_dark(self, event):
        current_item = self.currentItem()
        if current_item is not None:
            row = current_item.row()
            ncol = self.columnCount()
            first_col_item = self.item(row, 0)
            file = first_col_item.text()
            self.set_row_bkg(row, qt.QColor("blue"))
            first_col_item.setIcon(getQIcon('dark.ico'))
            current_dark_item = self.findItems(os.path.basename(str(self.darkFile)), qt.Qt.MatchExactly)
            # remove the previous empty cell icons
            if current_dark_item:
                self.set_row_bkg(current_dark_item[0].row(), qt.QColor("white"))
                current_dark_item[0].setIcon(qt.QIcon())
            # remove double reference
            fullfile = os.path.join(self.directory, file)
            self.darkFile = fullfile
            if fullfile == self.emptyCellFile:
                self.emptyCellFile = None
            elif fullfile == self.emptyBeamFile:
                self.emptyBeamFile = None
            elif fullfile == self.maskFile:
                self.maskFile = None
            elif fullfile in self.trashFiles:
                self.trashFiles.remove(fullfile)
            else:
                pass

    def _set_empty_beam(self):
        current_item = self.currentItem()
        if current_item is not None:
            row = current_item.row()
            ncol = self.columnCount()
            first_col_item = self.item(row, 0)
            file = first_col_item.text()
            self.set_row_bkg(row, qt.QColor("red"))
            first_col_item.setIcon(getQIcon('beam.ico'))

            current_eb_item = self.findItems(os.path.basename(str(self.emptyBeamFile)), qt.Qt.MatchExactly)
            # remove the previous empty cell icons
            if current_eb_item:
                self.set_row_bkg(current_eb_item[0].row(), qt.QColor("white"))
                current_eb_item[0].setIcon(qt.QIcon())
            # remove double reference
            fullfile = os.path.join(self.directory, file)
            self.emptyBeamFile = fullfile
            if fullfile == self.emptyCellFile:
                self.emptyCellFile = None
            elif fullfile == self.darkFile:
                self.darkFile = None
            elif fullfile == self.maskFile:
                self.maskFile = None
            elif fullfile in self.trashFiles:
                self.trashFiles.remove(fullfile)
            else:
                pass

    def _set_trash(self):
        # can be applied to the overall selection
        for current_item in self.selectedItems():
            if current_item is not None:
                row = current_item.row()
                ncol = self.columnCount()
                first_col_item = self.item(row, 0)
                file = first_col_item.text()
                self.set_row_bkg(row, qt.QColor("grey"))
                first_col_item.setIcon(getQIcon('cross.ico'))
                # remove double reference
                fullfile = os.path.join(self.directory, file)
                self.trashFiles.append(fullfile)
                if fullfile == self.emptyCellFile:
                    self.emptyCellFile = None
                elif fullfile == self.darkFile:
                    self.darkFile = None
                elif fullfile == self.emptyBeamFile:
                    self.emptyBeamFile = None
                elif fullfile == self.maskFile:
                    self.maskFile = None
                else:
                    pass

    def _set_mask(self):
        current_item = self.currentItem()
        if current_item is not None:
            row = current_item.row()
            ncol = self.columnCount()
            first_col_item = self.item(row, 0)
            file = first_col_item.text()
            first_col_item.setIcon(getQIcon('mask.ico'))
            # self.set_row_bkg(row, qt.QColor("cyan"))

            current_mask_item = self.findItems(os.path.basename(str(self.maskFile)), qt.Qt.MatchExactly)
            # remove the previous empty cell icons
            if current_mask_item:
                self.set_row_bkg(current_mask_item[0].row(), qt.QColor("white"))
                current_mask_item[0].setIcon(qt.QIcon())

            fullfile = os.path.join(self.directory, file)
            self.maskFile = fullfile
            # remove double reference
            if fullfile == self.darkFile:
                self.darkFile = None
            elif fullfile == self.emptyBeamFile:
                self.emptyBeamFile = None
            elif fullfile == self.emptyCellFile:
                self.maskFile = None
            elif fullfile in self.trashFiles:
                self.trashFiles.remove(fullfile)

    def get_sample_files(self):
        sampleList = []
        thicknessList = []
        for i in range(self.rowCount()):
            file = os.path.join(self.directory, self.item(i, 0).text())
            referenceList = [self.darkFile, self.emptyCellFile, self.emptyBeamFile, self.maskFile]
            referenceList += self.trashFiles
            if file not in referenceList:
                sampleList.append(file)
                try:
                    thicknessList.append(float(self.item(i, 2).text()))
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
        self.x0LineEdit.setValidator(qt.QDoubleValidator())
        self.y0LineEdit = qt.QLineEdit()
        self.y0LineEdit.setValidator(qt.QDoubleValidator())
        # sample to detector distance
        self.distanceLineEdit = qt.QLineEdit()
        self.distanceLineEdit.setValidator(qt.QDoubleValidator())
        # define the number of bins for azimutal averaging
        self.binsLineEdit = qt.QLineEdit('900')
        self.binsLineEdit.setValidator(qt.QIntValidator())
        # button to treat data
        self.treatButton = qt.QPushButton('treat now')
        self.treatButton.setIcon(getQIcon('gear.ico'))
        # parameter form layout
        formLayout = qt.QFormLayout()
        formLayout.addRow('x0 (pixels):', self.x0LineEdit)
        formLayout.addRow('y0 (pixels):', self.y0LineEdit)
        formLayout.addRow('distance (mm):', self.distanceLineEdit)
        formLayout.addRow('bins :', self.binsLineEdit)
        # general layout
        hlayout = qt.QHBoxLayout()
        # hlayout.addWidget(qt.QLabel('directory :'))
        # hlayout.addWidget(self.directoryLineEdit)
        # hlayout.addWidget(self.directoryPickerButton)
        # hlayout.addWidget(self.refreshButton)
        vlayout = qt.QVBoxLayout()
        vlayout.addLayout(formLayout)
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
            xeuss.set_beam_center(dark_file, x0=x0, y0=y0, new_entry=False)
            xeuss.azimutal_integration(dark_file, bins=nbins, mask=mask_file)
        # empty cell
        ec_file = self.table.emptyCellFile
        if ec_file is not None:
            nxlib.build_nexus_from_edf(ec_file)
            ec_file = ec_file.split('.')[0] + '.nxs'
            xeuss.set_beam_center(ec_file, x0=x0, y0=y0, new_entry=False)
            xeuss.azimutal_integration(ec_file, bins=nbins, mask=mask_file)
        # empty beam
        eb_file = self.table.emptyBeamFile
        if eb_file is not None:
            nxlib.build_nexus_from_edf(eb_file)
            eb_file = eb_file.split('.')[0] + '.nxs'
            xeuss.set_beam_center(eb_file, x0=x0, y0=y0, new_entry=False)
            xeuss.azimutal_integration(eb_file, bins=nbins, mask=mask_file)

        sampleList, thicknessList = self.table.get_sample_files()
        for sample, thickness in zip(sampleList, thicknessList):
            # nxfile = sample.split('.')[0]+'.nxs'
            try:
                nxlib.build_nexus_from_edf(sample)
                file = sample.split('.')[0]+'.nxs'
                xeuss.set_beam_center(file, x0=x0, y0=y0, new_entry=False)  # direct_beam_file=directbeam, new_entry=False)
                xeuss.azimutal_integration(file, bins=nbins, mask=mask_file)
                xeuss.resu(file, dark_file=dark_file, ec_file=ec_file, eb_file=eb_file,
                           distance=distance, thickness=thickness)
            except (KeyError, ValueError):
                print(('%s was ignored during the treatment') % (sample))
        self.table.refresh()

# TODO : save configuration in table
# TODO: tab for treated data inspection
# TODO : add checkbox when sample is treated
# TODO: pb with double legend


class FileSurvey(qt.QWidget):

    def __init__(self):
        super(FileSurvey, self).__init__()
        self.directoryLineEdit = qt.QLineEdit(parent=self)
        self.directoryPickerButton = qt.QPushButton()
        self.directoryPickerButton.setIcon(getQIcon('directory.ico'))
        self.refreshButton = qt.QPushButton()
        self.refreshButton.setIcon(getQIcon('refresh.ico'))

        self.tabWidget = qt.QTabWidget()
        self.edfTab = EdfTreatmentWidget()
        self.nxsTab = NexusTreatmentWidget()
        self.tabWidget.addTab(self.edfTab, 'raw data')
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
        self.setWindowTitle("Saxs Data Treatment")

        self.__asyncload = False
        central_wigdet = self.centralWidget()

        self.fileSurvey = FileSurvey()
        # self.nxsFileTable = NexusFileTable()

        # plot widget
        """Widget displaying information"""
        posInfo = [('X', lambda x, y: x),
                   ('Y', lambda x, y: y),
                   ('Data', self._zValue)]
        self.__plotWindow = PlotWindow(backend=None, resetzoom=True,
                                       autoScale=True, logScale=True,
                                       grid=True, curveStyle=True, colormap=True,
                                       aspectRatio=True, yInverted=True,
                                       copy=True, save=True, print_=True,
                                       control=True, position=posInfo,
                                       roi=False, mask=True, fit=True)
        legendWidget = self.__plotWindow.getLegendsDockWidget()
        # self.__plotWindow.setDockOptions()
        legendWidget.setAllowedAreas(qt.Qt.TopDockWidgetArea)
        self.profileTools = Profile.ProfileToolBar(parent=self.__plotWindow,
                                                   plot=self.__plotWindow)
        self.__plotWindow.addToolBar(self.profileTools)
        self.__plotWindow.setDefaultColormap(colors.Colormap(name='jet', normalization='log',
                                                             vmin=None, vmax=None, autoscaleMode='stddev3')
                                             )

        # treatment dock widget
        self.treatmentDock = qt.QDockWidget('treatment', self)
        self.editor = CommandTreatmentWidget(self)
        self.treatmentDock.setWidget(self.editor)
        self.treatmentDock.setFloating(False)
        self.addDockWidget(qt.Qt.RightDockWidgetArea, self.treatmentDock)
        self.treatmentDock.show()


        # directory picker layout
        spliter = qt.QSplitter(qt.Qt.Horizontal)
        spliter.addWidget(self.fileSurvey)
        spliter.addWidget(self.__plotWindow)
        spliter.setStretchFactor(1, 1)
        main_panel = qt.QWidget(self)
        layout = qt.QVBoxLayout()
        layout.addWidget(spliter)
        layout.setStretchFactor(spliter, 1)
        main_panel.setLayout(layout)
        self.setCentralWidget(main_panel)

        # connect signals
        self.fileSurvey.edfTab.edfSelectionChanged.connect(self.displayEdf)
        self.fileSurvey.nxsTab.nxsSelectionChanged.connect(self.displayNxs)
        self.editor.runClicked.connect(self.run_function)

    def run_function(self, cmdList):
        model = self.fileSurvey.nxsTab.treeWidget.treeview.findHdf5TreeModel()
        nrow = model.rowCount()
        files = []
        for n in range(nrow):
            index = model.index(n, 0, qt.QModelIndex())
            node = model.nodeFromIndex(index)
            filename = node.obj.filename
            model.removeH5pyObject(node.obj)
            node.obj.close()
            for script in cmdList:

                cmd = 'xeuss.' + script.replace('root', '\'' + filename + '\'')
                print(cmd)
                eval(cmd)
                # try:
                #     eval(cmd)
                # except:
                #     print('command'+cmd+'not performed on:' +filename)
            model.insertFile(filename, row=n)
            self.fileSurvey.nxsTab.treeWidget.treeview.expand(model.index(n, 0))
        self.fileSurvey.nxsTab.treeWidget.operationPerformed.emit()
        self.fileSurvey.nxsTab.tableWidget.on_selectionChanged()

    def displayEdf(self, file):
        """
        Called to update the dataviewer with the selected data.
        """
        self.__plotWindow.clear()
        self.__plotWindow.setKeepDataAspectRatio(True)
        self.__plotWindow.setXAxisLogarithmic(False)
        self.__plotWindow.setYAxisLogarithmic(False)
        if file.endswith('.edf'):
            dataObj = fabio.open(file)
            data = dataObj.data
            self.__plotWindow.clear()
            if 'Comments' in dataObj.header:
                legend = dataObj.header['Comments']
            else:
                legend = dataObj.filename
            self.__plotWindow.addImage(data, replace=True,
                                       legend=legend, xlabel='x [pixel]',
                                       ylabel='y [pixel]')

    def displayNxs(self, files):
        self.__plotWindow.clear()
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
                    self.__plotWindow.addCurve(nxd.axes[0], nxd.signal,
                                               yerror=nxd.errors,
                                               legend=legend,
                                               replace=False,
                                               color=col,
                                               xlabel=xlabel,
                                               ylabel=ylabel,
                                               resetzoom=True)
                    self.__plotWindow.setGraphXLabel(xlabel)
                    self.__plotWindow.setGraphYLabel(ylabel)
                    self.__plotWindow.setKeepDataAspectRatio(False)

    def _zValue(self, x, y):
        value = '-'
        valueZ = - float('inf')
        for image in self.__plotWindow.getAllImages():
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
        if os.path.isdir(directory):
            self.directory = directory
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


class NexusTreeWidget(qt.QWidget):
    operationPerformed = qt.pyqtSignal()

    def __init__(self):
        super(NexusTreeWidget, self).__init__()
        self.clear_last_btn = qt.QPushButton('clear last')
        self.clear_last_btn.setToolTip('Clear last performed calculation')
        self.clear_all_btn = qt.QPushButton('clear all')
        self.clear_all_btn.setToolTip('clear all performed calculations')
        self.sync_btn = qt.QPushButton('sync')
        self.sync_btn.setToolTip('synchronize the .nxs files')
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
        # layout
        hlayout = qt.QHBoxLayout()
        hlayout.addWidget(self.clear_last_btn)
        hlayout.addWidget(self.clear_all_btn)
        hlayout.addWidget(self.sync_btn)
        vlayout = qt.QVBoxLayout()
        vlayout.addLayout(hlayout)
        vlayout.addWidget(self.treeview)
        self.setLayout(vlayout)

        # connect signals
        self.sync_btn.clicked.connect(self.sync_all)
        self.clear_last_btn.clicked.connect(self.clear_last)
        self.clear_all_btn.clicked.connect(self.clear_all)

    def load_files(self, files):
        model = self.treeview.findHdf5TreeModel()
        model.clear()
        for file in files:
            model.insertFile(file, row=-1)
        # TODO : expand root
        # for row in range(model.rowCount()):
        #     # self.treeview.setExpanded(model.index(row, 0), True)
        #     index = model.index(row, 0)
        #     model.is


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
            self.treeview.expand(model.index(n, 0))
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
            self.treeview.expand(model.index(n, 0))
        self.operationPerformed.emit()

    def sync_all(self):
        model = self.treeview.findHdf5TreeModel()
        nrow = model.rowCount()
        for n in range(nrow):
            index = model.index(n, 0, qt.QModelIndex())
            model.synchronizeIndex(index)
            self.treeview.setExpanded(index, True)
        # self.treeview.expand(model.index(0, 0))
        # self.treeview.expandAll()
        self.operationPerformed.emit()


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
        spliter.setStretchFactor(0, 3)
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
        menu.addAction(concatAction)
        menu.addAction(convertAction)
        action = menu.exec_(self.mapToGlobal(event))

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
                xeuss.concat(firstFile, file=file, new_entry=True)
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

    def __init__(self, parent):
        super(CommandTreatmentWidget, self).__init__()
        self.run_btn = qt.QPushButton('run')
        self.runAll_btn = qt.QPushButton('run all')
        self.run_btn.clicked.connect(self.run)
        self.runAll_btn.clicked.connect(self.runAll)
        # self.remove_btn = qt.QPushButton('-')
        self.tabWidget = qt.QTabWidget()
        self.add_btn = qt.QPushButton('+')
        self.add_btn.setFixedSize(20,20)
        self.add_btn.clicked.connect(self.addTab)
        self.tabWidget.setCornerWidget(self.add_btn, corner=qt.Qt.TopLeftCorner)
        self.tabWidget.setTabsClosable(True)
        self.tabWidget.tabCloseRequested.connect(self.closeTabs)
        # widget = qt.QWidget(parent=self.tabWidget)
        # layoutTab = qt.QVBoxLayout()
        # layoutTab.addWidget(CodeEditor(self))
        # layoutTab.addStretch()
        # widget.setLayout(layoutTab)
        # self.tabWidget.addTab(widget, 'cmd1')
        self.tabWidget.addTab(CodeEditor(self), 'cmd1')
        # self.codeEditor = CodeEditor(self)
        # btn_layout = qt.QHBoxLayout()
        # btn_layout.addWidget(self.run_btn)
        # btn_layout.addWidget()
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
        self.tabWidget.addTab(CodeEditor(self), 'cmd'+str(count+1))
        # self.tabWidget.addTab(widget, 'cmd'+str(count+1))

    def run(self):
        widget = self.tabWidget.currentWidget()
        text = widget.text()
        self.runClicked.emit([text])

    def runAll(self):
        count = self.tabWidget.count()
        l = []
        for i in range(count):
            widget = self.tabWidget.widget(i)
            l.append(widget.text())
        self.runClicked.emit(l)


class CodeEditor(qt.QLineEdit):
    """
    QLineEdit widget with treatment function autocompletion
    """
    def __init__(self, parent=None):
        super().__init__()
        # self.setTabStopDistance(
            # qt.QFontMetricsF(self.font()).horizontalAdvance(' ') * 4)
        # self.highlighter = PythonHighlighter(self.document())
        completer = qt.QCompleter(COMPLETER_NAMES)
        self.setCompleter(completer)


def main():
    # unlock hdf5 files for file access during plotting
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    app = qt.QApplication([])
    # sys.excepthook = qt.exceptionHandler
    # warnings.filterwarnings("ignore", category=mplDeprecation)
    window = SaxsUtily()
    window.show()
    result = app.exec_()
    # remove ending warnings relative to QTimer
    app.deleteLater()
    sys.exit(result)


def testTree():
    app = qt.QApplication([])
    window = NexusTreatmentWidget()
    window.show()
    file = '/home/achennev/Bureau/PIL pour tiago/PIL NP/2021-07-21_TOC/2021-07-21_TOC_0_89050.nxs'
    file1 = '/home/achennev/Bureau/PIL pour tiago/PIL NP/2021-07-21_TOC/2021-07-21_TOC_0_89051.nxs'
    window.set_directory('/home/achennev/Bureau/PIL pour tiago/PIL NP/2021-07-21_TOC')
    result = app.exec_()
    # remove ending warnings relative to QTimer
    app.deleteLater()
    sys.exit(result)


if __name__ == "__main__":
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    app = qt.QApplication([])
    # warnings.filterwarnings("ignore", category=mplDeprecation)
    window = SaxsUtily()
    window.show()
    folder = '/home/achennev/Bureau/PIL pour tiago/PIL NP/2021-07-21_TOC'
    window.fileSurvey.directoryLineEdit.setText(folder)
    result = app.exec_()
    app.deleteLater()
    sys.exit(result)

