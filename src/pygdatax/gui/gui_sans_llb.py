import sys
import os
import warnings

from matplotlib import mplDeprecation
import numpy as np
import fabio
from silx.gui import qt, colors
import silx.io as sio
from silx.io.utils import is_group, is_dataset, is_file
from silx.io.nxdata import is_NXroot_with_default_NXdata, get_default
from numpy.random import randint

from pygdatax.icons import getQIcon
from pygdatax.instruments import sansllb
import pygdatax.gui as gui






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
        self.nxsTab = gui.NexusTreatmentWidget()

        # layout
        hlayout = qt.QHBoxLayout()
        hlayout.addWidget(qt.QLabel('directory :'))
        hlayout.addWidget(self.directoryLineEdit)
        hlayout.addWidget(self.directoryPickerButton)
        hlayout.addWidget(self.refreshButton)
        vlayout = qt.QVBoxLayout()
        vlayout.addLayout(hlayout)
        vlayout.addWidget(self.nxsTab)
        self.setLayout(vlayout)

        # connect signals
        self.directoryLineEdit.textChanged.connect(self.set_directory)
        self.directoryPickerButton.clicked.connect(self.choose_directory)
        self.refreshButton.clicked.connect(self.nxsTab.tableWidget.refresh)

    def set_directory(self):
        text = self.directoryLineEdit.text()
        self.nxsTab.set_directory(text)
        os.chdir(text)

    def choose_directory(self):
        basedir = os.path.expanduser("~")
        fname = qt.QFileDialog.getExistingDirectory(self, 'Select data directory', basedir,
                                                    options=qt.QFileDialog.DontUseNativeDialog)
        if fname:
            self.directoryLineEdit.setText(fname)
            # self.edfTab.table.setDirectory(fname)
            # self.nxsTab.setDirectory(fname)


class SansLlbMainWindow(qt.QMainWindow):
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
        self.plotWindow = gui.DataView()
        # treatment dock widget
        self.treatmentDock = qt.QDockWidget('command', self)
        # self.treatmentDock.setStyleSheet("border: 5px solid black")
        self.treatmentDock.setFeatures(qt.QDockWidget.DockWidgetFloatable |
                                       qt.QDockWidget.DockWidgetMovable)
        self.editor = gui.CommandTreatmentWidget(self, module=sansllb)
        self.treatmentDock.setWidget(self.editor)
        self.treatmentDock.setFloating(False)
        # replace the addTabbedwidget metho of the plot window
        self.plotWindow._dockWidgets.append(self.treatmentDock)
        self.plotWindow.addDockWidget(qt.Qt.RightDockWidgetArea, self.treatmentDock)
        # self.treatmentDock.show()
        # self.treatmentDock.setAllowedAreas(qt.Qt.BottomDockWidgetArea)
        # self.addDockWidget(qt.Qt.RightDockWidgetArea, self.treatmentDock)
        # Add function list widget
        self.functionListWidget = gui.FunctionListWidget(parent=self, module=sansllb)
        self.functionListDock = qt.QDockWidget('functions', self)
        self.functionListDock.setFeatures(qt.QDockWidget.DockWidgetFloatable |
                                          qt.QDockWidget.DockWidgetMovable)
        self.functionListDock.setFloating(False)
        self.functionListDock.setWidget(self.functionListWidget)
        self.plotWindow._dockWidgets.append(self.functionListDock)
        self.plotWindow.tabifyDockWidget(self.treatmentDock, self.functionListDock)
        # self.plotWindow.addDockWidget(qt.Qt.RightDockWidgetArea, self.functionListDock)
        self.functionListDock.show()
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
        # nxs table dispplay
        self.fileSurvey.nxsTab.nxsSelectionChanged.connect(self.displayNxs)
        # treatement run
        self.editor.runClicked.connect(self.run_function)
        self.fileSurvey.nxsTab.treeWidget.selectedNodeChanged.connect(self.displayNxTree)
        # function list ran
        self.functionListWidget.runFunction.connect(self.run_function)

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
        self.worker = gui.FunctionWorker(cmdList, selectedFiles)
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
        if type(cmdList) is not list:
            cmdList = [cmdList]
        selectedFiles = self.fileSurvey.nxsTab.tableWidget.get_selectedFiles()
        model = self.fileSurvey.nxsTab.treeWidget.treeview.findHdf5TreeModel()
        model.clear()
        for file in selectedFiles:
            for script in cmdList:
                for line in script.splitlines():
                    cmd = 'sansllb.' + line.replace('root', '\'' + file.replace('\\', '/') + '\'')
                    try:
                        eval(cmd)
                        # print(cmd)
                    except:
                        print('command : '+cmd+'not performed on:' + file)
            # model.insertFile(file)
        self.fileSurvey.nxsTab.treeWidget.operationPerformed.emit()
        self.fileSurvey.nxsTab.tableWidget.on_selectionChanged()
        self.fileSurvey.refreshButton.clicked.emit()

    def displayNxs(self, files):
        self.plotWindow.clear()
        c = ['blue', 'red', 'green', 'black', 'yellow', 'grey', 'magenta', 'cyan',
             'darkGreen', 'darkBrown', 'darkCyan', 'darkYellow', 'darkMagenta']

        for i, file in enumerate(files):
            h5file = sio.open(file)
            if is_NXroot_with_default_NXdata(h5file):
                nxd = get_default(h5file)
                legend = os.path.basename(file).split('.')[0] + '/'
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


def main():
    # unlock hdf5 files for file access during plotting
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    warnings.filterwarnings("ignore", category=mplDeprecation)
    app = qt.QApplication([])
    # sys.excepthook = qt.exceptionHandler
    window = SansLlbMainWindow()
    window.show()
    result = app.exec_()
    # remove ending warnings relative to QTimer
    app.deleteLater()
    sys.exit(result)


if __name__ == "__main__":
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    from pygdatax.instruments import sansllb
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
    window = SansLlbMainWindow()
    window.show()
    folder = '/home/achennev/python/pa20_psi'
    window.fileSurvey.directoryLineEdit.setText(folder)
    ####################################################################
    # test 3 detector view
    #####################################################################
    ############################################################################
    # w.set_function_decription(moduledescription.FunctionDescription(sansllb.make_reduction_package))
    # w.show()
    # w.get_commandLine()

    # splash.finish(window)
    result = app.exec_()
    app.deleteLater()
    sys.exit(result)

