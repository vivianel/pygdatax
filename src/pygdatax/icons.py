import pygdatax.resources
from silx.gui import qt
import os


def getQIcon(name):
    ans = None
    path = os.path.abspath(os.path.dirname(pygdatax.resources.__file__))
    fullname = os.path.join(path, name)
    icon = qt.QIcon(fullname)
    return icon
