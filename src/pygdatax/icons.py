import pygdatax.resources
from silx.gui import qt
import os


def getQIcon(name):
    ans = None
    path = os.path.abspath(os.path.dirname(pygdatax.resources.__file__))
    ans = os.path.join(path, name)
    print(ans)
    return qt.QIcon(ans)
