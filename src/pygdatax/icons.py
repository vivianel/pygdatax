from importlib.resources import path
import resources
from silx.gui import qt

def getQIcon(name):
    ans = None
    with path(resources, name) as a:
        ans = str(a.absolute())
    return qt.QIcon(ans)