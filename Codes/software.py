import enaml
from basic import *
from enaml.qt.qt_application import QtApplication


qm = QuineMcCluskey([])

if __name__ == '__main__':
    with enaml.imports():
        from gui import functionView
    app = QtApplication()
    view = functionView(function=qm)
    view.show()
    app.start()