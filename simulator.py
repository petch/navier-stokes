from PyQt4.QtGui import *
from PyQt4.QtCore import *
import sys
import os

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.splitter = QSplitter(Qt.Horizontal)
        self.splitter.addWidget(self.sideArea())
        self.splitter.addWidget(self.workArea())
        self.setCentralWidget(self.splitter)
        
        self.settings = QSettings('MMR', 'Simulator')
        self.restoreGeometry(self.settings.value('geometry', ''))
        self.restoreState(self.settings.value('state', ''))
        self.splitter.restoreState(self.settings.value('splitter', ''))
        self.setWindowTitle('Simulator')
        self.show()

    def closeEvent(self, event):
        self.settings.setValue('geometry', self.saveGeometry())
        self.settings.setValue('state', self.saveState())
        self.settings.setValue('splitter', self.splitter.saveState())
        QMainWindow.closeEvent(self, event)

    def sideArea(self):
        side = QWidget()
        box = QVBoxLayout()

        menu = QWidget()
        menuBox = QHBoxLayout()
        menuBox.addWidget(QLabel('Projects'), 1)
        menuBox.addWidget(QPushButton('New'))
        menuBox.addWidget(QPushButton('Clone'))
        menuBox.addWidget(QPushButton('Delete'))
        menu.setLayout(menuBox)
        box.addWidget(menu)

        tree = QTreeWidget()
        tree.setHeaderHidden(True)
        
        root, projects, files = next(os.walk('Projects'))
        for project in sorted(projects):
            item = QTreeWidgetItem([project])
            subroot, dirs, subfiles = next(os.walk(root + '/' + project))
            for d in sorted(dirs):
                child = QTreeWidgetItem([d])
                item.addChild(child)
            tree.addTopLevelItem(item)
        box.addWidget(tree)

        side.setLayout(box)
        return side


    def workArea(self):
        area = QWidget()
        vBox = QVBoxLayout()

        menu = QWidget()
        vBox.addWidget(menu)
        menuBox = QHBoxLayout()
        menu.setLayout(menuBox)
        menuBox.addWidget(QLabel('Project1'), 1)
        menuBox.addWidget(QPushButton('Run'))
        menuBox.addWidget(QPushButton('Test'))

        view = QFrame()
        view.setFrameShape(QFrame.StyledPanel)
        vBox.addWidget(view, 1)

        area.setLayout(vBox)
        return area

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec_())
 
if __name__ == '__main__':
    main()