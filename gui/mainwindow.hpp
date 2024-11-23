#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QList>
#include <QMainWindow>
#include <QTabWidget>
#include <QApplication>
#include <QColorDialog>
#include <QFileDialog>
#include <QImageWriter>
#include <QInputDialog>
#include <QMenuBar>
#include <QMessageBox>
#include <QCloseEvent>
#include <QSplitter>
#include <QSettings>
#include <QTabWidget>
#include <QStackedWidget>
#include <QVBoxLayout>
#include <QLabel>
#include <QScreen>
#include <QGuiApplication>
#include <chrono>

#include "drawgeometry.hpp"
#include "drawmesh.hpp"
#include "ConsoleOutput.hpp"

class DrawGeoArea;
class DrawMeshArea;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    // Constructor
    MainWindow(QWidget *parent = nullptr);

    // clear screen
    void clearScreen();

protected:
    // Event handler for close event
    void closeEvent(QCloseEvent *event) override;

private slots:
    // Slot for tab widget
    void updateSidebar(int index);

    // Slot for pen color
    void penColor();

    // Slot for pen width
    void penWidth();

    // Slot for about
    void about();

    // Handle tab change
    void handleTabChange(int index);

    
private:
    void createActions();
    void createMenus();

    void createSidebars();
    QTabWidget *tabWidget;

    std::shared_ptr<Me2sh_Geometry> geo;
    std::shared_ptr<Me2sh_Mesh> mesh;

    QStackedWidget *sidebarStack;
    QWidget *geometrySidebar;
    QWidget *meshSidebar;
    QWidget *simulationSidebar;

    DrawGeoArea *drawGeoArea;
    DrawMeshArea *drawMeshArea;
    ConsoleOutput *msgBox;
    

    QMenu *optionMenu;
    QMenu *helpMenu;

    QAction *exitAct;
    QAction *penColorAct;
    QAction *penWidthAct;
    QAction *clearScreenAct;
    QAction *aboutAct;
    QAction *aboutQtAct;
    QAction *clearAct;

    QScreen *screen;

    int screenWidth;
    int screenHeight;
};

#endif