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
#include "ConsoleOutput.hpp"

class DrawGeoArea;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    // Constructor
    MainWindow(QWidget *parent = nullptr);

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

    // Slot for clear screen
    void clearScreen();

    
private:
    void createActions();
    void createMenus();

    void createSidebars();
    QTabWidget *tabWidget;

    QStackedWidget *sidebarStack;
    QWidget *geometrySidebar;
    QWidget *meshSidebar;
    QWidget *simulationSidebar;

    DrawGeoArea *drawGeoArea;
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