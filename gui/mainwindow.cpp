#include "mainwindow.hpp"
#include "drawmesh.hpp"

#include <QApplication>
#include <QColorDialog>
#include <QFileDialog>
#include <QImageWriter>
#include <QInputDialog>
#include <QMenuBar>
#include <QMessageBox>
#include <QCloseEvent>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), drawMeshArea(new DrawMeshArea(this))
{   
    drawMeshArea->coords.cols = 2;
    drawMeshArea->segments.cols = 2;
    setCentralWidget(drawMeshArea);

    createActions();
    createMenus();

    setWindowTitle(tr("Me2sh"));
    resize(1000, 800);
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    event->accept();
}

void MainWindow::penColor()
{
    QColor newColor = QColorDialog::getColor(drawMeshArea->penColor());
    if (newColor.isValid())
        drawMeshArea->setPenColor(newColor);
}

void MainWindow::penWidth()
{
    bool ok;
    int newWidth = QInputDialog::getInt(this, tr("DrawMesh"),
                                        tr("Select pen width:"),
                                        drawMeshArea->penWidth(),
                                        1, 50, 1, &ok);
    if (ok)
        drawMeshArea->setPenWidth(newWidth);
}

void MainWindow::about()
{
    QMessageBox::about(this, tr("About Me2sh"),
            tr("<p>The <b>DrawMesh</b> example shows how to use QMainWindow as the "
               "base widget for an application, and how to reimplement some of "
               "QWidget's event handlers to receive the events generated for "
               "the application's widgets:</p><p> We reimplement the mouse event "
               "handlers to facilitate drawing, the paint event handler to "
               "update the application and the resize event handler to optimize "
               "the application's appearance. In addition we reimplement the "
               "close event handler to intercept the close events before "
               "terminating the application.</p><p> The example also demonstrates "
               "how to use QPainter to draw an image in real time, as well as "
               "to repaint widgets.</p>"));
}

void MainWindow::sethtarget(){
    bool ok;
    double h_target = QInputDialog::getDouble(this, tr("DrawMesh"),
                                        tr("Select target h value:"),
                                        drawMeshArea->h_target,
                                        0, 1, 4, &ok);
    if (ok)
        drawMeshArea->h_target = h_target;
}

void MainWindow::setsmoothingiters(){
    bool ok;
    int iters = QInputDialog::getInt(this, tr("DrawMesh"),
                                        tr("Choose Maximum number of smoothing iterations:"),
                                        drawMeshArea->max_smoothing_iters,
                                        1, 1000, 1, &ok);
    if (ok)
        drawMeshArea->max_smoothing_iters = iters;
}

void MainWindow::triangulate(){
    drawMeshArea->mesh.Triangulate(drawMeshArea->coords);
    drawMeshArea->showMesh();
}

void MainWindow::constrainedTriangulate(){
    drawMeshArea->mesh.Triangulate(drawMeshArea->coords, drawMeshArea->segments);
    drawMeshArea->showMesh();
}

void MainWindow::refineMesh(){
    drawMeshArea->mesh.Refine(drawMeshArea->h_target);
    drawMeshArea->showMesh();
}

void MainWindow::smoothMesh(){
    drawMeshArea->mesh.Smooth(drawMeshArea->max_smoothing_iters);
    drawMeshArea->showMesh();
}

void MainWindow::computeVolumeLengthMetric(){
    drawMeshArea->mesh.Compute_volume_length_metric();
    drawMeshArea->showQuality();
}

void MainWindow::createActions()
{
    penColorAct = new QAction(tr("&Pen Color..."), this);
    connect(penColorAct, &QAction::triggered, this, &MainWindow::penColor);

    penWidthAct = new QAction(tr("Pen &Width..."), this);
    connect(penWidthAct, &QAction::triggered, this, &MainWindow::penWidth);

    clearScreenAct = new QAction(tr("&Clear Screen"), this);
    clearScreenAct->setShortcut(tr("Ctrl+L"));
    connect(clearScreenAct, &QAction::triggered,
            drawMeshArea, &DrawMeshArea::clearImage);

    aboutAct = new QAction(tr("&About"), this);
    connect(aboutAct, &QAction::triggered, this, &MainWindow::about);

    aboutQtAct = new QAction(tr("About &Qt"), this);
    connect(aboutQtAct, &QAction::triggered, qApp, &QApplication::aboutQt);

    // mesh actions
    sethtargetAct = new QAction(tr("&Set h target"), this);
    connect(sethtargetAct, &QAction::triggered, this, &MainWindow::sethtarget);

    setsmoothingitersAct = new QAction(tr("&Set Max Smoothing Iterations"), this);
    connect(setsmoothingitersAct, &QAction::triggered, this, &MainWindow::setsmoothingiters);

    triangulateAct = new QAction(tr("&Triangulate"), this);
    triangulateAct->setShortcut(tr("Ctrl+T"));
    connect(triangulateAct, &QAction::triggered, this, &MainWindow::triangulate);

    constrainedTriangulateAct = new QAction(tr("&Constrained Triangulation"), this);
    constrainedTriangulateAct->setShortcut(tr("Ctrl+Shift+T"));
    connect(constrainedTriangulateAct, &QAction::triggered, this, &MainWindow::constrainedTriangulate);

    refineMeshAct = new QAction(tr("&Refine Mesh"), this);
    refineMeshAct->setShortcut(tr("Ctrl+R"));
    connect(refineMeshAct, &QAction::triggered, this, &MainWindow::refineMesh);

    smoothMeshAct = new QAction(tr("&Smooth Mesh"), this);
    smoothMeshAct->setShortcut(tr("Ctrl+S"));
    connect(smoothMeshAct, &QAction::triggered, this, &MainWindow::smoothMesh);

    ComputeVolumeLengthMetricAct = new QAction(tr("&Compute Volume Length Metric"), this);
    ComputeVolumeLengthMetricAct->setShortcut(tr("Ctrl+A"));
    connect(ComputeVolumeLengthMetricAct, &QAction::triggered, this, &MainWindow::computeVolumeLengthMetric);
}

void MainWindow::createMenus()
{
    optionMenu = new QMenu(tr("&Options"), this);
    optionMenu->addAction(penColorAct);
    optionMenu->addAction(penWidthAct);
    optionMenu->addSeparator();
    optionMenu->addAction(clearScreenAct);

    helpMenu = new QMenu(tr("&Help"), this);
    helpMenu->addAction(aboutAct);
    helpMenu->addAction(aboutQtAct);

    meshMenu = new QMenu(tr("&Mesh"), this);
    meshMenu->addAction(sethtargetAct);
    meshMenu->addAction(setsmoothingitersAct);
    meshMenu->addSeparator();
    meshMenu->addAction(triangulateAct);
    meshMenu->addAction(constrainedTriangulateAct);
    meshMenu->addSeparator();
    meshMenu->addAction(refineMeshAct);
    meshMenu->addAction(smoothMeshAct);
    meshMenu->addSeparator();
    meshMenu->addAction(ComputeVolumeLengthMetricAct);

    

    menuBar()->addMenu(optionMenu);
    menuBar()->addMenu(helpMenu);
    menuBar()->addMenu(meshMenu);
}
