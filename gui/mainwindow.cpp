#include "mainwindow.hpp"

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
#include <chrono>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), drawMeshArea(new DrawMeshArea(this))
{   
    drawMeshArea->coords.cols = 2;
    drawMeshArea->segments.cols = 2;
    drawMeshArea->Cpoints.cols = 2;

    setCentralWidget(drawMeshArea);
    msgBox = new ConsoleOutput(this);

    QSettings settings;
    int height = settings.value("Split", 600).toInt();
    QSplitter *centralWidget = new QSplitter(Qt::Vertical, this);
    centralWidget->addWidget(drawMeshArea);
    centralWidget->addWidget(msgBox);
    centralWidget->setSizes(QList<int>{height, 800-height});

    setCentralWidget(centralWidget);
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

void MainWindow::clearScreen(){
    drawMeshArea->clearImage();
    msgBox->clear();
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
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    if (drawMeshArea->coords.nrows() > 0){
        drawMeshArea->mesh.Triangulate(drawMeshArea->coords);
        drawMeshArea->showMesh();
    } else if (drawMeshArea->spline.npoints() > 0){
        drawMeshArea->spline.create_segments(drawMeshArea->h_target, drawMeshArea->coords, drawMeshArea->segments);
        drawMeshArea->mesh.Triangulate(drawMeshArea->coords);
        drawMeshArea->showMesh();
    } else if (drawMeshArea->bezier.npoints() > 0){
        drawMeshArea->bezier.create_segments(drawMeshArea->h_target, drawMeshArea->coords, drawMeshArea->segments);
        drawMeshArea->mesh.Triangulate(drawMeshArea->coords);
        drawMeshArea->showMesh();
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    msgBox->addMessage(QString::fromStdString("Triangulation time: " + std::to_string(elapsed_seconds.count()) + " seconds"));
    
}

void MainWindow::constrainedTriangulate(){
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    if (drawMeshArea->coords.nrows() > 0){
        drawMeshArea->mesh.Triangulate(drawMeshArea->coords, drawMeshArea->segments);
        drawMeshArea->showMesh();
    } else if (drawMeshArea->spline.npoints() > 0){
        drawMeshArea->spline.create_segments(drawMeshArea->h_target, drawMeshArea->coords, drawMeshArea->segments);
        drawMeshArea->mesh.Triangulate(drawMeshArea->coords, drawMeshArea->segments);
        drawMeshArea->showMesh();
    } else if (drawMeshArea->bezier.npoints() > 0){
        drawMeshArea->bezier.create_segments(drawMeshArea->h_target, drawMeshArea->coords, drawMeshArea->segments);
        drawMeshArea->mesh.Triangulate(drawMeshArea->coords, drawMeshArea->segments);
        drawMeshArea->showMesh();
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    msgBox->addMessage(QString::fromStdString("Constrained Triangulation time: " + std::to_string(elapsed_seconds.count()) + " seconds"));
}

void MainWindow::refineMesh(){
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    drawMeshArea->mesh.Refine(drawMeshArea->h_target);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    msgBox->addMessage(QString::fromStdString("Mesh Refinement time: " + std::to_string(elapsed_seconds.count()) + " seconds"));
    drawMeshArea->showMesh();
}

void MainWindow::smoothMesh(){
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    drawMeshArea->mesh.Smooth(drawMeshArea->max_smoothing_iters);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    msgBox->addMessage(QString::fromStdString("Mesh Smoothing time: " + std::to_string(elapsed_seconds.count()) + " seconds"));
    drawMeshArea->showMesh();
}

void MainWindow::computeVolumeLengthMetric(){
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    drawMeshArea->mesh.Compute_volume_length_metric();
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    msgBox->addMessage(QString::fromStdString("Compute Mesh Quality: " + std::to_string(elapsed_seconds.count()) + " seconds"));
    drawMeshArea->showQuality();
}

void MainWindow::makeSpline(){
    drawMeshArea->bezier.reset();
    drawMeshArea->spline.reset();
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    drawMeshArea->spline.init(drawMeshArea->Cpoints);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    msgBox->addMessage(QString::fromStdString("Computed Spline: " + std::to_string(elapsed_seconds.count()) + " seconds"));
    drawMeshArea->showSpline();
}

void MainWindow::makeBezier(){
    drawMeshArea->bezier.reset();
    drawMeshArea->spline.reset();
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    drawMeshArea->bezier.init(drawMeshArea->Cpoints);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    msgBox->addMessage(QString::fromStdString("Computed Bezier: " + std::to_string(elapsed_seconds.count()) + " seconds"));
    drawMeshArea->showBezier();
}

void MainWindow::solvePoisson(){
    drawMeshArea->fem.reset();
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    drawMeshArea->mesh.compute_boundary_nodes();
    drawMeshArea->fem.init(drawMeshArea->mesh, 0.0, 1.0, 1.0);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    msgBox->addMessage(QString::fromStdString("Initialized Laplace Equation: " + std::to_string(elapsed_seconds.count()) + " seconds"));


    start = std::chrono::system_clock::now();
    drawMeshArea->fem.assemble();
    drawMeshArea->fem.apply_dbc();
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    msgBox->addMessage(QString::fromStdString("Assembled Laplace Equation: " + std::to_string(elapsed_seconds.count()) + " seconds"));

    start = std::chrono::system_clock::now();
    drawMeshArea->fem.solve();
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    msgBox->addMessage(QString::fromStdString("Solved Laplace Equation: " + std::to_string(elapsed_seconds.count()) + " seconds"));

    drawMeshArea->showPoisson();
}

void MainWindow::solveEikonal(){

    bool ok;
    double alpha = QInputDialog::getDouble(this, tr("DrawMesh"),
                                        tr("Select alpha value for Eikonal Equation\n|\\grad(u)|^2= 1, using \nalpha^2 \\laplacian(v)-v = 0 \nu = -alpha*log(v):"),
                                        0.001,
                                        0, 1, 5, &ok);
    if (ok){
        drawMeshArea->fem_eikonal.reset();
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        drawMeshArea->mesh.compute_boundary_nodes();
        drawMeshArea->fem_eikonal.init(drawMeshArea->mesh, alpha);
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        msgBox->addMessage(QString::fromStdString("Initialized Eikonal Equation: " + std::to_string(elapsed_seconds.count()) + " seconds"));

        start = std::chrono::system_clock::now();
        drawMeshArea->fem_eikonal.assemble();
        drawMeshArea->fem_eikonal.apply_dbc();
        end = std::chrono::system_clock::now();
        elapsed_seconds = end-start;
        msgBox->addMessage(QString::fromStdString("Assembled Eikonal Equation: " + std::to_string(elapsed_seconds.count()) + " seconds"));

        start = std::chrono::system_clock::now();
        drawMeshArea->fem_eikonal.solve();
        end = std::chrono::system_clock::now();
        elapsed_seconds = end-start;
        msgBox->addMessage(QString::fromStdString("Solved Eikonal Equation: " + std::to_string(elapsed_seconds.count()) + " seconds"));

        drawMeshArea->showEikonal();
    }
}

void MainWindow::createActions()
{
    penColorAct = new QAction(tr("&Pen Color..."), this);
    connect(penColorAct, &QAction::triggered, this, &MainWindow::penColor);

    penWidthAct = new QAction(tr("Pen &Width..."), this);
    connect(penWidthAct, &QAction::triggered, this, &MainWindow::penWidth);

    clearAct = new QAction(tr("&Clear Screen"), this);
    clearAct->setShortcut(tr("Ctrl+L"));
    connect(clearAct, &QAction::triggered,
            this, &MainWindow::clearScreen);

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

    // spline actions
    makeSplineAct = new QAction(tr("&Make Cubic Spline"), this);
    makeSplineAct->setShortcut(tr("Ctrl+1"));
    connect(makeSplineAct, &QAction::triggered, this, &MainWindow::makeSpline);

    makeBezierAct = new QAction(tr("&Make Bezier Spline"), this);
    makeBezierAct->setShortcut(tr("Ctrl+2"));
    connect(makeBezierAct, &QAction::triggered, this, &MainWindow::makeBezier);

    // simulation actions
    solvePoissonAct = new QAction(tr("&Solve Laplace Equation"), this);
    solvePoissonAct->setShortcut(tr("Ctrl+Shift+1"));
    connect(solvePoissonAct, &QAction::triggered, this, &MainWindow::solvePoisson);

    solveEikonalAct = new QAction(tr("&Solve Eikonal Equation"), this);
    solveEikonalAct->setShortcut(tr("Ctrl+Shift+2"));
    connect(solveEikonalAct, &QAction::triggered, this, &MainWindow::solveEikonal);
}

void MainWindow::createMenus()
{
    // options menu
    optionMenu = new QMenu(tr("&Options"), this);
    optionMenu->addAction(penColorAct);
    optionMenu->addAction(penWidthAct);
    optionMenu->addSeparator();
    optionMenu->addAction(clearAct);

    // help menu
    helpMenu = new QMenu(tr("&Help"), this);
    helpMenu->addAction(aboutAct);
    helpMenu->addAction(aboutQtAct);

    // mesh menu
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

    // spline menu
    splineMenu = new QMenu(tr("&Spline"), this);
    splineMenu->addAction(makeSplineAct);
    splineMenu->addAction(makeBezierAct);

    // simulation menu
    simMenu = new QMenu(tr("&Simulation"), this);
    simMenu->addAction(solvePoissonAct);
    simMenu->addAction(solveEikonalAct);
    

    menuBar()->addMenu(optionMenu);
    menuBar()->addMenu(helpMenu);
    menuBar()->addMenu(meshMenu);
    menuBar()->addMenu(splineMenu);
    menuBar()->addMenu(simMenu);
}
