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
#include <QTabWidget>
#include <QStackedWidget>
#include <QVBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <chrono>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{ 
    geo = std::make_shared<Me2sh_Geometry>();
    mesh = std::make_shared<Me2sh_Mesh>();
    PythonTerminal = new ConsoleOutput(this);
    drawGeoArea = new DrawGeoArea(PythonTerminal, geo, this);
    drawMeshArea = new DrawMeshArea(PythonTerminal, geo, mesh, this);

    screen = QGuiApplication::primaryScreen();
    QRect screenGeometry = screen->geometry();
    screenWidth = screenGeometry.width()*0.9;
    screenHeight = screenGeometry.height()*0.9;

    tabWidget = new QTabWidget(this);
    tabWidget->addTab(drawGeoArea, tr("Geometry"));
    tabWidget->addTab(drawMeshArea, tr("Mesh"));
    tabWidget->addTab(new QWidget(), tr("Simulation"));

    QSettings settings;
    int height = settings.value("Split", 600).toInt();
    QSplitter *centralWidget = new QSplitter(Qt::Vertical, this);
    centralWidget->addWidget(tabWidget);
    centralWidget->addWidget(PythonTerminal);
    centralWidget->setSizes(QList<int>{int(0.8*screenHeight), int(0.2*screenHeight)});

    setCentralWidget(centralWidget);
    createActions();
    createMenus();
    createSidebars();

    setWindowTitle(QString("Me2sh - Version %1").arg(PROJECT_VERSION));
    resize(1600, 800);

    connect(tabWidget, &QTabWidget::currentChanged, this, &MainWindow::updateSidebar);
    connect(tabWidget, &QTabWidget::currentChanged, this, &MainWindow::handleTabChange);

    drawGeoArea->clearImage();
    
}

void MainWindow::handleTabChange(int index)
{
    if (tabWidget->widget(index) == drawMeshArea) {
        // if (!drawMeshArea->hasMesh){
        //     drawMeshArea->drawgeometry();
        // } else {
        //     drawMeshArea->drawgeometry();
        //     drawMeshArea->displayMesh();
        // }
    }
}

void MainWindow::createSidebars()
{
    sidebarStack = new QStackedWidget(this);

    geometrySidebar = new QWidget(this);
    meshSidebar = new QWidget(this);
    simulationSidebar = new QWidget(this);

    // Add widgets to sidebars as needed
    QVBoxLayout *geometryLayout = new QVBoxLayout;
    geometrySidebar->setLayout(geometryLayout);

    // Add the label at the top
    QLabel *geometryLabel = new QLabel("Geometry Sidebar", this);
    geometryLayout->addWidget(geometryLabel, 0, Qt::AlignTop);

    QFrame *line = new QFrame(this);
    line->setFrameShape(QFrame::HLine);
    line->setFrameShadow(QFrame::Sunken);
    geometryLayout->addWidget(line);

    // Add label for geometry creation
    QLabel *geometryCreationLabel = new QLabel("Create", this);
    geometryLayout->addWidget(geometryCreationLabel);

    // Add the push buttons right under the label
    QPushButton *drawCircleButton = new QPushButton(tr("Create Circle"), this);
    geometryLayout->addWidget(drawCircleButton);
    connect(drawCircleButton, &QPushButton::clicked, drawGeoArea, &DrawGeoArea::drawCircle);

    QPushButton *drawEllipseButton = new QPushButton(tr("Create Ellipse"), this);
    geometryLayout->addWidget(drawEllipseButton);
    connect(drawEllipseButton, &QPushButton::clicked, drawGeoArea, &DrawGeoArea::drawEllipse);

    QPushButton *drawRectangleButton = new QPushButton(tr("Create Rectangle"), this);
    geometryLayout->addWidget(drawRectangleButton);
    connect(drawRectangleButton, &QPushButton::clicked, drawGeoArea, &DrawGeoArea::drawRectangle);

    QPushButton *drawSplineButton = new QPushButton(tr("Create Spline"), this);
    geometryLayout->addWidget(drawSplineButton);
    connect(drawSplineButton, &QPushButton::clicked, drawGeoArea, &DrawGeoArea::drawSpline);

    QPushButton *drawBSplineButton = new QPushButton(tr("Create BSpline"), this);
    geometryLayout->addWidget(drawBSplineButton);
    connect(drawBSplineButton, &QPushButton::clicked, drawGeoArea, &DrawGeoArea::drawBSpline);

    #ifdef USE_GEO
    QPushButton *drawBezierButton = new QPushButton(tr("Draw Bezier"), this);
    geometryLayout->addWidget(drawBezierButton);
    connect(drawBezierButton, &QPushButton::clicked, drawGeoArea, &DrawGeoArea::drawBezier);
    #endif

    QLabel *geometryModificationLabel = new QLabel("Modify", this);
    geometryLayout->addWidget(geometryModificationLabel);

    QPushButton *fuseButton = new QPushButton(tr("Fuse Overlapping"), this);
    geometryLayout->addWidget(fuseButton);
    connect(fuseButton, &QPushButton::clicked, drawGeoArea, &DrawGeoArea::FuseAll);

    QPushButton *makeExteriorGeometryButton = new QPushButton(tr("Make Exterior Geometry"), this);
    geometryLayout->addWidget(makeExteriorGeometryButton);
    connect(makeExteriorGeometryButton, &QPushButton::clicked, drawGeoArea, &DrawGeoArea::MakeExteriorGeometry);

    // Add a stretch at the end to push everything to the top
    geometryLayout->addStretch(1);




    // Add widgets to sidebars as needed
    QVBoxLayout *meshLayout = new QVBoxLayout;
    meshSidebar->setLayout(meshLayout);

    // Add the label at the top
    QLabel *meshLabel = new QLabel("Mesh Sidebar", this);
    meshLayout->addWidget(meshLabel, 0, Qt::AlignTop);


    QFrame *line_mesh = new QFrame(this);
    line_mesh->setFrameShape(QFrame::HLine);
    line_mesh->setFrameShadow(QFrame::Sunken);
    meshLayout->addWidget(line_mesh);

    // Add label for mesh options
    QLabel *MeshParamLabel = new QLabel("Mesh Parameters", this);
    meshLayout->addWidget(MeshParamLabel);

    // Add the push buttons right under the label
    QPushButton *changeMeshSizeButton = new QPushButton(tr("Change Mesh Size"), this);
    meshLayout->addWidget(changeMeshSizeButton);
    connect(changeMeshSizeButton, &QPushButton::clicked, drawMeshArea, &DrawMeshArea::chanegeMeshSize);

    QPushButton *changeElementTypeButton = new QPushButton(tr("Change Element Type"), this);
    meshLayout->addWidget(changeElementTypeButton);
    connect(changeElementTypeButton, &QPushButton::clicked, drawMeshArea, &DrawMeshArea::changeElementType);

    QPushButton *changeMeshAlgorithmButton = new QPushButton(tr("Change Mesh Algorithm"), this);
    meshLayout->addWidget(changeMeshAlgorithmButton);
    connect(changeMeshAlgorithmButton, &QPushButton::clicked, drawMeshArea, &DrawMeshArea::changeMeshAlgorithm);

    QFrame *line_mesh2 = new QFrame(this);
    line_mesh->setFrameShape(QFrame::HLine);
    line_mesh->setFrameShadow(QFrame::Sunken);
    meshLayout->addWidget(line_mesh2);

    // Add label for mesh options
    QLabel *MeshToolLabel = new QLabel("Meshing Tools", this);
    meshLayout->addWidget(MeshToolLabel);

    QPushButton *generateMeshButton = new QPushButton(tr("Generate Mesh"), this);
    meshLayout->addWidget(generateMeshButton);
    connect(generateMeshButton, &QPushButton::clicked, drawMeshArea, &DrawMeshArea::generateMesh);

    meshLayout->addStretch(1);






    simulationSidebar->setLayout(new QVBoxLayout);
    simulationSidebar->layout()->addWidget(new QLabel("Simulation Sidebar"));





    sidebarStack->addWidget(geometrySidebar);
    sidebarStack->addWidget(meshSidebar);
    sidebarStack->addWidget(simulationSidebar);

    // Update layout to include sidebarStack
    QSplitter *mainSplitter = new QSplitter(Qt::Horizontal, this);
    mainSplitter->addWidget(sidebarStack);
    mainSplitter->addWidget(centralWidget());
    mainSplitter->setSizes(QList<int>{int(0.1*screenWidth), int(0.9*screenWidth)});
    setCentralWidget(mainSplitter);
}

void MainWindow::updateSidebar(int index)
{
    sidebarStack->setCurrentIndex(index);
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    event->accept();
}

void MainWindow::penColor()
{
    QColor newColor = QColorDialog::getColor(drawGeoArea->penColor());
    if (newColor.isValid())
        drawGeoArea->setPenColor(newColor);
}

void MainWindow::penWidth()
{
    bool ok;
    int newWidth = QInputDialog::getInt(this, tr("me2sh"),
                                        tr("Select pen width:"),
                                        drawGeoArea->penWidth(),
                                        1, 50, 1, &ok);
    if (ok)
        drawGeoArea->setPenWidth(newWidth);
}

void MainWindow::about()
{
    QMessageBox::about(this, tr("About Me2sh"),
            tr("Write Something better Later"));
}

void MainWindow::clearScreen(){
    drawGeoArea->clearImage();
    drawMeshArea->clearImage();
    PythonTerminal->clear();
}

void MainWindow::createActions()
{
    connect(tabWidget, &QTabWidget::currentChanged, this, &MainWindow::updateSidebar);

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


    menuBar()->addMenu(optionMenu);
    menuBar()->addMenu(helpMenu);
}
