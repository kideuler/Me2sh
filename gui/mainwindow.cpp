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
    : QMainWindow(parent), drawGeoArea(new DrawGeoArea(this))
{ 
    screen = QGuiApplication::primaryScreen();
    QRect screenGeometry = screen->geometry();
    screenWidth = screenGeometry.width()*0.9;
    screenHeight = screenGeometry.height()*0.9;

    tabWidget = new QTabWidget(this);
    tabWidget->addTab(drawGeoArea, tr("Geometry"));
    tabWidget->addTab(new QWidget(), tr("Mesh"));
    tabWidget->addTab(new QWidget(), tr("Simulation"));

    msgBox = new ConsoleOutput(this);

    QSettings settings;
    int height = settings.value("Split", 600).toInt();
    QSplitter *centralWidget = new QSplitter(Qt::Vertical, this);
    centralWidget->addWidget(tabWidget);
    centralWidget->addWidget(msgBox);
    centralWidget->setSizes(QList<int>{int(0.8*screenHeight), int(0.2*screenHeight)});

    setCentralWidget(centralWidget);
    createActions();
    createMenus();
    createSidebars();

    setWindowTitle(tr("Me2sh"));
    resize(1600, 800);

    connect(tabWidget, &QTabWidget::currentChanged, this, &MainWindow::updateSidebar);

    drawGeoArea->clearImage();
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
    QPushButton *drawCircleButton = new QPushButton(tr("Draw Circle"), this);
    geometryLayout->addWidget(drawCircleButton);
    connect(drawCircleButton, &QPushButton::clicked, drawGeoArea, &DrawGeoArea::drawCircle);

    QPushButton *drawEllipseButton = new QPushButton(tr("Draw Ellipse"), this);
    geometryLayout->addWidget(drawEllipseButton);
    connect(drawEllipseButton, &QPushButton::clicked, drawGeoArea, &DrawGeoArea::drawEllipse);

    QPushButton *drawSplineButton = new QPushButton(tr("Draw Spline"), this);
    geometryLayout->addWidget(drawSplineButton);
    connect(drawSplineButton, &QPushButton::clicked, drawGeoArea, &DrawGeoArea::drawSpline);

    QPushButton *drawBSplineButton = new QPushButton(tr("Draw BSpline"), this);
    geometryLayout->addWidget(drawBSplineButton);
    connect(drawBSplineButton, &QPushButton::clicked, drawGeoArea, &DrawGeoArea::drawBSpline);

    #ifdef USE_GEO
    QPushButton *drawBezierButton = new QPushButton(tr("Draw Bezier"), this);
    geometryLayout->addWidget(drawBezierButton);
    connect(drawBezierButton, &QPushButton::clicked, drawGeoArea, &DrawGeoArea::drawBezier);
    #endif

    // Add a stretch at the end to push everything to the top
    geometryLayout->addStretch(1);




    meshSidebar->setLayout(new QVBoxLayout);
    meshSidebar->layout()->addWidget(new QLabel("Mesh Sidebar"));






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
    msgBox->clear();
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
