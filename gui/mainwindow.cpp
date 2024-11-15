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
#include <chrono>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), drawGeoArea(new DrawGeoArea(this))
{   
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
    centralWidget->setSizes(QList<int>{height, 800-height});

    setCentralWidget(centralWidget);
    createActions();
    createMenus();
    createSidebars();

    setWindowTitle(tr("Me2sh"));
    resize(1600, 800);

    connect(tabWidget, &QTabWidget::currentChanged, this, &MainWindow::updateSidebar);
}

void MainWindow::createSidebars()
{
    sidebarStack = new QStackedWidget(this);

    geometrySidebar = new QWidget(this);
    meshSidebar = new QWidget(this);
    simulationSidebar = new QWidget(this);

    // Add widgets to sidebars as needed
    geometrySidebar->setLayout(new QVBoxLayout);
    geometrySidebar->layout()->addWidget(new QLabel("Geometry Sidebar"));

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
    mainSplitter->setSizes(QList<int>{100, width()-100});
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
