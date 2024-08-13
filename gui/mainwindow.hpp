#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QList>
#include <QMainWindow>

#include "drawmesh.hpp"
#include "ConsoleOutput.hpp"

class DrawMeshArea;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);

protected:
    void closeEvent(QCloseEvent *event) override;

private slots:
    void penColor();
    void penWidth();
    void about();
    void clearScreen();

    // meshing slots
    void sethtarget();
    void setsmoothingiters();
    void triangulate();
    void constrainedTriangulate();
    void refineMesh();
    void smoothMesh();
    void computeVolumeLengthMetric();

    // splining slots
    void makeSpline();
    void makeBezier();

    // simulation slots
    void solvePoisson();
    void solveEikonal();

private:
    void createActions();
    void createMenus();

    DrawMeshArea *drawMeshArea;
    ConsoleOutput *msgBox;

    QMenu *optionMenu;
    QMenu *helpMenu;
    QMenu *meshMenu;
    QMenu *splineMenu;
    QMenu *simMenu;

    QAction *exitAct;
    QAction *penColorAct;
    QAction *penWidthAct;
    QAction *printAct;
    QAction *clearScreenAct;
    QAction *aboutAct;
    QAction *aboutQtAct;
    QAction *clearAct;

    // meshing actions
    QAction *sethtargetAct;
    QAction *setsmoothingitersAct;
    QAction *triangulateAct;
    QAction *constrainedTriangulateAct;
    QAction *refineMeshAct;
    QAction *smoothMeshAct;
    QAction *ComputeVolumeLengthMetricAct;

    // splining actions
    QAction *makeSplineAct;
    QAction *makeBezierAct;

    // simulation actions
    QAction *solvePoissonAct;
    QAction *solveEikonalAct;
};

#endif