#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QList>
#include <QMainWindow>

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

    // meshing slots
    void sethtarget();
    void setsmoothingiters();
    void triangulate();
    void constrainedTriangulate();
    void refineMesh();
    void smoothMesh();
    void computeVolumeLengthMetric();

private:
    void createActions();
    void createMenus();

    DrawMeshArea *drawMeshArea;

    QMenu *optionMenu;
    QMenu *helpMenu;
    QMenu *meshMenu;

    QAction *exitAct;
    QAction *penColorAct;
    QAction *penWidthAct;
    QAction *printAct;
    QAction *clearScreenAct;
    QAction *aboutAct;
    QAction *aboutQtAct;

    // meshing actions
    QAction *sethtargetAct;
    QAction *setsmoothingitersAct;
    QAction *triangulateAct;
    QAction *constrainedTriangulateAct;
    QAction *refineMeshAct;
    QAction *smoothMeshAct;
    QAction *ComputeVolumeLengthMetricAct;
};

#endif