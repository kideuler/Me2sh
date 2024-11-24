#ifndef DRAWMESH_H
#define DRAWMESH_H

#include "ConsoleOutput.hpp"

#include <QColor>
#include <QImage>
#include <QTimer>
#include <QPoint>
#include <QWidget>
#include <memory>

#include "geometry.hpp"
#include "meshing.hpp"

class DrawMeshArea : public QWidget
{
    Q_OBJECT

public:
    DrawMeshArea(ConsoleOutput *PyTerm, std::shared_ptr<Me2sh_Geometry> geometry, std::shared_ptr<Me2sh_Mesh> mesh, QWidget *parent = nullptr);
    void clearImage();

    void startAnimation();

    bool isModified() const { return modified; }

    std::shared_ptr<Me2sh_Geometry> geo;
    std::shared_ptr<Me2sh_Mesh> mesh;
    ConsoleOutput *PythonTerminal;

    bool hasMesh = false;

    void drawgeometry();

    void chanegeMeshSize();
    void changeElementType();
    void changeMeshAlgorithm();

    void generateMesh();

    void displayMesh();
    void ShowMeshColorGradient();

    double h_target = 0.01;
    int GmshMeshAlgorithm = 6;
    int TriMeshAlgorithm = 0;
    int QuadMeshAlgorithm = 0;
    int GmshRecombinationAlgorithm = 3;
    int ElementType = 1;

public slots:

protected:
    void paintEvent(QPaintEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;

private slots:
    void updateAnimation();

private:
    void resizeImage(QImage *image, const QSize &newSize);
    void drawShape();
    void clearTemporaryLayer();
    void drawAxis(QPainter &painter); 
    void drawAxisLabels(QPainter &painter); 
    void drawPoint(const QPoint &P, QImage *img = nullptr);
    void drawPolygon(const std::vector<QPoint> &points, QImage *img = nullptr);

    bool modified = false;
    QImage image;
    void drawLineTo(const QPoint &endPoint, QImage *img = nullptr, QColor PenColor = nullptr);
    QImage tempImage;
    QPoint lastPoint;
    QPoint firstPoint;
    int firstPointIndex = -1;
    bool firstpointinit = false;
    QColor myPenColor = Qt::black; 
    int myPenWidth = 1;

    QTimer *animationTimer;
    double animationPhase = 0.0;
};

#endif