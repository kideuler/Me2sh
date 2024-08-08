#ifndef DRAWMESH_H
#define DRAWMESH_H

#include <QColor>
#include <QImage>
#include <QPoint>
#include <QWidget>

#include "TriMesh.hpp"
#include "Spline.hpp"

class DrawMeshArea : public QWidget
{
    Q_OBJECT

public:
    DrawMeshArea(QWidget *parent = nullptr);

    void setPenColor(const QColor &newColor);
    void setPenWidth(int newWidth);

    // drawing mesh
    void showMesh();
    void showQuality();

    bool isModified() const { return modified; }
    QColor penColor() const { return myPenColor; }
    int penWidth() const { return myPenWidth; }


    // trimesh arrays
    Matrix<double> coords;
    Matrix<int> segments;
    Mesh2D mesh;
    double h_target = 0.05;
    int max_smoothing_iters = 100;

    // splines
    Spline2D spline;

public slots:
    void clearImage();
    void print();

protected:
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void paintEvent(QPaintEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;

private:
    void drawLineTo(const QPoint &endPoint);
    void resizeImage(QImage *image, const QSize &newSize);

    bool modified = false;
    bool scribbling = false;
    int myPenWidth = 2;
    QColor myPenColor = Qt::blue;
    QImage image;
    QPoint lastPoint;
    QPoint firstPoint;
    int firstPointIndex = -1;
    bool firstpointinit = false;
};

#endif