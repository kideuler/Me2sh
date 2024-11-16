#ifndef DRAWMESH_H
#define DRAWMESH_H

#include <QColor>
#include <QImage>
#include <QPoint>
#include <QWidget>

#include "geometry.hpp"

class DrawMeshArea : public QWidget
{
    Q_OBJECT

public:
    DrawMeshArea(Me2sh_Geometry *geometry, QWidget *parent = nullptr);
    void clearImage();

    bool isModified() const { return modified; }

    Me2sh_Geometry *geo = new Me2sh_Geometry();

    void drawgeometry();

    void chanegeMeshSize();
    void changeElementType();

public slots:

protected:
    void paintEvent(QPaintEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;
private:
    void resizeImage(QImage *image, const QSize &newSize);
    void drawShape();
    void clearTemporaryLayer();
    void drawAxis(QPainter &painter); // Add this line
    void drawAxisLabels(QPainter &painter); // Add this line

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

    double h_target = 0.05;
    int GmshMeshAlgorithm = 6;
    int GmshRecombinationAlgorithm = 3;
    bool ElementType = 0;
};

#endif