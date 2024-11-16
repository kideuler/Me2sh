#ifndef DRAWGEOMETRY_H
#define DRAWGEOMETRY_H

#include <QColor>
#include <QImage>
#include <QPoint>
#include <QWidget>

#include "geometry.hpp"

class DrawGeoArea : public QWidget
{
    Q_OBJECT

public:
    DrawGeoArea(Me2sh_Geometry *geometry, QWidget *parent = nullptr);
    void clearImage();

    void setPenColor(const QColor &newColor);
    void setPenWidth(int newWidth);

    void drawSpline();
    void drawBSpline();
    void drawCircle();
    void drawEllipse();
    void drawRectangle();

    #ifdef USE_GEO
    void drawBezier();
    #endif

    bool isModified() const { return modified; }
    QColor penColor() const { return myPenColor; }
    int penWidth() const { return myPenWidth; }

    Me2sh_Geometry *geo = new Me2sh_Geometry();

    double h_target = 0.02;
    int max_smoothing_iters = 100;
    bool selectingCenter = false;
    bool selectingRect = false;
    double rx = 0.0;
    double ry = 0.0;

public slots:

protected:
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void paintEvent(QPaintEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;
    bool event(QEvent *event) override;
private:
    void drawPoint(const QPoint &P, QImage *img = nullptr);
    void resizeImage(QImage *image, const QSize &newSize);
    void drawTemporaryCircle(const QPoint &center);
    void drawTemporarRectangle(const QPoint &center);
    void drawShape();
    void clearTemporaryLayer();
    void drawAxis(QPainter &painter); // Add this line
    void drawAxisLabels(QPainter &painter); // Add this line

    bool modified = false;
    bool scribbling = false;
    int myPenWidth = 1;
    QColor myPenColor = Qt::blue;
    QImage image;
    void drawLineTo(const QPoint &endPoint, QImage *img = nullptr, QColor PenColor = nullptr);
    QImage tempImage;
    QPoint lastPoint;
    QPoint firstPoint;
    int firstPointIndex = -1;
    bool firstpointinit = false;
};

#endif