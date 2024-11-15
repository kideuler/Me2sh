#ifndef DRAWGEOMETRY_H
#define DRAWGEOMETRY_H

#include <QColor>
#include <QImage>
#include <QPoint>
#include <QWidget>

class DrawGeoArea : public QWidget
{
    Q_OBJECT

public:
    DrawGeoArea(QWidget *parent = nullptr);
    void clearImage();

    void setPenColor(const QColor &newColor);
    void setPenWidth(int newWidth);

    bool isModified() const { return modified; }
    QColor penColor() const { return myPenColor; }
    int penWidth() const { return myPenWidth; }

    double h_target = 0.02;
    int max_smoothing_iters = 100;

public slots:

protected:
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void paintEvent(QPaintEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;

private:
    void drawLineTo(const QPoint &endPoint);
    void drawPoint(const QPoint &P);
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