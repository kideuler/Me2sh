#include "drawgeometry.hpp"

#include <QMouseEvent>
#include <QPainter>
#include <QApplication>

#if defined(QT_PRINTSUPPORT_LIB)
#include <QtPrintSupport/qtprintsupportglobal.h>
#if QT_CONFIG(printdialog)
#include <QPrinter>
#include <QPrintDialog>
#endif
#endif

DrawGeoArea::DrawGeoArea(QWidget *parent)
    : QWidget(parent)
{
    setAttribute(Qt::WA_StaticContents);
}

void DrawGeoArea::setPenColor(const QColor &newColor)
{
    myPenColor = newColor;
}

void DrawGeoArea::setPenWidth(int newWidth)
{
    myPenWidth = newWidth;
}

void DrawGeoArea::clearImage()
{
    image.fill(qRgb(255, 255, 255));
    modified = true;
    
    update();
}

void DrawGeoArea::mousePressEvent(QMouseEvent *event)
{
    if (Qt::ControlModifier == QApplication::keyboardModifiers()){
        if (event->button() == Qt::LeftButton) {
            lastPoint = event->position().toPoint();
            firstPoint = event->position().toPoint();
            firstPointIndex = 0;
            scribbling = true;
            firstpointinit = true;
        }
    } else if (Qt::ShiftModifier == QApplication::keyboardModifiers()){
        if (event->button() == Qt::LeftButton) {
            double x = (double)event->position().x()/(double)width();
            double y = 1.0-(double)event->position().y()/(double)height();
            drawPoint(event->position().toPoint());
            update();
        }
    } else{
        if (event->button() == Qt::LeftButton) {
            lastPoint = event->position().toPoint();
            firstPoint = event->position().toPoint();
            scribbling = true;
            firstpointinit = true;
        }
    }
}

void DrawGeoArea::mouseMoveEvent(QMouseEvent *event)
{
    if (Qt::ControlModifier == QApplication::keyboardModifiers()){
        if ((event->buttons() & Qt::LeftButton) && scribbling){
            QPoint p = event->position().toPoint();
            
            if (!(p.x()==lastPoint.x() || p.y() == lastPoint.y())) {
                double x = (double)lastPoint.x()/(double)width();
                double y = 1.0-(double)lastPoint.y()/(double)height();

                drawLineTo(p);
            }
        }
    } else {
        if ((event->buttons() & Qt::LeftButton) && scribbling){
            QPoint p = event->position().toPoint();
            
            if (!(p.x()==lastPoint.x() || p.y() == lastPoint.y())) {
                double x = (double)lastPoint.x()/(double)width();
                double y = 1.0-(double)lastPoint.y()/(double)height();

                if (firstpointinit){
                    firstpointinit = false;
                }else {
                }

                drawLineTo(p);
            }
        }
    }
}

void DrawGeoArea::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton && scribbling) {
        drawLineTo(event->position().toPoint());
        drawLineTo(firstPoint);
        scribbling = false;
        
        
    }
}

void DrawGeoArea::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    QRect dirtyRect = event->rect();
    painter.drawImage(dirtyRect, image, dirtyRect);
}

void DrawGeoArea::resizeEvent(QResizeEvent *event)
{
    if (width() > image.width() || height() > image.height()) {
        int newWidth = qMax(width() + 128, image.width());
        int newHeight = qMax(height() + 128, image.height());
        resizeImage(&image, QSize(newWidth, newHeight));
        update();
    }
    QWidget::resizeEvent(event);
}

void DrawGeoArea::drawLineTo(const QPoint &endPoint)
{
    QPainter painter(&image);
    painter.setPen(QPen(myPenColor, myPenWidth, Qt::SolidLine, Qt::RoundCap,
                        Qt::RoundJoin));
    painter.drawLine(lastPoint, endPoint);
    modified = true;

    int rad = (myPenWidth / 2) + 2;
    update(QRect(lastPoint, endPoint).normalized()
                                     .adjusted(-rad, -rad, +rad, +rad));
    lastPoint = endPoint;
}

void DrawGeoArea::drawPoint(const QPoint &P){
    QPainter painter(&image);
    painter.setBrush(Qt::black);
    painter.drawEllipse(P, myPenWidth, myPenWidth);
    modified = true;
    update();
}

void DrawGeoArea::resizeImage(QImage *image, const QSize &newSize)
{
    if (image->size() == newSize)
        return;

    QImage newImage(newSize, QImage::Format_RGB32);
    newImage.fill(qRgb(255, 255, 255));
    QPainter painter(&newImage);
    painter.drawImage(QPoint(0, 0), *image);
    *image = newImage;
}