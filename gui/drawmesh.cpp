#include "drawmesh.hpp"

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

DrawMeshArea::DrawMeshArea(QWidget *parent)
    : QWidget(parent)
{
    setAttribute(Qt::WA_StaticContents);
}

void DrawMeshArea::setPenColor(const QColor &newColor)
{
    myPenColor = newColor;
}

void DrawMeshArea::setPenWidth(int newWidth)
{
    myPenWidth = newWidth;
}

void DrawMeshArea::clearImage()
{
    image.fill(qRgb(255, 255, 255));
    modified = true;
    mesh.reset();
    coords = Matrix<double>(0,0);
    coords.cols = 2;
    segments = Matrix<int>(0,0);
    segments.cols = 2;
    Cpoints = Matrix<double>(0,0);
    Cpoints.cols = 2;
    spline.reset();
    update();
}

void DrawMeshArea::mousePressEvent(QMouseEvent *event)
{
    if (Qt::ShiftModifier == QApplication::keyboardModifiers()){
        if (event->button() == Qt::LeftButton) {
            double x = (double)event->position().x()/(double)width();
            double y = 1.0-(double)event->position().y()/(double)height();
            Cpoints.arr.push_back(x);
            Cpoints.arr.push_back(y);
            Cpoints.rows++;
            drawPoint(event->position().toPoint());
            update();
        }
    } else{
        if (event->button() == Qt::LeftButton) {
            lastPoint = event->position().toPoint();
            firstPoint = event->position().toPoint();
            firstPointIndex = coords.nrows();
            scribbling = true;
            firstpointinit = true;
        }
    }
}

void DrawMeshArea::mouseMoveEvent(QMouseEvent *event)
{
    if ((event->buttons() & Qt::LeftButton) && scribbling){
        QPoint p = event->position().toPoint();
        
        if (!(p.x()==lastPoint.x() || p.y() == lastPoint.y())) {
            double x = (double)lastPoint.x()/(double)width();
            double y = 1.0-(double)lastPoint.y()/(double)height();
            coords.arr.push_back(x);
            coords.arr.push_back(y);
            coords.rows++;

            if (firstpointinit){
                firstpointinit = false;
            }else {
                segments.arr.push_back(coords.nrows()-2);
                segments.arr.push_back(coords.nrows()-1);
                segments.rows++;
            }

            drawLineTo(p);
        }
    }
}

void DrawMeshArea::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton && scribbling) {
        drawLineTo(event->position().toPoint());
        drawLineTo(firstPoint);
        segments.arr.push_back(coords.nrows()-1);
        segments.arr.push_back(firstPointIndex);
        segments.rows++;
        scribbling = false;
        
        
    }
}

void DrawMeshArea::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    QRect dirtyRect = event->rect();
    painter.drawImage(dirtyRect, image, dirtyRect);
}

void DrawMeshArea::resizeEvent(QResizeEvent *event)
{
    if (width() > image.width() || height() > image.height()) {
        int newWidth = qMax(width() + 128, image.width());
        int newHeight = qMax(height() + 128, image.height());
        resizeImage(&image, QSize(newWidth, newHeight));
        update();
    }
    QWidget::resizeEvent(event);
}

void DrawMeshArea::drawLineTo(const QPoint &endPoint)
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

void DrawMeshArea::drawPoint(const QPoint &P){
    QPainter painter(&image);
    painter.setBrush(Qt::black);
    painter.drawEllipse(P, myPenWidth, myPenWidth);
    modified = true;
    update();
}

void DrawMeshArea::resizeImage(QImage *image, const QSize &newSize)
{
    if (image->size() == newSize)
        return;

    QImage newImage(newSize, QImage::Format_RGB32);
    newImage.fill(qRgb(255, 255, 255));
    QPainter painter(&newImage);
    painter.drawImage(QPoint(0, 0), *image);
    *image = newImage;
}

void DrawMeshArea::print()
{
#if defined(QT_PRINTSUPPORT_LIB) && QT_CONFIG(printdialog)
    QPrinter printer(QPrinter::HighResolution);

    QPrintDialog printDialog(&printer, this);
    if (printDialog.exec() == QDialog::Accepted) {
        QPainter painter(&printer);
        QRect rect = painter.viewport();
        QSize size = image.size();
        size.scale(rect.size(), Qt::KeepAspectRatio);
        painter.setViewport(rect.x(), rect.y(), size.width(), size.height());
        painter.setWindow(image.rect());
        painter.drawImage(0, 0, image);
    }
#endif // QT_CONFIG(printdialog)
}


void DrawMeshArea::showMesh(){
    image.fill(qRgb(255, 255, 255));
    QPainter painter(&image);
    painter.setPen(QPen(Qt::black, 1, Qt::SolidLine));

    for (int i = 0; i<mesh.elems.nrows(); i++){
        int i0 = mesh.elems(i,0);
        int i1 = mesh.elems(i,1);
        int i2 = mesh.elems(i,2);

        int x0 = mesh.coords(i0,0)*width();
        int y0 = (1.0-mesh.coords(i0,1))*height();
        int x1 = mesh.coords(i1,0)*width();
        int y1 = (1.0-mesh.coords(i1,1))*height();
        int x2 = mesh.coords(i2,0)*width();
        int y2 = (1.0-mesh.coords(i2,1))*height();

        painter.drawLine(x0,y0,x1,y1);
        painter.drawLine(x1,y1,x2,y2);
        painter.drawLine(x2,y2,x0,y0);
    }
    update();
}

void DrawMeshArea::showQuality(){
    image.fill(qRgb(255, 255, 255));
    QPainter painter(&image);
    painter.setBrush(QBrush(Qt::blue));

    for (int i = 0; i<mesh.elems.nrows(); i++){
        int i0 = mesh.elems(i,0);
        int i1 = mesh.elems(i,1);
        int i2 = mesh.elems(i,2);

        int x0 = mesh.coords(i0,0)*width();
        int y0 = (1.0-mesh.coords(i0,1))*height();
        int x1 = mesh.coords(i1,0)*width();
        int y1 = (1.0-mesh.coords(i1,1))*height();
        int x2 = mesh.coords(i2,0)*width();
        int y2 = (1.0-mesh.coords(i2,1))*height();
        QPoint points[3] = {QPoint(x0,y0),QPoint(x1,y1),QPoint(x2,y2)};

        double q = mesh.quality[i];
        painter.setBrush(QColor(255*(1-q),0,255*q));

        painter.drawPolygon(points, 3);
    }
    update();
}

void DrawMeshArea::showSpline(){

    int npoints = 500;
    double dt = 1.0/(double)(npoints-1);

    std::array<double,2> xy = spline.eval(0.0);
    int x = xy[0]*width();
    int y = (1.0-xy[1])*height();
    lastPoint = QPoint(x,y);

    for (double n = 1; n<npoints; ++n){
        double t = n*dt;
        std::array<double,2> xy = spline.eval(t);
        x = xy[0]*width();
        y = (1.0-xy[1])*height();
        QPoint p = QPoint(x,y);
        drawLineTo(p);
        lastPoint = p;
    }

    npoints = 50;
    double dh = 0.01;
    dt = 1.0/(double)(npoints);
    int nx,ny;
    for (double n = 1; n<npoints; ++n){
        double t = n*dt;
        std::array<double,2> xy = spline.eval(t);
        x = xy[0]*width();
        y = (1.0-xy[1])*height();
        lastPoint = QPoint(x,y);

        std::array<double,2> N = spline.normal(t);
        N[0] = xy[0] + dh*N[0];
        N[1] = xy[1] + dh*N[1];
        nx = N[0]*width();
        ny = (1.0-N[1])*height();
        QPoint np = QPoint(nx,ny);
        drawLineTo(np);
    }

    update();
}