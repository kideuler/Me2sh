#include "drawgeometry.hpp"

#include <QMouseEvent>
#include <QPainter>
#include <QApplication>
#include <QInputDialog>
#include <QHoverEvent>

#if defined(QT_PRINTSUPPORT_LIB)
#include <QtPrintSupport/qtprintsupportglobal.h>
#if QT_CONFIG(printdialog)
#include <QPrinter>
#include <QPrintDialog>
#endif
#endif

DrawGeoArea::DrawGeoArea(std::shared_ptr<Me2sh_Geometry> geometry, QWidget *parent)
    : QWidget(parent), geo(geometry)
{
    setAttribute(Qt::WA_StaticContents);
    setMouseTracking(true);
    setAttribute(Qt::WA_Hover);
    // Initialize images with the correct size
    image = QImage(size(), QImage::Format_RGB32);
    image.fill(qRgb(255, 255, 255)); // Ensure the background is white

    tempImage = QImage(size(), QImage::Format_ARGB32_Premultiplied);
    tempImage.fill(Qt::transparent); // Ensure the temporary image is transparent

    update();
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
    tempImage.fill(Qt::transparent);
    image.fill(qRgb(255, 255, 255));
    clearTemporaryLayer();
    modified = true;
    exterior_geom = false;
    geo->clear();
    update();
}

void DrawGeoArea::clearTemporaryLayer()
{
    tempImage.fill(Qt::transparent); // Clear only the temporary image
    geo->cleartemp();
    update();
}

void DrawGeoArea::mousePressEvent(QMouseEvent *event)
{
    if (selectingCenter) {
        if (event->button() == Qt::LeftButton) {
            double scale = (double)std::max(width(), height());
            double centerX = (double)event->position().x() / scale;
            double centerY = 1.0 - (double)event->position().y() / scale;
            geo->addEllipse(centerX, centerY, rx, ry);
            drawShape();
            selectingCenter = false;
            clearTemporaryLayer();
        }
    } else if (selectingRect) {
        if (event->button() == Qt::LeftButton) {
            double scale = (double)std::max(width(), height());
            double centerX = (double)event->position().x() / scale;
            double centerY = 1.0 - (double)event->position().y() / scale;
            geo->addRectangle(centerX, centerY, rx, ry);
            drawShape();
            selectingRect = false;
            clearTemporaryLayer();
        }
    } else if (Qt::ControlModifier == QApplication::keyboardModifiers()){
        if (event->button() == Qt::LeftButton) {
            lastPoint = event->position().toPoint();
            firstPoint = event->position().toPoint();
            firstPointIndex = 0;
            scribbling = true;
            firstpointinit = true;
        }
    } else if (Qt::ShiftModifier == QApplication::keyboardModifiers()){
        if (event->button() == Qt::LeftButton) {
            double scale = (double)std::max(width(), height());
            double x = (double)event->position().x()/scale;
            double y = 1.0-(double)event->position().y()/scale;
            drawPoint(event->position().toPoint(), &tempImage);
            geo->points.push_back({x, y, 0.0});
            update();
        }
    } else{
        if (event->button() == Qt::LeftButton) {
            double scale = (double)std::max(width(), height());
            lastPoint = event->position().toPoint();
            firstPoint = event->position().toPoint();
            scribbling = true;
            firstpointinit = true;
        }
    }
}

void DrawGeoArea::drawTemporaryCircle(const QPoint &center)
{
    int scale = std::max(width(), height());
    int rxi = (int)(rx*scale);
    int ryi = (int)(ry*scale);
    QPainter painter(&tempImage);
    painter.setPen(QPen(Qt::red, myPenWidth, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin));
    painter.drawEllipse(center, rxi, ryi);
    modified = true;
    update();
}

void DrawGeoArea::drawTemporarRectangle(const QPoint &center)
{
    int scale = std::max(width(), height());
    int rxi = (int)(rx*scale);
    int ryi = (int)(ry*scale);
    QPainter painter(&tempImage);
    painter.setPen(QPen(Qt::red, myPenWidth, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin));
    painter.drawRect(center.x()-rxi, center.y()-ryi, 2*rxi, 2*ryi);
    modified = true;
    update();
}

void DrawGeoArea::mouseMoveEvent(QMouseEvent *event)
{
    if (Qt::ControlModifier == QApplication::keyboardModifiers()){
        if ((event->buttons() & Qt::LeftButton) && scribbling){
            QPoint p = event->position().toPoint();
            
            if (!(p.x()==lastPoint.x() || p.y() == lastPoint.y())) {
                double scale = (double)std::max(width(), height());
                double x = (double)lastPoint.x()/scale;
                double y = 1.0-(double)lastPoint.y()/scale;
                drawPoint(event->position().toPoint());
                geo->points.push_back({x, y, 0.0});
                drawLineTo(p);
            }
        }
    } else {
        if ((event->buttons() & Qt::LeftButton) && scribbling){
            QPoint p = event->position().toPoint();
            
            if (!(p.x()==lastPoint.x() || p.y() == lastPoint.y())) {
                double scale = (double)std::max(width(), height());
                double x = (double)lastPoint.x()/scale;
                double y = 1.0-(double)lastPoint.y()/scale;
                geo->points.push_back({x, y, 0.0});
                drawLineTo(p, &tempImage);

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
        drawLineTo(firstPoint, &tempImage);
        scribbling = false;
        
        
    }
}

bool DrawGeoArea::event(QEvent *event)
{
    if (event->type() == QEvent::HoverMove && selectingCenter) {
        QHoverEvent *hoverEvent = static_cast<QHoverEvent *>(event);
        QPoint p = hoverEvent->position().toPoint();
        clearTemporaryLayer(); // Clear the previous temporary circle
        drawTemporaryCircle(p); // Draw the new temporary circle
        return true;
    } else if (event->type() == QEvent::HoverMove && selectingRect) {
        QHoverEvent *hoverEvent = static_cast<QHoverEvent *>(event);
        QPoint p = hoverEvent->position().toPoint();
        clearTemporaryLayer(); // Clear the previous temporary circle
        drawTemporarRectangle(p); // Draw the new temporary circle
        return true;
    }
    return QWidget::event(event);
}

void DrawGeoArea::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    QRect dirtyRect = event->rect();
    painter.drawImage(dirtyRect, image, dirtyRect); // Draw the main image
    painter.drawImage(dirtyRect, tempImage, dirtyRect); // Draw the temporary image

    // Draw the axis and labels
    drawAxis(painter);
    drawAxisLabels(painter);
}

void DrawGeoArea::resizeEvent(QResizeEvent *event)
{
    if (width() > image.width() || height() > image.height()) {
        int newWidth = qMax(width() + 128, image.width());
        int newHeight = qMax(height() + 128, image.height());
        resizeImage(&image, QSize(newWidth, newHeight));
        resizeImage(&tempImage, QSize(newWidth, newHeight)); // Resize the temporary image
        update();
    }
    QWidget::resizeEvent(event);
}

void DrawGeoArea::drawLineTo(const QPoint &endPoint, QImage *img, QColor PenColor)
{   
    if (PenColor == nullptr) {
        PenColor = myPenColor;
    }
    if (img == nullptr) {
        img = &image;
    }
    QPainter painter(img);
    painter.setPen(QPen(PenColor, myPenWidth, Qt::SolidLine, Qt::RoundCap,
                        Qt::RoundJoin));
    painter.drawLine(lastPoint, endPoint);
    modified = true;

    int rad = (myPenWidth / 2) + 2;
    update(QRect(lastPoint, endPoint).normalized()
                                     .adjusted(-rad, -rad, +rad, +rad));
    lastPoint = endPoint;
}

void DrawGeoArea::drawPoint(const QPoint &P, QImage *img){
    if (img == nullptr) {
        img = &image;
    }
    QPainter painter(img);
    painter.setBrush(Qt::black);
    painter.drawEllipse(P, myPenWidth, myPenWidth);
    modified = true;
    update();
}

void DrawGeoArea::resizeImage(QImage *image, const QSize &newSize)
{
    if (image->size() == newSize)
        return;

    QImage newImage(newSize, image->format());
    newImage.fill(image->format() == QImage::Format_RGB32 ? qRgb(255, 255, 255) : Qt::transparent);
    QPainter painter(&newImage);
    painter.drawImage(QPoint(0, 0), *image);
    *image = newImage;
}

void DrawGeoArea::drawShape(){
    // make last point the geo->plot_points[0]
    int scale = std::max(width(), height());
    lastPoint = {(int)(geo->plot_points[0][0]*scale), (int)((1.0-geo->plot_points[0][1])*scale)};
    firstPoint = lastPoint;
    for (int i = 1; i < geo->plot_points.size(); i++){
        QPoint p = {(int)(geo->plot_points[i][0]*scale), (int)((1.0-geo->plot_points[i][1])*scale)};
        drawLineTo(p, &image, Qt::black);
    }
    //drawLineTo(firstPoint);
    update();
}

void DrawGeoArea::drawSpline(){
    if (geo->points.size()-geo->firstPointIndex <= 1){return;}
    geo->addSpline(geo->firstPointIndex, geo->points.size());

    drawShape();
    clearTemporaryLayer();
}

void DrawGeoArea::drawBSpline(){
    if (geo->points.size()-geo->firstPointIndex <= 1){return;}
    geo->addBSpline(geo->firstPointIndex, geo->points.size());

    drawShape();
    clearTemporaryLayer();
}

#ifdef USE_GEO
void DrawGeoArea::drawBezier(){
    geo->addBezier(geo->firstPointIndex, geo->points.size());

    drawShape();
    clearTemporaryLayer();
}
#endif

void DrawGeoArea::drawCircle()
{
    bool ok;
    rx = QInputDialog::getDouble(this, tr("Input Radius"),
                                           tr("Radius:"), 0.10, 0, 100, 2, &ok);
    if (!ok) return;
    ry = rx;

    selectingCenter = true;
}

void DrawGeoArea::drawEllipse()
{
    bool ok;
    rx = QInputDialog::getDouble(this, tr("Input Major-Radius"),
                                           tr("Major-Radius:"), 0.10, 0, 100, 2, &ok);
    if (!ok) return;

    ry = QInputDialog::getDouble(this, tr("Input Minor-Radius"),
                                           tr("Minor-Radius:"), 0.10, 0, 100, 2, &ok);
    if (!ok) return;
    double rM = std::max(rx, ry);
    double rm = std::min(rx, ry);

    rx = rM;
    ry = rm;

    selectingCenter = true;
}

void DrawGeoArea::drawRectangle() {
    bool ok;
    rx = QInputDialog::getDouble(this, tr("Input Length"),
                                           tr("Length:"), 0.10, 0, 100, 2, &ok);
    if (!ok) return;

    ry = QInputDialog::getDouble(this, tr("Input Height"),
                                           tr("Height:"), 0.10, 0, 100, 2, &ok);
    if (!ok) return;

    selectingRect = true;
}

void DrawGeoArea::FuseAll(){
    geo->FuseOverlapping(exterior_geom);
    tempImage.fill(Qt::transparent);
    image.fill(qRgb(255, 255, 255));
    drawgeometry();
    clearTemporaryLayer();
}

void DrawGeoArea::MakeExteriorRectangle(){
    // Make rectangle that is screen size
    double scale = (double)std::max(width(), height());
    double w = (double)width()/scale; 
    double cx = w/2.0;
    double dx = cx;
    double h = (double)height()/scale;
    double cy = 1.0 - h/2.0;
    double dy = 1.0 - cy;
    geo->MakeRectangleAndCut(cx, cy, dx-0.005, dy-0.005);
    tempImage.fill(Qt::transparent);
    image.fill(qRgb(255, 255, 255));
    drawgeometry();
    clearTemporaryLayer();
}

void DrawGeoArea::MakeExteriorGeometry(){
    if (exterior_geom){return;}

    exterior_geom = true;
    MakeExteriorRectangle();
}

void DrawGeoArea::drawAxis(QPainter &painter)
{
    painter.setPen(QPen(Qt::black, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

    // Draw X axis
    painter.drawLine(0, height() / 2, width(), height() / 2);

    // Draw Y axis
    painter.drawLine(width() / 2, 0, width() / 2, height());
}

void DrawGeoArea::drawAxisLabels(QPainter &painter)
{
    painter.setPen(QPen(Qt::black));
    painter.setFont(QFont("Arial", 10));
    double scale = (double)std::max(width(), height());

    // Draw X axis labels
    for (int i = 0; i <= width(); i += width() / 10) {
        int x = i - width();
        painter.drawText(i, height() / 2 + 15, QString::number(x/scale));
    }

    // Draw Y axis labels
    for (int i = 0; i <= height(); i += height() / 10) {
        int y = height() - i;
        painter.drawText(width() / 2 + 5, i, QString::number(y/scale));
    }
}


void DrawGeoArea::drawgeometry(){
    // draw all goemetries in geo filled in with
    for (int i = 0; i < geo->curveTags.size(); i++){
        int stag = geo->curveTags[i];
        std::vector<double> min, max;
        gmsh::model::getParametrizationBounds(1, stag, min, max);

        std::vector<double> params;
        double h = (max[0]-min[0])/500.0;
        for (double t = min[0]; t < max[0]; t += h){
            params.push_back(t);
        }
        std::vector<double> coords;

        gmsh::model::getValue(1,stag,params,coords);
        geo->plot_points.clear();
        for (int i = 0; i < coords.size(); i+=3){
            geo->plot_points.push_back({coords[i], coords[i+1]});
        }
        drawShape();
    }
}