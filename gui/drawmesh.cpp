#include "drawmesh.hpp"

#include <QMouseEvent>
#include <QPainter>
#include <QApplication>
#include <QInputDialog>
#include <QHoverEvent>

DrawMeshArea::DrawMeshArea(Me2sh_Geometry *geometry, QWidget *parent)
    : QWidget(parent), geo(geometry)
{
    setAttribute(Qt::WA_StaticContents);
    setMouseTracking(true);
    setAttribute(Qt::WA_Hover);
    image = QImage(size(), QImage::Format_RGB32);
    image.fill(qRgb(255, 255, 255));
    //tempImage = QImage(size(), QImage::Format_ARGB32_Premultiplied);
    //tempImage.fill(Qt::transparent);
}

void DrawMeshArea::clearImage()
{
    tempImage.fill(Qt::transparent);
    image.fill(qRgb(255, 255, 255));
    modified = true;
    geo->clear();
    update();
}

void DrawMeshArea::clearTemporaryLayer()
{
    tempImage.fill(Qt::transparent); // Clear only the temporary image
    geo->cleartemp();
    update();
}

void DrawMeshArea::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    QRect dirtyRect = event->rect();
    painter.drawImage(dirtyRect, image, dirtyRect); // Draw the main image
    painter.drawImage(dirtyRect, tempImage, dirtyRect); // Draw the temporary image

    // Draw the axis and labels
    drawAxis(painter);
    drawAxisLabels(painter);
}

void DrawMeshArea::resizeEvent(QResizeEvent *event)
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

void DrawMeshArea::drawLineTo(const QPoint &endPoint, QImage *img, QColor PenColor)
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

void DrawMeshArea::resizeImage(QImage *image, const QSize &newSize)
{
    if (image->size() == newSize)
        return;

    QImage newImage(newSize, image->format());
    newImage.fill(image->format() == QImage::Format_RGB32 ? qRgb(255, 255, 255) : Qt::transparent);
    QPainter painter(&newImage);
    painter.drawImage(QPoint(0, 0), *image);
    *image = newImage;
}

void DrawMeshArea::drawShape(){
    // make last point the geo->plot_points[0]
    int scale = std::max(width(), height());
    lastPoint = {(int)(geo->plot_points[0][0]*scale), (int)((1.0-geo->plot_points[0][1])*scale)};
    firstPoint = lastPoint;
    for (int i = 1; i < geo->plot_points.size(); i++){
        QPoint p = {(int)(geo->plot_points[i][0]*scale), (int)((1.0-geo->plot_points[i][1])*scale)};
        drawLineTo(p, &image, Qt::black);
    }
    drawLineTo(firstPoint);
    update();
}

void DrawMeshArea::drawAxis(QPainter &painter)
{
    painter.setPen(QPen(Qt::black, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

    // Draw X axis
    painter.drawLine(0, height() / 2, width(), height() / 2);

    // Draw Y axis
    painter.drawLine(width() / 2, 0, width() / 2, height());
}

void DrawMeshArea::drawAxisLabels(QPainter &painter)
{
    painter.setPen(QPen(Qt::black));
    painter.setFont(QFont("Arial", 10));
    double scale = (double)std::max(width(), height());

    // Draw X axis labels
    for (int i = 0; i <= width(); i += width() / 10) {
        int x = i - width() / 2;
        painter.drawText(i, height() / 2 + 15, QString::number(x/scale));
    }

    // Draw Y axis labels
    for (int i = 0; i <= height(); i += height() / 10) {
        int y = height() / 2 - i;
        painter.drawText(width() / 2 + 5, i, QString::number(y/scale));
    }
}

void DrawMeshArea::drawgeometry(){
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

void DrawMeshArea::chanegeMeshSize(){
    bool ok;
    double h = QInputDialog::getDouble(this, tr("Input Mesh Size"),
                                           tr("Mesh Size:"), 0.05, 0, 100, 2, &ok);
    if (!ok) return;
    this->h_target = h;
    gmsh::option::setNumber("Mesh.MeshSizeMin", h_target);
    gmsh::option::setNumber("Mesh.MeshSizeMax", 1.3*h_target);
}

void DrawMeshArea::changeElementType(){
    bool ok;
    QStringList items;
    items << tr("Triangles") << tr("Quadrilaterals");

    QString item = QInputDialog::getItem(this, tr("Select Element Type"),
                                         tr("Element Type:"), items, 0, false, &ok);
    if (!ok) return;

    if (item == "Triangles") {
        this->ElementType = 1;
    } else if (item == "Quadrilaterals") {
        this->ElementType = 2;
    }
}