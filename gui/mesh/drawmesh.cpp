#include "drawmesh.hpp"

#include <QMouseEvent>
#include <QPainter>
#include <QLinearGradient>
#include <QApplication>
#include <QInputDialog>
#include <QHoverEvent>
#include <chrono>

DrawMeshArea::DrawMeshArea(ConsoleOutput *PyTerm, std::shared_ptr<Me2sh_Geometry> geometry, std::shared_ptr<Me2sh_Mesh> mesh, QWidget *parent)
    : QWidget(parent), PythonTerminal(PyTerm), geo(geometry), mesh(mesh), animationPhase(0.0)
{
    setAttribute(Qt::WA_StaticContents);
    setMouseTracking(true);
    setAttribute(Qt::WA_Hover);
    image = QImage(size(), QImage::Format_RGB32);
    image.fill(qRgb(255, 255, 255));
    tempImage = QImage(size(), QImage::Format_ARGB32_Premultiplied);
    tempImage.fill(Qt::transparent);

    // Set up the animation timer
    animationTimer = new QTimer(this);
    connect(animationTimer, &QTimer::timeout, this, &DrawMeshArea::updateAnimation);
}

void DrawMeshArea::startAnimation()
{
    animationTimer->start(16); // Update at ~60 FPS
}

void DrawMeshArea::clearImage()
{
    tempImage.fill(Qt::transparent);
    image.fill(qRgb(255, 255, 255));
    modified = true;
    hasMesh = false;
    //geo->clear();
    mesh->clear();
    update();
}

void DrawMeshArea::clearTemporaryLayer()
{
    tempImage.fill(Qt::transparent); // Clear only the temporary image
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
    if (image->size() == newSize){return;}
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
        int x = i;
        painter.drawText(i, height() / 2 + 15, QString::number(x/scale));
    }

    // Draw Y axis labels
    for (int i = 0; i <= height(); i += height() / 10) {
        int y = i;
        painter.drawText(width() / 2 + 5, i, QString::number(1.0 - y/scale));
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
    bool ok = false;
    double h = 0.05;

    try {
        h = QInputDialog::getDouble(this, tr("Input Mesh Size"),tr("Mesh Size:"), h_target, 0, 100, 4, &ok);
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        h = 0.05;
    }  
    
    if (!ok){ 
        std::cerr << "Not OKAY" << std::endl;
        return;
    }
    h_target = h;
}

void DrawMeshArea::changeElementType(){
    bool ok;
    QStringList items;
    items.append(tr("Triangles"));
    items.append(tr("Quadrilaterals"));

    QString item = QInputDialog::getItem(this, tr("Select Element Type"),
                                         tr("Element Type:"), items, ElementType-1, true, &ok);
    if (!ok) {
        return;
    }
    int index = items.indexOf(item);

    std::cout << item.toStdString() << " " << index << std::endl;
    int type = 0;
    if (index == 0 ) {
        type = 1;
    }
    if (index == 1) {
        type = 2;
    }
    ElementType = type;
    std::cout << "ElementType " << ElementType << std::endl;
}

void DrawMeshArea::changeMeshAlgorithm(){
    bool ok;
    QStringList items_tri;
    items_tri.append(tr("Frontal-Delaunay"));
    items_tri.append(tr("Delaunay"));
    
    QStringList items_quad;
    items_quad.append(tr("Blossom"));
    items_quad.append(tr("Frontal-Delaunay for Quads"));
    items_quad.append(tr("Packing-Parallellelogram"));
    items_quad.append(tr("Crossfield"));

    QString item;
    if (ElementType == 1){
        item = QInputDialog::getItem(this, tr("Select Mesh Algorithm"),
                                         tr("Mesh Algorithm:"), items_tri, GmshMeshAlgorithm-5, true, &ok);
    } else {
        item = QInputDialog::getItem(this, tr("Select Mesh Algorithm"),
                                         tr("Mesh Algorithm:"), items_quad, GmshMeshAlgorithm-5, true, &ok);
    }
    if (!ok) {
        return;
    }
    if (ElementType == 1){
        int index = items_tri.indexOf(item);
        TriMeshAlgorithm = index;
    } else {
        int index = items_quad.indexOf(item);
        QuadMeshAlgorithm = index;
    }

}

void DrawMeshArea::generateMesh(){

    clearImage();
    int TriAlgos[2] = {6,5};
    int QuadAlgos[4] = {6,8,9,11};
    if (ElementType == 1){
        GmshMeshAlgorithm = TriAlgos[TriMeshAlgorithm];
    } else {
        GmshMeshAlgorithm = QuadAlgos[QuadMeshAlgorithm];
    }
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    mesh->generate(geo->planeTags, GmshMeshAlgorithm, ElementType, GmshRecombinationAlgorithm, h_target);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    PythonTerminal->printString("Mesh Generated in " + QString::number(elapsed_seconds.count()) + " seconds");
    if (ElementType == 1){
        int numTri = 0;
        for (int i = 0; i < mesh->triMeshes.size(); i++){
            numTri += mesh->triMeshes[i].first.size();
        }
        int numNodes = 0;
        for (int i = 0; i < mesh->triMeshes.size(); i++){
            numNodes += mesh->triMeshes[i].second.size();
        }

        PythonTerminal->printString("Number of Triangles: " + QString::number(numTri) + " Number of Nodes: " + QString::number(numNodes));
    } else {
        int numQuad = 0;
        for (int i = 0; i < mesh->quadMeshes.size(); i++){
            numQuad += mesh->quadMeshes[i].first.size();
        }
        int numNodes = 0;
        for (int i = 0; i < mesh->quadMeshes.size(); i++){
            numNodes += mesh->quadMeshes[i].second.size();
        }

        PythonTerminal->printString("Number of Quads: " + QString::number(numQuad) + " Number of Nodes: " + QString::number(numNodes));
    }
    PythonTerminal->appendPrompt();

    displayMesh();
    hasMesh = true;
}

void DrawMeshArea::drawPoint(const QPoint &P, QImage *img){
    if (img == nullptr) {
        img = &image;
    }
    QPainter painter(img);
    painter.setBrush(Qt::black);
    painter.drawEllipse(P, myPenWidth, myPenWidth);
    modified = true;
    update();
}

void DrawMeshArea::drawPolygon(const std::vector<QPoint> &points, QImage *img) {
    if (img == nullptr) {
        img = &image;
    }
    QPainter painter(img);
    painter.setPen(QPen(Qt::black, myPenWidth, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    painter.drawPolygon(points.data(), points.size());
    modified = true;
    update();
}

void DrawMeshArea::displayMesh(){
    // first, draw all points in mesh
    int scale = std::max(width(), height());
    for (int i = 0; i < mesh->triMeshes.size(); i++){
        // draw all points in the mesh
        for (int k = 0; k < mesh->triMeshes[i].second.size(); k++){
            QPoint p = {(int)(mesh->triMeshes[i].second[k][0]*scale), (int)((1.0-mesh->triMeshes[i].second[k][1])*scale)};
            drawPoint(p);
        }

        // draw all triangles in the mesh
        for (int k = 0; k < mesh->triMeshes[i].first.size(); k++){
            std::vector<QPoint> points;
            for (int l = 0; l < 3; l++){
                int index = mesh->triMeshes[i].first[k][l];
                QPoint p = {(int)(mesh->triMeshes[i].second[index][0]*scale), (int)((1.0-mesh->triMeshes[i].second[index][1])*scale)};
                points.push_back(p);
            }
            drawPolygon(points);
        }
    }

    for (int i = 0; i < mesh->quadMeshes.size(); i++){
        // draw all points in the mesh
        for (int k = 0; k < mesh->quadMeshes[i].second.size(); k++){
            QPoint p = {(int)(mesh->quadMeshes[i].second[k][0]*scale), (int)((1.0-mesh->quadMeshes[i].second[k][1])*scale)};
            drawPoint(p);
        }

        // draw all triangles in the mesh
        for (int k = 0; k < mesh->quadMeshes[i].first.size(); k++){
            std::vector<QPoint> points;
            for (int l = 0; l < 4; l++){
                int index = mesh->quadMeshes[i].first[k][l];
                QPoint p = {(int)(mesh->quadMeshes[i].second[index][0]*scale), (int)((1.0-mesh->quadMeshes[i].second[index][1])*scale)};
                points.push_back(p);
            }
            drawPolygon(points);
        }
    }

}

void DrawMeshArea::updateAnimation() {
    animationPhase += 0.1; // Increment the phase for the wave
    if (animationPhase > 2 * M_PI) {
        animationPhase -= 2 * M_PI;
    }
    ShowMeshColorGradient();
}

void DrawMeshArea::ShowMeshColorGradient() {
    if (mesh->triMeshes.empty()) return;

    tempImage.fill(Qt::transparent);
    image.fill(qRgb(255, 255, 255));
    update();

    QPainter painter(&image);
    painter.setRenderHint(QPainter::Antialiasing);

    // Find the min and max x coordinates to normalize the gradient
    int minX = std::numeric_limits<int>::max();
    int maxX = std::numeric_limits<int>::lowest();
    int scale = std::max(width(), height());
    for (const auto &v : mesh->triMeshes[0].second) {
        QPoint node = {(int)(v[0]*scale), (int)((1.0 - v[1])*scale)};
        if (node.x() < minX) minX = node.x();
        if (node.x() > maxX) maxX = node.x();
    }

    // Draw each triangle with a color gradient
    for (const auto &triangle : mesh->triMeshes[0].first) {
        QPolygon polygon;
        QVector<QColor> colors;

        for (const auto &index : triangle) {
            const auto &v = mesh->triMeshes[0].second[index];
            QPoint node = {(int)(v[0]*scale), (int)((1.0 - v[1])*scale)};
            
            polygon << node;
            
            // Normalize the x coordinate to [0, 1]
            double normalizedX = double(node.x() - minX) / double(maxX - minX);

            // Apply a sine wave to the normalized x coordinate
            double wave = 0.5 * (1 + std::sin(normalizedX * 2 * M_PI + animationPhase));

            // Interpolate color from red to blue
            int c1 = static_cast<int>(255 * (1 - wave));
            int c2 = static_cast<int>(255 * wave);
            QColor color(c1, c2, 0);
            colors << color;
        }

        // Create a gradient for the triangle
        QLinearGradient gradient(polygon.boundingRect().topLeft(), polygon.boundingRect().bottomRight());
        gradient.setColorAt(0.0, colors[0]);
        gradient.setColorAt(0.5, colors[1]);
        gradient.setColorAt(1.0, colors[2]);

        painter.setBrush(gradient);
        painter.drawPolygon(polygon);
    }
}