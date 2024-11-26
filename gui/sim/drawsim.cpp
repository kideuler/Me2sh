#include "drawsim.hpp"

DrawSimArea::DrawSimArea(ConsoleOutput *PyTerm, std::shared_ptr<Me2sh_Geometry> geometry, std::shared_ptr<Me2sh_Mesh> Mesh, std::shared_ptr<Me2sh_Simulation> simulation, QWidget *parent)
    : QWidget(parent), PythonTerminal(PyTerm), geo(geometry), mesh(Mesh), sim(simulation)
{
    setAttribute(Qt::WA_StaticContents);
    setMouseTracking(true);
    setAttribute(Qt::WA_Hover);
    // Initialize images with the correct size
    image = QImage(size(), QImage::Format_RGB32);
    image.fill(qRgb(255, 255, 255)); // Ensure the background is white

    //tempImage = QImage(size(), QImage::Format_ARGB32_Premultiplied);
    //tempImage.fill(Qt::transparent); // Ensure the temporary image is transparent

    // Create the play button and progress bar and label
    playButton = new QPushButton("Play", this);
    playButton->setStyleSheet("color: black;");
    slider = new QSlider(Qt::Horizontal, this);
    slider->setRange(0, time_steps);
    slider->setValue(0);
    sliderLabel = new QLabel("Time: 0", this);
    sliderLabel->setStyleSheet("color: black;");

    // Create the top layout and add the play button and progress bar
    topLayout = new QHBoxLayout();
    topLayout->addWidget(playButton);
    topLayout->addWidget(slider);
    topLayout->addWidget(sliderLabel);

    // Set the layout for the widget
    QVBoxLayout *mainLayout = new QVBoxLayout(this);
    mainLayout->addLayout(topLayout);
    mainLayout->addStretch();
    setLayout(mainLayout);

    setTimeSim(true, time_steps, dt);

    // Connect the slider to the sliderValueChanged slot
    connect(slider, &QSlider::valueChanged, this, &DrawSimArea::sliderValueChanged);

    update();
}

void DrawSimArea::clearImage()
{
    tempImage.fill(Qt::transparent);
    image.fill(qRgb(255, 255, 255));
    modified = true;
    update();
}

void DrawSimArea::setTimeSim(bool value, int num_steps, double dt)
{
    TimeSim = value;
    playButton->setEnabled(TimeSim);
    slider->setEnabled(TimeSim);
    slider->setRange(0, num_steps);
    sliderLabel->setText("Time: 0");
}

void DrawSimArea::sliderValueChanged(int value)
{
    sliderLabel->setText(QString("Time: %1").arg(value*dt));
}

void DrawSimArea::clearTemporaryLayer()
{
    tempImage.fill(Qt::transparent); // Clear only the temporary image
    update();
}

void DrawSimArea::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    QRect dirtyRect = event->rect();
    painter.drawImage(dirtyRect, image, dirtyRect); // Draw the main image
    painter.drawImage(dirtyRect, tempImage, dirtyRect); // Draw the temporary image

    // Draw the axis and labels
    drawAxis(painter);
    drawAxisLabels(painter);
}

void DrawSimArea::resizeEvent(QResizeEvent *event)
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

void DrawSimArea::drawLineTo(const QPoint &endPoint, QImage *img, QColor PenColor)
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

void DrawSimArea::resizeImage(QImage *image, const QSize &newSize)
{
    if (image->size() == newSize){return;}
    QImage newImage(newSize, image->format());
    newImage.fill(image->format() == QImage::Format_RGB32 ? qRgb(255, 255, 255) : Qt::transparent);
    QPainter painter(&newImage);
    painter.drawImage(QPoint(0, 0), *image);
    *image = newImage;
}

void DrawSimArea::drawShape(){
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

void DrawSimArea::drawAxis(QPainter &painter)
{
    painter.setPen(QPen(Qt::black, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

    // Draw X axis
    painter.drawLine(0, height() / 2, width(), height() / 2);

    // Draw Y axis
    painter.drawLine(width() / 2, 0, width() / 2, height());
}

void DrawSimArea::drawAxisLabels(QPainter &painter)
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

void DrawSimArea::drawgeometry(){
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

void DrawSimArea::drawPoint(const QPoint &P, QImage *img){
    if (img == nullptr) {
        img = &image;
    }
    QPainter painter(img);
    painter.setBrush(Qt::black);
    painter.drawEllipse(P, myPenWidth, myPenWidth);
    modified = true;
    update();
}

void DrawSimArea::drawPolygon(const std::vector<QPoint> &points, QImage *img) {
    if (img == nullptr) {
        img = &image;
    }
    QPainter painter(img);
    painter.setPen(QPen(Qt::black, myPenWidth, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    painter.drawPolygon(points.data(), points.size());
    modified = true;
    update();
}

void DrawSimArea::displayMesh(){
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