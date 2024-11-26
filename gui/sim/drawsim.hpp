#ifndef DRAWSIM_H
#define DRAWSIM_H

#include "ConsoleOutput.hpp"

#include <QColor>
#include <QImage>
#include <QTimer>
#include <QPoint>
#include <QWidget>
#include <QMouseEvent>
#include <QPainter>
#include <QLinearGradient>
#include <QApplication>
#include <QInputDialog>
#include <QHoverEvent>
#include <QPushButton>
#include <QProgressBar>
#include <QSlider>
#include <QLabel>
#include <QHBoxLayout>
#include <chrono>
#include <memory>

#include "geometry.hpp"
#include "meshing.hpp"

#include "simulation.hpp"

class DrawSimArea : public QWidget
{
    Q_OBJECT
public:
    DrawSimArea(ConsoleOutput *PyTerm, std::shared_ptr<Me2sh_Geometry> geometry, std::shared_ptr<Me2sh_Mesh> Mesh, std::shared_ptr<Me2sh_Simulation> simulation, QWidget *parent = nullptr);
    void clearImage();
    void setTimeSim(bool value, int num_steps, double dt);

    std::shared_ptr<Me2sh_Geometry> geo;
    std::shared_ptr<Me2sh_Mesh> mesh;
    std::shared_ptr<Me2sh_Simulation> sim;
    ConsoleOutput *PythonTerminal;

    mfem::Mesh *mfem_mesh;

    void drawgeometry();
    void displayMesh();

    int time_steps = 100;
    double dt = 0.02;

protected:
    void paintEvent(QPaintEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;

private slots:
    void sliderValueChanged(int value);

private:
    void resizeImage(QImage *image, const QSize &newSize);
    void drawShape();
    void clearTemporaryLayer();
    void drawAxis(QPainter &painter); 
    void drawAxisLabels(QPainter &painter); 
    void drawPoint(const QPoint &P, QImage *img = nullptr);
    void drawPolygon(const std::vector<QPoint> &points, QImage *img = nullptr);

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

    QPushButton *playButton;
    QSlider *slider;
    QLabel *sliderLabel;
    QHBoxLayout *topLayout;

    bool TimeSim = false;
};


#endif // DRAWSIM_H