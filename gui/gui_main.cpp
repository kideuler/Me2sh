#include <QApplication>

#include "mainwindow.hpp"

int main(int argc, char *argv[])
{   

    // Create the application
    QApplication app(argc, argv);

    // Create the main window
    MainWindow window;

    // Show the main window
    window.show();

    // Run the application
    return app.exec();
}