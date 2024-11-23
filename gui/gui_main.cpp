#include "mainwindow.hpp"
#include <QApplication>

int main(int argc, char *argv[])
{   

    // Create the application
    QApplication app(argc, argv);
    gmsh::initialize();
    gmsh::model::add("me2sh");

    // Create the main window
    MainWindow window;

    // Show the main window
    window.show();
    window.clearScreen(); // ensure the screen is clear
    // Run the application
    window.PythonTerminal->executeCommand("print('Python interpreter initialized')");
    window.PythonTerminal->appendPrompt();
    app.exec();
    gmsh::finalize();
    return 0;
}