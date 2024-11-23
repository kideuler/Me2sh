#ifndef CONSOLEOUTPUT_HPP
#define CONSOLEOUTPUT_HPP

#include <Python.h>
#include <QTextEdit>
#include <QWidget>

class ConsoleOutput : public QTextEdit {
    Q_OBJECT

public:
    ConsoleOutput(QWidget *parent = nullptr);
    void clear() { QTextEdit::clear(); }
    void printString(const QString &text);
    void executeCommand(const QString &command);
    void appendPrompt();

protected:
    void keyPressEvent(QKeyEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;

private:
    
    void showPreviousCommand();

    QString prompt = ">>> ";
    QString currentCommand;
    QStringList commandHistory;
    int historyIndex = -1;
};

#endif // CONSOLEOUTPUT_HPP