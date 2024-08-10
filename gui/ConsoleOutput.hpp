#ifndef CONSOLEOUTPUT_HPP
#define CONSOLEOUTPUT_HPP

#include <QWidget>
#include <QTextEdit>

class ConsoleOutput : public QWidget {
    Q_OBJECT

    public:
        ConsoleOutput(QWidget *parent = nullptr);
        void clear() {msgBox->clear();}
    
    public Q_SLOTS:

        void addMessage(const QString &text);

    private:
        QTextEdit *msgBox = nullptr;
};

#endif