#include <QApplication>
#include <QSettings>
#include <QColor>
#include <QHBoxLayout>
#include <QMessageBox>
#include <QFontDatabase>

#include "ConsoleOutput.hpp"

ConsoleOutput::ConsoleOutput(QWidget *parent) : QWidget(parent)
{
    msgBox = new QTextEdit(this);
    msgBox->setReadOnly(true);
    msgBox->setLineWrapMode(QTextEdit::NoWrap);
    msgBox->setStyleSheet("background-color: white; color: black;");
    msgBox->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    QHBoxLayout *layout = new QHBoxLayout;
    layout->addWidget(msgBox);
    setLayout(layout);
}

void ConsoleOutput::addMessage(const QString &text)
{
    msgBox->append(">>> " + text);
}