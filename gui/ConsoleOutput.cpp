#include "ConsoleOutput.hpp"
#include <QVBoxLayout>
#include <QScrollBar>
#include <QKeyEvent>
#include <QFont>
#include <string>
#include <iostream>

class PythonInterpreter {
public:
    PythonInterpreter() {
        Py_Initialize();
        main_module = PyImport_AddModule("__main__");
        main_dict = PyModule_GetDict(main_module);

        // Redirect stdout and stderr to custom streams
        PyRun_SimpleString(
            "import sys\n"
            "class StreamCapturer:\n"
            "    def __init__(self):\n"
            "        self.data = ''\n"
            "    def write(self, s):\n"
            "        self.data += s\n"
            "    def flush(self):\n"
            "        pass\n"
            "sys.stdout = StreamCapturer()\n"
            "sys.stderr = StreamCapturer()\n"
        );
    }

    ~PythonInterpreter() {
        Py_Finalize();
    }

    QString executeCommand(const QString &command) {
        PyObject *result = PyRun_String(command.toStdString().c_str(), Py_single_input, main_dict, main_dict);
        if (result == nullptr) {
            PyErr_Print();
        } else {
            Py_DECREF(result);
        }

        // Retrieve the captured output
        PyObject *sys_module = PyImport_ImportModule("sys");
        PyObject *stdout_obj = PyObject_GetAttrString(sys_module, "stdout");
        PyObject *stderr_obj = PyObject_GetAttrString(sys_module, "stderr");
        PyObject *stdout_data = PyObject_GetAttrString(stdout_obj, "data");
        PyObject *stderr_data = PyObject_GetAttrString(stderr_obj, "data");

        QString output = QString::fromUtf8(PyUnicode_AsUTF8(stdout_data)) + QString::fromUtf8(PyUnicode_AsUTF8(stderr_data));

        // Clear the captured output
        PyObject_SetAttrString(stdout_obj, "data", PyUnicode_FromString(""));
        PyObject_SetAttrString(stderr_obj, "data", PyUnicode_FromString(""));

        Py_DECREF(sys_module);
        Py_DECREF(stdout_obj);
        Py_DECREF(stderr_obj);
        Py_DECREF(stdout_data);
        Py_DECREF(stderr_data);

        return output;
    }

private:
    PyObject *main_module;
    PyObject *main_dict;
};

ConsoleOutput::ConsoleOutput(QWidget *parent) : QTextEdit(parent) {
    setReadOnly(false);
    setAcceptRichText(false);
    setUndoRedoEnabled(false);
    // Set the background color to white and text color to black
    setStyleSheet("background-color: white; color: black;");

    // Set the font
    QFont font("Courier", 14);
    setFont(font);

    appendPrompt();
}

void ConsoleOutput::keyPressEvent(QKeyEvent *event) {
    if (event->key() == Qt::Key_Return || event->key() == Qt::Key_Enter) {
        moveCursor(QTextCursor::End);
        QString command = toPlainText().section('\n', -1).mid(prompt.length());
        executeCommand(command);
        commandHistory.append(command);
        historyIndex = commandHistory.size();
        appendPrompt();
    } else if (event->key() == Qt::Key_Up) {
        showPreviousCommand();
    } else if (event->key() == Qt::Key_Backspace || event->key() == Qt::Key_Delete) {
        QTextCursor cursor = textCursor();
        if (cursor.position() > toPlainText().lastIndexOf(prompt) + prompt.length()) {
            QTextEdit::keyPressEvent(event);
        }
    } else if (event->key() == Qt::Key_Left || event->key() == Qt::Key_Right) {
        QTextCursor cursor = textCursor();
        if (cursor.position() > toPlainText().lastIndexOf(prompt) + prompt.length()) {
            QTextEdit::keyPressEvent(event);
        }
    } else {
        QTextEdit::keyPressEvent(event);
    }
}

void ConsoleOutput::mousePressEvent(QMouseEvent *event) {
    // Ignore mouse press events to prevent cursor movement
    event->ignore();
}

void ConsoleOutput::mouseReleaseEvent(QMouseEvent *event) {
    // Ignore mouse release events to prevent cursor movement
    event->ignore();
}

void ConsoleOutput::executeCommand(const QString &command) {
    static PythonInterpreter interpreter;
    QString output = interpreter.executeCommand(command);
    insertPlainText("\n" + output);
    moveCursor(QTextCursor::End);
}

void ConsoleOutput::appendPrompt() {
    insertPlainText("" + prompt);
    moveCursor(QTextCursor::End);
}

void ConsoleOutput::showPreviousCommand() {
    if (historyIndex > 0) {
        historyIndex--;
        QString previousCommand = commandHistory.at(historyIndex);
        QTextCursor cursor = textCursor();
        cursor.movePosition(QTextCursor::End, QTextCursor::MoveAnchor);
        cursor.select(QTextCursor::LineUnderCursor);
        cursor.removeSelectedText();
        cursor.insertText(prompt + previousCommand);
        setTextCursor(cursor);
    }
}

void ConsoleOutput::printString(const QString &text) {
    insertPlainText(text);
    moveCursor(QTextCursor::End);
}