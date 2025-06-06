#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QAction>
#include <QMenuBar>
#include <QMessageBox>
#include "window.h"

int main (int argc, char *argv[])
{
  QApplication app (argc, argv);

  QMainWindow *window = new QMainWindow;
  QMenuBar *tool_bar = new QMenuBar (window);
  Window *graph_area = new Window (window);
  QAction *action = new QAction();

  if (graph_area->parse_command_line (argc, argv))
    {
      QMessageBox::warning (0, "Wrong input arguments!",
                            "Usage: program a b n k\n"
                            "a - left bound (double)\n"
                            "b - right bound (double)\n"
                            "n - initial number of interpolation points (int)\n"
                            "k - initial function index [0-6] (int)");
      delete window;
      delete graph_area;
      delete tool_bar;
      delete action;
      return -1;
    }
  else
    {
      graph_area->set_func_by_id();
    }

  action = tool_bar->addAction ("&Change function (0)", graph_area, SLOT (change_func ()));
  action->setShortcut (QString ("0"));

  action = tool_bar->addAction ("Change &display mode (1)", graph_area, SLOT (change_display_mode()));
  action->setShortcut (QString ("1"));

  action = tool_bar->addAction ("E&xit", window, SLOT (close ()));
  action->setShortcut (QString ("Ctrl+X"));

  tool_bar->setMaximumHeight (30);

  window->setMenuBar (tool_bar);
  window->setCentralWidget (graph_area);
  window->setWindowTitle ("Approximation 1D");

  window->show ();
  app.exec ();
  delete window;
  return 0;
}
