#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QAction>
#include <QMenuBar>
#include <QMessageBox>

#include "window.h"

int main(int argc, char** argv) { 
	FILE *file;
	int i, j = 0, n, c;
	double a, b;
	
	file = NULL;
	cout<<"1.Formula\n2.File"<<endl;
	cin>>i;
	
	if (i == 1) {
		cout<<"[a, b]:\n"; 		cin>>a>>b;
		cout<<"n:\n"; 	cin>>n;
		cout<<"f:"<<endl;		cin>>j;
	}
	else if (i == 2) {
			file = fopen("input.txt", "r");
			c = fscanf(file, "%d", &n);
			cout<<c<<endl;
		}
		else {
			cout<<"ERROR!"<<endl;
			return 0;
		}
		
	if (n < 3) {
		cout<<"Слишком маленькое n"<<endl;
		return 0;
	}
	if (a >= b) {
		cout<<"b <= a"<<endl;
		return 0;
	}
	
	QApplication app (argc, argv);
	QMainWindow *window = new QMainWindow;
	QMenuBar *tool_bar = new QMenuBar (window);
	Window *graph_area = new Window (window, a, b, n, j, file);
	QAction *action;
	
	 if (graph_area->parse_command_line (argc, argv))
    {
      QMessageBox::warning (0, "Wrong input arguments!", 
                            "Wrong input arguments!");
      return -1;
    }
      
    
    action = tool_bar->addAction ("&First method", graph_area, SLOT (method1()));
    action->setShortcut (QString ("Ctrl+C"));
    action = tool_bar->addAction ("&Second method", graph_area, SLOT (method2()));
    action->setShortcut (QString ("Ctrl+V"));
    action = tool_bar->addAction ("&Discrepancy", graph_area, SLOT (disc()));
    action->setShortcut (QString ("Ctrl+Z"));
    action = tool_bar->addAction ("&All", graph_area, SLOT (all()));
    action->setShortcut (QString ("Ctrl+A"));
    action = tool_bar->addAction ("&+const", graph_area, SLOT (plus()));
    action->setShortcut (QString ("Ctrl+S"));
    action = tool_bar->addAction ("&-const", graph_area, SLOT (minus()));
    action->setShortcut (QString ("Ctrl+D"));
    action = tool_bar->addAction ("&Zoom in", graph_area, SLOT (zoom1()));
    action->setShortcut (QString ("Ctrl+N"));
	action = tool_bar->addAction ("&Zoom out", graph_area, SLOT (zoom2()));
    action->setShortcut (QString ("Ctrl+M"));
    
    if (i == 1) {
		action = tool_bar->addAction ("&Resize", graph_area, SLOT (reload()));
		action->setShortcut (QString ("Ctrl+F"));
	}
	action = tool_bar->addAction ("E&xit", window, SLOT (close ()));
	action->setShortcut (QString ("Ctrl+X"));
	
	tool_bar->setMaximumHeight (30);

	window->setMenuBar (tool_bar);
	window->setCentralWidget (graph_area);
	window->setWindowTitle ("Graph");

	window->show ();
	return app.exec ();
}
