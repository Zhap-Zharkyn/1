
#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include "mylib.h"

class Window : public QWidget
{
  Q_OBJECT

private:
	int n;
	int pos;
	double a;
	double b;
	double A;
	double B;
	int num;
public:
	first v;
	second w;
	Window (QWidget *parent, double x1, double x2, int n1, int j, FILE *file);
	QSize minimumSizeHint () const;
	QSize sizeHint () const;
	int parse_command_line (int argc, char *argv[]);

public slots:
	void method1 ();
	void method2 ();
	void all ();
	void disc ();
	void minus ();
	void plus ();
	void zoom1 ();
	void zoom2 ();
	void reload ();

protected:
	void paintEvent (QPaintEvent *event);
};

#endif
