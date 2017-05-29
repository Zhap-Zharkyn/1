#include <QPainter>
#include <stdio.h>

#include "window.h"

Window::Window (QWidget *parent, double x1, double x2, int n1, int j, FILE *file) : QWidget (parent), v(n1, 5, j), w(n1, j)  { 
	a = x1; b = x2;
	A = a; B = b;
	n = n1; pos = 0;
	num = j;
	if (j == 0) {
		for (int i = 0; i < n; i++) { 
			j = fscanf(file, "%lf", &v.x[i]);
			j = fscanf(file, "%lf", &v.f[i]);
			w.x[i] = v.x[i];
			w.f[i] = v.f[i];
		}
		fclose(file);
	}
	else {
		int i;
		double h = (b - a) / n;
		v.x[n - 1] = b;
		w.x[n - 1] = b;
		for (i = 0; i < n - 1; i++) {
			v.x[i] = a + i * h;
			w.x[i] = v.x[i];
		}
		for (i = 0; i < n; i++) {
			v.f[i] = formula(num, v.x[i]);
			w.f[i] = v.f[i];
		}
	}
	v.method_init();
	cout<<v.discrepancy_res(width())<<endl;
	w.method_init();
	cout<<w.discrepancy_res(width())<<endl;
	update ();
}

QSize Window::minimumSizeHint () const
{
  return QSize (100, 100);
}

QSize Window::sizeHint () const
{
  return QSize (1000, 1000);
}

int Window::parse_command_line (int argc, char *argv[])
{
  if (argc == 1)
    return 0;

  if (argc == 2)
    return -1;

  if (   sscanf (argv[1], "%lf", &a) != 1
      || sscanf (argv[2], "%lf", &b) != 1
      || b - a < 1.e-6
      || (argc > 3 && sscanf (argv[3], "%d", &n) != 1)
      || n <= 0)
    return -2;

  return 0;
}

/// change current function for drawing
void Window::method1() {
	pos = 1;
	update();
}

void Window::method2() {
	pos = 2;
	update();
}

void Window::all() {
	pos = 3;
	update();
}

void Window::disc() {
	pos = 4;
	update();
}

void Window::minus() {
	int i;
	double s;
	cout<<"Enter i for funstion minus:\n";
	cin>>i;
	cout<<"Enter constant:"<<endl;
	cin>>s;
	v.f[i] -= s;
	w.f[i] -= s;
	v.method_init();
	w.method_init();
	cout<<v.discrepancy_res(width())<<endl;
	cout<<w.discrepancy_res(width())<<endl; 
}

void Window::plus() {
	int i;
	double s;
	cout<<"Enter i for funstion plus:\n";
	cin>>i;
	cout<<"Enter constant:"<<endl;
	cin>>s;
	v.f[i] += s;
	w.f[i] += s;
	v.method_init();
	w.method_init(); 
	cout<<v.discrepancy_res(width())<<endl;
	cout<<w.discrepancy_res(width())<<endl;
}

void Window::reload() {
	int i;
	double h;
	cout<<"Enter new n:\n";
	cin>>n;
	v.resize(n);
	w.resize(n);
	h = (b - a) / n;
	v.x[0] = w.x[0] = a;
	v.x[n - 1] = w.x[n - 1] = b;
	for (i = 1; i < n - 1; i++) {
		v.x[i] = v.x[0] + i * h;
		w.x[i] = v.x[i];
	}
	for (i = 0; i < n; i++) {
		v.f[i] = formula(num, v.x[i]);
		w.f[i] = v.f[i];
	}
	v.method_init();
	w.method_init(); 
	cout<<v.discrepancy_res(width())<<endl;
	cout<<w.discrepancy_res(width())<<endl;
	
}



void Window::zoom1 () { 
	double c = (a + b) / 2;
	A = A + 0.125 * (b - a);
	B = B - 0.125 * (b - a);
	if (fabs(A) < EPS || (A - c) >= 0) {
		A = A - 0.125 * (b - a);
		B = B + 0.125 * (b - a);
		cout<<"ERROR"<<endl;
	}
	cout<<"a = "<<A<<", b = "<<B<<endl;
	update();
}

void Window::zoom2 () {
	if (fabs(a - A) > EPS) {
		A = A - 0.125 * (b - a);
		B = B + 0.125 * (b - a);
	}
	cout<<"a = "<<A<<", b = "<<B<<endl;
	update();
}

/// render graph
void Window::paintEvent (QPaintEvent * ) {  
	QPainter painter (this);
	double x1, x2, y1, y2, z1, z2;
	double max_y, min_y;
	double delta_y, delta_x = (b - a) / width ();
	
	painter.setPen ("black");
	x1 = 0;
	z1 = z2 = 0.0;
	y1 = y2 = 0.0;
	max_y = min_y = 0;
	
	if (pos == 0) 
		for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
			y1 = formula(num, x1);
			if (y1 < min_y) min_y = y1;
			if (y1 > max_y) max_y = y1;
		}
	if (pos <= 2 && pos > 0) {
		for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
			if (pos == 1) y1 = v.method_compute(x1);
			if (pos == 2) y1 = w.method_compute(x1);
			if (y1 < min_y) min_y = y1;
			if (y1 > max_y) max_y = y1;
		}
	}
	else if (pos == 3) 
			for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
				y1 = v.method_compute(x1);
				if (y1 < min_y) min_y = y1;
				if (y1 > max_y) max_y = y1;
				z1 = w.method_compute(x1);
				if (z1 < min_y) min_y = z1;
				if (z1 > max_y) max_y = z1;
			}
		else if (pos == 4) 
				for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
					y1 = v.discrepancy(x1);
					if (y1 < min_y) min_y = y1;
					if (y1 > max_y) max_y = y1;
					z1 = w.discrepancy(x1);
					if (z1 < min_y) min_y = z1;
					if (z1 > max_y) max_y = z1;
				}
	if (x1 < b && pos == 4) {
		x1 = b;
		y1 = v.discrepancy(x1);
		if (y1 < min_y) min_y = y1;
		if (y1 > max_y) max_y = y1;
		z1 = w.discrepancy(x1);
		if (z1 < min_y) min_y = z1;
		if (z1 > max_y) max_y = z1;
	}
	delta_y = 0.01 * (max_y - min_y);
	min_y -= delta_y; max_y += delta_y;
	painter.save ();
	
	painter.translate (0.5 * width (), 0.5 * height ());
	painter.scale (width () / (B - A), -height () / (max_y - min_y));
	painter.translate (-0.5*(A + B), -0.5*(min_y + max_y));

	x1 = a;
	if (pos == 0) y1 = formula(num, x1);
	if (pos == 1) y1 = v.method_compute(x1);
	if (pos == 2) y1 = w.method_compute(x1);
	if (pos == 3) { y1 = v.method_compute(x1); z1 = w.method_compute(x1);}
	if (pos == 4) { y1 = v.discrepancy(x1); z1 = w.discrepancy(x1);}
	
	for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
		if (pos == 0) y2 = formula(num, x2);
		if (pos == 1) y2 = v.method_compute(x2);
		if (pos == 2) {
			y2 = w.method_compute(x2);
		}
		if (pos == 3) {	y2 = v.method_compute(x2); z2 = w.method_compute(x2); }
		if (pos == 4) { y2 = v.discrepancy(x2); z2 = w.discrepancy(x2); }
		
		painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
		if (pos >= 3) {
			painter.drawLine (QPointF (x1, z1), QPointF (x2, z2));
			z1 = z2;
		}
		x1 = x2, y1 = y2; 
	}
		x2 = b;
		if (pos == 0) y2 = formula(num, x2);
		if (pos == 1) y2 = v.method_compute(x2);
		if (pos == 2) y2 = w.method_compute(x2);
		if (pos == 3) { y2 = v.method_compute(x2); z2 = w.method_compute(x2); }
		if (pos == 4) { y2 = v.discrepancy(x2); z2 = w.discrepancy(x2); }
	
		painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
		if (pos >= 3) painter.drawLine (QPointF (x1, z1), QPointF (x2, z2));

	painter.setPen ("blue");
	painter.drawLine (a, 0, b, 0);
	painter.drawLine (0, max_y, 0, min_y);

	painter.restore ();
}
