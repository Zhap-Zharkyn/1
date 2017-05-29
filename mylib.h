#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#define EPS 1e-16

using namespace std;

class first {
	int N;
	int n;
	int num;
	double *y;
	double *c;
	
	public: 
		double *x;
		double *f;
		first(int N, int n, int j) : N(N), n(n), num(j) {
			x = new double[N];
			y = new double[N];
			f = new double[N];
			c = new double[n];
			memset(c, 0, sizeof(double) * n);

		}
		
		~first(){
			delete[] x;
			delete[] y;
			delete[] f;
			delete[] c;
		}
			
		double discrepancy_res(int M);
		double discrepancy(double z);
		void method_init();
		double method_compute(double z);
		void change();
		double a_ij(int i, int j);
		double b_ij(int i, int j);
		double c_ij(int i, int j);
		double d_ij(int i, int j);
		double u_ij(int i, int j);
		void resize(int m);
};

double a(double y, int i);
double b(double y, int i);
double T(double y, int i);
double U(double y, int i);

class second {
	int n;
	int num;
	double *d;
	double *c; 
	double *A; 
	double *b; 
	
	public: 
		double *x;
		double *f;
		second(int n, int j) : n(n), num(j) {
			x = new double[n];
			f = new double[n];
			d = new double[n];
			c = new double[(n - 1) * 4];
			A = new double[n * n];
			b = new double[n];
			memset(A, 0, sizeof(double) * n * n);
		}
		
		~second(){
			delete[] x;
			delete[] f;
			delete[] d;
			delete[] c;
			delete[] A;
			delete[] b;
		}
		
		double discrepancy_res(int N);
		double discrepancy(double z);
		void method_init();
		double method_compute(double y);
		void fill_in();
		void method_Gauss();
		double Q(double z);
		double R(double z);
		void resize(int m);		
};

double formula(int j, double y);

int search(double y, double *x, int n);

double f1(double x);
double f2(double x);
double f3(double x);

