#include "mylib.h"

int search(double y, double *x, int n) {
	int left = 0, right = n, midd = 0;
    if (fabs(x[n - 1] -y) < EPS) return n - 2;
    if (fabs(x[0] - y) < EPS) return 0;
	while (true) { 
		midd = (left + right) / 2;
		if (y >= x[midd] && y < x[midd + 1]) 
			return midd;   
		else if (y < x[midd])  
			right = midd;      
			else if (y > x[midd])  
				left = midd;	   
	}
}

double formula(int j, double x) { 
	switch (j) {
		case 1 : 
			return f1(x);
		case 2 :
			return f2(x);
		default : 
			cout<<"ERROR"<<endl;
			return 0.0;
	}
}

double f_formula(double x, double x1, double x2, double y1, double y2) {
	return ((x1 * y2 - x2 * y1) - (y1 - y2) * x) / (x2 - x1);
}

double f1(double /*x*/) {
	return 1;
}

double f2(double x) { 
	return x * x;
}

////////////////////////////////////////////////////////////////////////

double first::discrepancy_res(int M) {
	int k = 1, i;
	double s, h, max = 0.0;
	h = (x[N - 1] - x[0]) / M; 
	if (num == 0) {
		max = fabs(method_compute(x[0]) - f_formula(x[0], x[0], x[1], f[0], f[1]));
		for (i = 1; i < M; i++) {
			if(x[0] + i * h > x[i]) k++;
			s = fabs(method_compute(x[0] + i * h) - f_formula(x[0] + i * h, x[k - 1], x[k], f[k - 1], f[k]));
			if (s > max) max = s;
		}
	}
	else {
		max = fabs(method_compute(x[0]) - formula(num, x[0]));
		for (i = 1; i < M; i++) {
			s = fabs(method_compute(x[0] + i * h) - formula(num, x[0] + i * h));
			if (s > max) max = s;
		}
	}
	return max;
} 

double first::discrepancy(double z) {
	int k;
	if (num == 0) {
		k = search(z, x, N);
		return method_compute(z) - f_formula(z, x[k], x[k + 1], f[k], f[k + 1]);
	}
	return method_compute(z) - formula(num, z);
} 

void first::method_init() { //n количество точек
	int i, j;
	change();
	for (j = 0; j < N; j++)
		for (i = 0; i < n; i++) 
			c[i] += u_ij(i, j) * f[j];
	c[0] /=  M_PI;
	for (i = 1; i < n; i++) 
		c[i] *= 2 / M_PI;
}

double first::method_compute(double z) { //n - 1 
	double res = 0.0;
	for(int i = 0; i < n; i++) 
		res += c[i] * T((2 * z - (x[N - 1] + x[0])) / (x[N - 1] - x[0]), i);
	return res ; 
}

void first::change(){
	y[0] = -1;
	y[N - 1] = 1;
	for (int i = 1; i < N - 1; i++) 
		y[i] = (2 * x[i] - (x[N - 1] + x[0])) / (x[N - 1] - x[0]);
}

double a(double z, int i) {
	if (i == 0) return acos(z);
	else return U(z, i - 1) / i;
}

double b(double z, int i) {
	if (i == 0) return U(z, i);
	else if (i == 1) return (acos(z) + U(z, i) / 2) / 2;
		else return (U(z, i - 2) / (i - 1) + U(z, i) / (i + 1)) / 2;  
} 

double T(double z, int i) { 
	if (z > 1) z = 1;
	if (z < -1) z = -1;
	return cos(i * acos(z));
}

double U(double z, int i) {
	return sin((i + 1) * acos(z));
}

double first::a_ij(int i, int j) {
	return -(a(y[j + 1], i) - a(y[j], i));
}

double first::b_ij(int i, int j) {
	return -(b(y[j + 1], i) - b(y[j], i));
}

double first::c_ij(int i, int j) {
	return (a_ij(i, j) * y[j + 1] - b_ij(i, j)) / (y[j + 1] - y[j]);
}

double first::d_ij(int i, int j) {
	return (b_ij(i, j) - a_ij(i, j) * y[j]) / (y[j + 1] - y[j]);
}


double first::u_ij(int i, int j) {
	 if(j == 0) return c_ij(i, j);
	 else if(j == N - 1) return d_ij(i, j - 1);
		  else return c_ij(i, j) + d_ij(i, j - 1);
}

////////////////////////////////////////////////////////////////////////

double second::discrepancy_res(int N) {
	int k = 1;
	double s, h, max = 0.0;
	h = (x[n - 1] - x[0]) / N; 
	if (num == 0) {
		max = fabs(method_compute(x[0]) - f_formula(x[0], x[0], x[1], f[0], f[1]));
		for (int i = 1; i < N; i++) {
			if(x[0] + i * h > x[i]) k++;
			s = fabs(method_compute(x[0] + i * h) - f_formula(x[0] + i * h, x[k - 1], x[k], f[k - 1], f[k]));
			if (s > max) max = s;
		}
	}
	else {
		max = fabs(method_compute(x[0]) - formula(num, x[0]));
		for (int i = 1; i < N; i++) {
			s = fabs(method_compute(x[0] + i * h) - formula(num, x[0] + i * h));
			if (s > max) max = s;
		}
	}
	return max;
} 

double second::discrepancy(double z) {
	int k = 1;
	double s;
	if (num == 0) {
		k = search(z, x, n);
		s = method_compute(z) - f_formula(z, x[k], x[k + 1], f[k], f[k + 1]);
	}
	s = method_compute(z) - formula(num, z);
	return s;
} 

void second::method_init() {
	double v, w;
	fill_in();
	method_Gauss();
	for (int i = 0; i < n - 1; i++) { 
		v = x[i + 1] - x[i];
		w = f[i + 1] - f[i];
		c[i * 4] = f[i];
		c[i * 4 + 1] = d[i];
		c[i * 4 + 2] = (3 * w / v - 2 * d[i] - d[i + 1]) / v;
		c[i * 4 + 3] = (d[i] + d[i + 1] - 2 * w / v) / (v * v);
	}
		
}

double second::method_compute(double y) { 
	int i;
	double res, s;
	if (y > x[n - 1]) y = x[n - 1];
	if (y < x[0]) y = x[0];
	i = search(y, x, n);
	 if (fabs(x[n - 1] -y) < EPS) i = n - 2;
	s = y - x[i];
	res = c[i * 4] + c[i * 4 + 1] * s + c[ i * 4 + 2] * s * s + c[i * 4 + 3] * s * s * s;
	return res;
}

void second::fill_in() {
	double v, w;
	A[0] = 1;
	A[(n - 1) * n + (n - 1)] = 1;
	b[0] = Q(x[0]);
	b[n - 1] = R(x[n - 1]);
	for (int i = 1; i < n - 1; i++) {
		v = (x[i] - x[i - 1]);
		w = (x[i + 1] - x[i]);
		A[i * n + (i - 1)] = x[i + 1] - x[i];
		A[i * n + i] = 2 * (x[i + 1] - x[i - 1]);
		A[i * n + (i + 1)] = x[i] - x[i - 1];
		b[i] = 3 * (f[i] - f[i - 1]) * w / v + 3 * (f[i + 1] - f[i]) * v / w;
	}
}
			
void second::method_Gauss() { 
	int k, i, j;
	double v;
	for (k = 0; k < n; k++) 
		for (i = k; i < n; i++) {
			v = b[k] / A[k * n + k];
			if (i == k) b[i] = v;
			else b[i] -= v * A[k * n + i];
			v = A[i * n + k ] / A[k * n + k];
			for (j =  k; j < n; j++) 
				if (i == k) A[i * n + j] /= A[i * n + i];
				else A[i * n + j] -=  v * A[k * n + i];
		}
	for (i = n - 1; i >= 0; i--) {
		v = 0.0;
        for (j = n - 1; j > i; j--)
            v += d[j] * A[i * n + j];
        if (fabs(A[i * n + i]) > EPS)
        d[i] = (b[i] - v) / A[i * n + i];
        else d[i] = 0.0;
    }
}

double second::Q(double z) {
	double res = 0.0;
	res += f[1] * (2 * z - x[1] - x[2]) / ((x[0] - x[1]) * (x[0] - x[2]));
	res += f[2] * (2 * z - x[0] - x[2]) / ((x[1] - x[0]) * (x[1] - x[2]));
	res += f[3] * (2 * z - x[0] - x[1]) / ((x[2] - x[0]) * (x[2] - x[1]));
	return res;
}

double second::R(double z) {
	double res = 0.0;
	res += f[n - 1] * (2 * z - x[n - 2] - x[n - 3]) / ((x[n - 1] - x[n - 2]) * (x[n - 1] - x[n - 3]));
	res += f[n - 2] * (2 * z - x[n - 1] - x[n - 3]) / ((x[n - 2] - x[n - 1]) * (x[n - 2] - x[n - 3]));
	res += f[n - 3] * (2 * z - x[n - 1] - x[n - 2]) / ((x[n - 3] - x[n - 1]) * (x[n - 3] - x[n - 2]));
	return res;
}
////////////////////////////////////////////////////////////////////////


void first::resize(int m) {
	delete[] x;
	delete[] y;
	delete[] f;
	delete[] c;
	N = m;
	c = new double[n];
	x = new double[N];
	y = new double[N];
	f = new double[N];
	memset(c, 0, sizeof(double) * n);
}

void second::resize(int m) {
	delete[] x;
	delete[] f;
	delete[] d;
	delete[] c;
	delete[] A;
	delete[] b;
	n = m;
	x = new double[n];
	f = new double[n];
	d = new double[n];
	c = new double[(n - 1) * 4];
	A = new double[n * n];
	b = new double[n];
	memset(A, 0, sizeof(double) * n * n);
}



