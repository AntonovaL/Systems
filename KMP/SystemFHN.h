#ifndef SYSTEMFHN_H
#define SYSTEMFHN_H


#include <vector>
#include <math.h>

using namespace std;

class SystemFHN {
	float a, eps;
	int n;
	float h, t0, tn;
public:
	vector<double> x;
	vector<double> y;
	SystemFHN() {
		a = 1.1;
		eps = 0.05;
		t0 = 0; tn = 100;
		h = 0.1;
		n = (tn - t0) / h;
		x.resize(n);
		y.resize(n);
		for (int i = 0; i < n; i++) {
			x[i] = 0;
			y[i] = 0;
		}
	}

	void setParametrs(float _a, float _eps) {
		a = _a;
		eps = _eps;
	}

	void set(float _h, float _t0, float _tn) {
		if (_h != h && t0 != _t0 && tn != _tn) {
			h = _h;
			t0 = _t0;
			tn = _tn;
			n = (tn - t0) / h;
			x.resize(n);
			y.resize(n);
			for (int i = 0; i < n; i++) {
				x[i] = 0;
				y[i] = 0;
			}
		}
	}

	void setInitialCondition(double _x0, double _y0) {
		x[0] = _x0;
		y[0] = _y0;
	}

	int getN() {
		return n;
	}
	double f(double x,double y) {
		return x - pow(x, 3)/3 - y;
	}

	double g(double x) {
		return eps*(x - a);
	}
	void methodRungeKutta() {
		double k1, k2, k3, k4;
		double l1, l2, l3, l4;
		for (int i = 0; i < n-1; i++) {
			k1 = f(x[i], y[i]);
			l1 = g(x[i]);
			k2 = f(x[i] + h*k1 / 2, y[i] + h*l1 / 2);
			l2 = g(x[i] + h*k1 / 2);

			k3 = f(x[i] + h*k2 / 2, y[i] + h*l2 / 2);
			l3 = g(x[i] + h*k2 / 2);

			k4 = f(x[i] + h*k3, y[i] + h*l3);
			l4 = g(x[i] + h*k3);

			x[i + 1] = x[i] + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
			y[i + 1] = y[i] + h*(l1 + 2 * l2 + 2 * l3 + l4) / 6;
		}
	}
};

#endif // !