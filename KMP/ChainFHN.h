#ifndef CHAINFHN_H
#define CHAINFHN_H
#include <vector>
#include <math.h>
#include "SystemFHN.h"
#include<iostream>
using namespace std;
class ChainFHN {
	int n;
	float t0, tn, h;
	int m; //количество элементов в цепочке
public:
	vector<SystemFHN> sys;
	vector<float> d;
	ChainFHN() {
		t0 = 0; tn = 100;
		h = 0.1;
		n = (tn - t0) /h;
		m = 100;
		sys.resize(m);
		d.resize(m);
		for (int i = 0; i < m; i++) {
			sys[i].set(h, t0, tn);
			d[i] = 1;
		}
	}
	
	void set(float _h, float _t0, float _tn,int _m) {
		if (_h != h&& t0 != _t0 && _tn != tn && _m != m) {
			t0 = _t0;
			tn = _tn;
			h = _h;
			n = (tn - t0) / h;
			m = _m;
			sys.resize(m);
			d.resize(m);
			for (int i = 0; i < m; i++) {
				sys[i].set(h, t0, tn);
				d[i] = 1;
			}
		}
	}

	void setInitialCondition(double _x, double _y) {
		for (int i = 0; i < m; i++) {
			sys[i].setInitialCondition(_x, _y);
		}
	}

	void setParametrs(float _a, float _eps) {
		for (int i = 0; i < m; i++) {
			sys[i].setParametrs(_a, _eps);
		}
	}

	void setD(float _d, int i) {
		if (i < m)
			d[i] = _d;
	}

	double fStartEnd(double xi1, double xi, double yi, int i) {
		double temp = xi - pow(xi, 3) / 3 - yi;
		temp = temp + d[i] * (xi1 - xi);
		return temp;
	}
	double fMid(double xi1, double xi, double xi2, double yi, int i) {
		double temp = xi - pow(xi, 3) / 3 - yi;
		temp = temp + d[i] * (xi2 + xi1 - 2 * xi);
		return temp;
	}
	void methodRungeKutta(int j) {
		double k1, k2, k3, k4;
		double l1, l2, l3, l4;
		k1 = fStartEnd(sys[2].x[j], sys[1].x[j], sys[1].y[j], 1);
		l1 = sys[1].g(sys[1].x[j]);
		k2 = fStartEnd(sys[2].x[j] + h * k1 / 2, sys[1].x[j] + h * k1 / 2, sys[1].y[j] + h * l1 / 2, 1);
		l2 = sys[1].g(sys[1].x[j] + h * k1 / 2);
		k3 = fStartEnd(sys[2].x[j] + h * k2 / 2, sys[1].x[j] + h * k2 / 2, sys[1].y[j] + h * l2 / 2, 1);
		l3 = sys[1].g(sys[1].x[j] + h * k2 / 2);
		k4 = fStartEnd(sys[2].x[j] + h * k3, sys[1].x[j] + h * k3, sys[1].y[j] + h * l3 / 2, 1);
		l4 = sys[1].g(sys[1].x[j] + h * k3 / 2);
		sys[1].x[j + 1] = sys[1].x[j] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		sys[1].y[j + 1] = sys[1].y[j] + h * (l1 + 2 * l2 + 2 * l3 + l4) / 6;
		sys[0].x[j + 1] = sys[1].x[j + 1];
		sys[0].y[j + 1] = sys[1].y[j + 1];
		for (int i = 1; i < m - 1; i++) {
			k1 = fMid(sys[i - 1].x[j], sys[i].x[j], sys[i + 1].x[j], sys[i].y[j], i);
			l1 = sys[i].g(sys[i].x[j]);
			k2 = fMid(sys[i - 1].x[j] + h * k1 / 2, sys[i].x[j] + h * k1 / 2, sys[i + 1].x[j] + h * k1 / 2, sys[i].y[j] + h * l1 / 2, i);
			l2 = sys[i].g(sys[i].x[j] + h * k1 / 2);
			k3 = fMid(sys[i - 1].x[j] + h * k2 / 2, sys[i].x[j] + h * k2 / 2., sys[i + 1].x[j] + h * k2 / 2, sys[i].y[j] + h * l2 / 2., i);
			l3 = sys[i].g(sys[i].x[j] + h * k2 / 2);
			k4 = fMid(sys[i - 1].x[j] + h * k3, sys[i].x[j] + h * k3, sys[i + 1].x[j] + h * k3, sys[i].y[j] + h * l3, i);
			l4 = sys[i].g(sys[i].x[j] + h * k3);
			sys[i].x[j + 1] = sys[i].x[j] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
			sys[i].y[j + 1] = sys[i].y[j] + h * (l1 + 2 * l2 + 2 * l3 + l4) / 6.;
		}
		k1 = fStartEnd(sys[m - 3].x[j], sys[m - 2].x[j], sys[m - 2].y[j], m - 2);
		l1 = sys[m - 2].g(sys[m - 2].x[j]);
		k2 = fStartEnd(sys[m - 3].x[j] + h * k1 / 2., sys[m - 2].x[j] + h * k1 / 2., sys[m - 2].y[j] + h * l1 / 2., m - 2);
		l2 = sys[m - 2].g(sys[m - 2].x[j] + h * k1 / 2);
		k3 = fStartEnd(sys[m - 3].x[j] + h * k2 / 2., sys[m - 2].x[j] + h * k2 / 2., sys[m - 2].y[j] + h * l2 / 2., m - 2);
		l3 = sys[m - 2].g(sys[m - 2].x[j] + h * k2 / 2);
		k4 = fStartEnd(sys[m - 3].x[j] + h * k3, sys[m - 2].x[j] + h * k3, sys[m - 2].y[j] + h * l3, m - 2);
		l4 = sys[m - 2].g(sys[m - 2].x[j] + h * k3);
		sys[m - 2].x[j + 1] = sys[m - 2].x[j] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
		sys[m - 2].y[j + 1] = sys[m - 2].y[j] + h * (l1 + 2 * l2 + 2 * l3 + l4) / 6.;

		sys[m - 1].x[j + 1] = sys[m - 2].x[j + 1];
		sys[m - 1].y[j + 1] = sys[m - 2].y[j + 1];
	}
	void Calculate() {
		for (int j = 0; j < n - 1; j++) {
			methodRungeKutta(j);
			
		}
		
	}
};
#endif // !