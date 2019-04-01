#ifndef CHAINAP_H
#define CHAINAP_H
#include <vector>
#include <math.h>
#include "SystemAP.h"
#include<iostream>
using namespace std;
class ChainAP {
	int n;
	float t0, tn, h;
	int m; //количество элементов в цепочке
public:
	vector<SystemAP> sys;
	vector<float> d;
	ChainAP() {
		t0 = 0; tn = 100;
		h = 0.01;
		n = (tn - t0) / h;
		m = 100;
		sys.resize(m);
		d.resize(m);
		for (int i = 0; i < m; i++) {
			sys[i].set(h, t0, tn);
			d[i] = 1;
		}
	}

	void set(float _h, float _t0, float _tn, int _m) {
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

	void setInitialCondition( double _u,  double _v, double _ta) {
		for (int i = 0; i < m; i++) {
			sys[i].setInitialCondition(_u, _v,_ta);
		}
	}

	void setParametrs(float _a, float _k, float _kt, float _eps1, float _eps0) {
		for (int i = 0; i < m; i++) {
			sys[i].setParametrs(_a, _k, _kt, _eps1, _eps0);
		}
	}

	void setD(float _d, int i) {
		if (i < m)
			d[i] = _d;
	}

	 double fStartEnd( double ui1,  double ui,  double vi, int i) {
		double temp = sys[i].f(ui,vi,i);
		temp = temp + d[i] * (ui1 - ui);
		return temp;
	}
	 double fMid( double ui1,  double ui, double ui2,  double vi, int i) {
		 double temp = sys[i].f(ui, vi, i);
		temp = temp + d[i] * (ui2 + ui1 - 2 * ui);
		return temp;
	}

	void methodRungeKutta(int j) {
		 double k1, k2, k3, k4;
		 double l1, l2, l3, l4;
		 double m1, m2, m3, m4;
		k1 = fStartEnd(sys[2].u[j], sys[1].u[j], sys[1].v[j], 1);
		l1 = sys[1].g(sys[1].u[j], sys[1].v[j]);
		m1 = sys[1].th(sys[1].u[j], sys[1].ta[j]);
		k2 = fStartEnd(sys[2].u[j] + h * k1 / 2, sys[1].u[j] + h * k1 / 2, sys[1].v[j] + h * l1 / 2, 1);
		l2 = sys[1].g(sys[1].u[j] + h * k1 / 2, sys[1].v[j] + h * l1 / 2);
		m2 = sys[1].th(sys[1].u[j] + h*k1, sys[1].ta[j] + h*m1);
		k3 = fStartEnd(sys[2].u[j] + h * k2 / 2, sys[1].u[j] + h * k2 / 2, sys[1].v[j] + h * l2 / 2, 1);
		l3 = sys[1].g(sys[1].u[j] + h * k2 / 2, sys[1].v[j] + h * l2 / 2);
		m3 = sys[1].th(sys[1].u[j] + h * k2, sys[1].ta[j] + h * m2);
		k4 = fStartEnd(sys[2].u[j] + h * k3, sys[1].u[j] + h * k3, sys[1].v[j] + h * l3 / 2, 1);
		l4 = sys[1].g(sys[1].u[j] + h * k3 / 2, sys[1].v[j] + h * l3 / 2);
		m4 = sys[1].th(sys[1].u[j] + h * k3, sys[1].ta[j] + h * m3);
		sys[1].u[j + 1] = sys[1].u[j] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
		sys[1].v[j + 1] = sys[1].v[j] + h * (l1 + 2 * l2 + 2 * l3 + l4) / 6.;
		sys[1].ta[j + 1] = sys[1].ta[j] + h * (m1 + 2 * m2 + 2 * m3 + m4) / 6.;
		sys[0].u[j + 1] = sys[1].u[j + 1];
		sys[0].v[j + 1] = sys[1].v[j + 1];
		sys[0].ta[j + 1] = sys[1].ta[j + 1];
		for (int i = 1; i < m - 1; i++) {
			k1 = fMid(sys[i - 1].u[j], sys[i].u[j], sys[i + 1].u[j], sys[i].v[j], i);
			l1 = sys[i].g(sys[i].u[j], sys[i].v[j]);
			m1 = sys[i].th(sys[i].u[j], sys[i].ta[j]);
			k2 = fMid(sys[i - 1].u[j] + h * k1 / 2, sys[i].u[j] + h * k1 / 2, sys[i + 1].u[j] + h * k1 / 2, sys[i].v[j] + h * l1 / 2, i);
			l2 = sys[i].g(sys[i].u[j] + h * k1 / 2, sys[i].v[j] + h * l1 / 2);
			m2 = sys[i].th(sys[i].u[j] + h*k1 / 2, sys[i].ta[j] + h*m1 / 2);
			k3 = fMid(sys[i - 1].u[j] + h * k2 / 2, sys[i].u[j] + h * k2 / 2, sys[i + 1].u[j] + h * k2 / 2, sys[i].v[j] + h * l2 / 2, i);
			l3 = sys[i].g(sys[i].u[j] + h * k2 / 2, sys[i].v[j] + h * l2 / 2);
			m3 = sys[i].th(sys[i].u[j] + h * k2 / 2, sys[i].ta[j] + h * m2 / 2);
			k4 = fMid(sys[i - 1].u[j] + h * k3, sys[i].u[j] + h * k3, sys[i + 1].u[j] + h * k3, sys[i].v[j] + h * l3, i);
			l4 = sys[i].g(sys[i].u[j] + h * k3, sys[i].v[j] + h * l3);
			m4 = sys[i].th(sys[i].u[j] + h * k3 / 2, sys[i].ta[j] + h * m3 / 2);
			sys[i].u[j + 1] = sys[i].u[j] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
			sys[i].v[j + 1] = sys[i].v[j] + h * (l1 + 2 * l2 + 2 * l3 + l4) / 6.;
			sys[i].ta[j + 1] = sys[i].ta[j] + h * (m1 + 2 * m2 + 2 * m3 + m4) / 6.;
		}
		k1 = fStartEnd(sys[m - 3].u[j], sys[m - 2].u[j], sys[m - 2].v[j], m - 2);
		l1 = sys[m - 2].g(sys[m - 2].u[j], sys[m - 2].v[j]);
		m1 = sys[m - 2].th(sys[m - 1].u[j], sys[m - 1].ta[j]);

		k2 = fStartEnd(sys[m - 3].u[j] + h * k1 / 2, sys[m - 2].u[j] + h * k1 / 2, sys[m - 2].v[j] + h * l1 / 2, m - 2);
		l2 = sys[m - 2].g(sys[m - 2].u[j] + h * k1 / 2, sys[m - 2].v[j] + h*l1 / 2);
		m2 = sys[m - 2].th(sys[m - 1].u[j] + h * k1 / 2, sys[m - 1].ta[j] + h*m1 / 2);
		k3 = fStartEnd(sys[m - 3].u[j] + h * k2 / 2, sys[m - 2].u[j] + h * k2 / 2, sys[m - 2].v[j] + h * l2 / 2, m - 2);
		l3 = sys[m - 2].g(sys[m - 2].u[j] + h * k2 / 2, sys[m - 2].v[j] + h * l2 / 2);
		m3 = sys[m - 2].th(sys[m - 1].u[j] + h * k2 / 2, sys[m - 1].ta[j] + h * m2 / 2);
		k4 = fStartEnd(sys[m - 3].u[j] + h * k3, sys[m - 2].u[j] + h * k3, sys[m - 2].v[j] + h * l3, m - 2);
		l4 = sys[m - 2].g(sys[m - 2].u[j] + h * k3, sys[m - 2].v[j] + h * l3);
		m4 = sys[m - 2].th(sys[m - 1].u[j] + h * k3 / 2, sys[m - 1].ta[j] + h * m3 / 2);
		sys[m - 2].u[j + 1] = sys[m - 2].u[j] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
		sys[m - 2].v[j + 1] = sys[m - 2].v[j] + h * (l1 + 2 * l2 + 2 * l3 + l4) / 6.;
		sys[m - 2].ta[j + 1] = sys[m - 2].ta[j] + h * (m1 + 2 * m2 + 2 * m3 + m4) / 6.;
		sys[m - 1].u[j + 1] = sys[m - 2].u[j + 1];
		sys[m - 1].v[j + 1] = sys[m - 2].v[j + 1];
		sys[m - 1].ta[j + 1] = sys[m - 2].ta[j + 1];
	}
	void Calculate() {
		for (int j = 0; j < n - 1; j++) {
			methodRungeKutta(j);

		}

	}
};
#endif // !
