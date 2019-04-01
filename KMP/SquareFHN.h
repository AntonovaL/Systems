#ifndef SQUAREFHN_H
#define SQUAREFHN_H
#include <vector>
#include <math.h>
#include "SystemFHN.h"
#include "ChainFHN.h"
#include<iostream>
using namespace std; 

class SquareFHN {
	int m;
	int k;
	int n;
	float t0, tn, h;
public:
	vector<ChainFHN> sys;
	SquareFHN() {
		t0 = 0; tn = 100;
		h = 0.1;
		n = (tn - t0) / h;
		m = 100;
		k = 100;
		sys.resize(k);
		for (int i = 0; i < k; i++) {
			sys[i].set(h, t0, tn,m);
		}
	}
	void set(float _h, float _t0, float _tn, int _m, int _k) {
		if (_h != h&& t0 != _t0 && _tn != tn && _m != m && k!=_k) {
			t0 = _t0;
			tn = _tn;
			h = _h;
			n = (tn - t0) / h;
			m = _m;
			k = _k;
			sys.resize(k);
			for (int i = 0; i < k; i++) {
				sys[i].set(h, t0, tn, m);
			}
		}
	}
	void setInitialCondition(double _x, double _y) {
		for (int i = 0; i < k; i++) {
			for (int j=0; j<m; k++)
				sys[i].sys[j].setInitialCondition(_x, _y);
		}
	}
	void setParametrs(float _a, float _eps) {
		for (int i = 0; i < k; i++) {
			for(int j=0; j<m; j++)
				sys[i].sys[j].setParametrs(_a, _eps);
		}
	}
	void setD(float _d, int i, int j) {
		if (i < k && j < m)
			sys[i].d[j] = _d;
	}

	double d_corners(double x1, double x2, double x3, int i, int j) {
		double d = sys[i].d[j] * (x1 + x2 - 2 * x3);
		return d;
	}
	double d_edges(double x1, double x2, double x3, double x4, int i, int j) {
		double d = sys[i].d[j] * (x1 + x2 + x3 - 3 * x4);
		return d;
	}
	double d_inner(double x1, double x2, double x3, double x4, double x5, int i, int j) {
		double d = sys[i].d[j] * (x1 + x2 + x3 + x4 - 4 * x5);
		return d;
	}
	double x_corners(double x1, double x2, double x3, double y,int i, int j) {
		double res = sys[i].sys[j].f(x3, y);
		res += d_corners(x1, x2, x3, i, j);
		return res;
	}
	double x_edges(double x1, double x2, double x3, double x4, double y,int i, int j) {
		double res = sys[i].sys[j].f(x4, y);
		res += d_edges(x1, x2, x3, x4, i, j);
		return res;
	}
	double x_inner(double x1, double x2, double x3, double x4, double x5, double y, int i, int j) {
		double res;
		res = sys[i].sys[j].f(x5, y);
		res += d_inner(x1, x2, x3, x4, x5, i, j);
		return res;
	}

	void MethodRunge_Kutta_Corners(int i, int j, int z)
	{

		double k1, k2, k3, k4;
		double l1, l2, l3, l4;
		if (i == k - 1 && j == 0)
		{
			k1 = x_corners(sys[i - 1].sys[j].x[z], sys[i].sys[j + 1].x[z], sys[i].sys[j].x[z], sys[i].sys[j].y[z], i, j);
			l1 = sys[i].sys[j].g(sys[i].sys[j].x[z]);

			k2 = x_corners(sys[i - 1].sys[j].x[z] + h*k1 / 2, sys[i].sys[j + 1].x[z] + h*k1 / 2, sys[i].sys[j].x[z] + h*k1 / 2, sys[i].sys[j].y[z] + h*l1 / 2, i, j);
			l2 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k1 / 2);

			k3 = x_corners(sys[i - 1].sys[j].x[z] + h*k2 / 2, sys[i].sys[j + 1].x[z] + h*k2 / 2, sys[i].sys[j].x[z] + h*k2 / 2, sys[i].sys[j].y[z] + h*l2 / 2, i, j);
			l3 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k2 / 2);

			k4 = x_corners(sys[i - 1].sys[j].x[z] + h*k3, sys[i].sys[j + 1].x[z] + h*k3, sys[i].sys[j].x[z] + h*k3, sys[i].sys[j].y[z] + h*l3, i, j);
			l4 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k3);
		}
		if (i == 0 && j == m - 1)
		{
			k1 = x_corners(sys[i].sys[j - 1].x[z], sys[i + 1].sys[j].x[z], sys[i].sys[j].x[z], sys[i].sys[j].y[z], i, j);
			l1 = sys[i].sys[j].g(sys[i].sys[j].x[z]);

			k2 = x_corners(sys[i].sys[j - 1].x[z] + h*k1 / 2, sys[i + 1].sys[j].x[z] + h*k1 / 2, sys[i].sys[j].x[z] + h*k1 / 2, sys[i].sys[j].y[z] + h*l1 / 2, i, j);
			l2 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k1 / 2);

			k3 = x_corners(sys[i].sys[j - 1].x[z] + h*k2 / 2, sys[i + 1].sys[j].x[z] + h*k2 / 2, sys[i].sys[j].x[z] + h*k2 / 2, sys[i].sys[j].y[z] + h*l2 / 2, i, j);
			l3 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k2 / 2);

			k4 = x_corners(sys[i].sys[j - 1].x[z] + h*k3, sys[i + 1].sys[j].x[z] + h*k3, sys[i].sys[j].x[z] + h*k3, sys[i].sys[j].y[z] + h*l3, i, j);
			l4 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k3);
		}
		if (i == 0 && j == 0)
		{
			k1 = x_corners(sys[0].sys[1].x[z], sys[1].sys[0].x[z], sys[0].sys[0].x[z], sys[0].sys[0].y[z], 0, 0);
			l1 = sys[0].sys[0].g(sys[0].sys[0].x[z]);

			k2 = x_corners(sys[0].sys[1].x[z] + h*k1 / 2, sys[1].sys[0].x[z] + h*k1 / 2, sys[0].sys[0].x[z] + h*k1 / 2, sys[0].sys[0].y[z] + h*l1 / 2, 0, 0);
			l2 = sys[0].sys[0].g(sys[0].sys[0].x[z] + h*k1 / 2);

			k3 = x_corners(sys[0].sys[1].x[z] + h*k2 / 2, sys[1].sys[0].x[z] + h*k2 / 2, sys[0].sys[0].x[z] + h*k2 / 2, sys[0].sys[0].y[z] + h*l2 / 2, 0, 0);
			l3 = sys[0].sys[0].g(sys[0].sys[0].x[z] + h*k2 / 2);


			k4 = x_corners(sys[0].sys[1].x[z] + h*k3, sys[1].sys[0].x[z] + h*k3, sys[0].sys[0].x[z] + h*k3, sys[0].sys[0].y[z] + h*l3, 0, 0);
			l4 = sys[0].sys[0].g(sys[0].sys[0].x[z] + h*k3);
		}
		if (i == k - 1 && j == m - 1)
		{
			k1 = x_corners(sys[k - 2].sys[m - 1].x[z], sys[k - 1].sys[m - 2].x[z], sys[k - 1].sys[m - 1].x[z], sys[k - 1].sys[m - 1].y[z], k - 1, m - 1);
			l1 = sys[k - 1].sys[m - 1].g(sys[k - 1].sys[m - 1].x[z]);

			k2 = x_corners(sys[k - 2].sys[m - 1].x[z] + h*k1 / 2, sys[k - 1].sys[m - 2].x[z] + h*k1 / 2, sys[k - 1].sys[m - 1].x[z] + h*k1 / 2, sys[k - 1].sys[m - 1].y[z] + h*l1 / 2, k - 1, m - 1);
			l2 = sys[k - 1].sys[m - 1].g(sys[k - 1].sys[m - 1].x[z] + h*k1 / 2);

			k3 = x_corners(sys[k - 2].sys[m - 1].x[z] + h*k2 / 2, sys[k - 1].sys[m - 2].x[z] + h*k2 / 2, sys[k - 1].sys[m - 1].x[z] + h*k2 / 2, sys[k - 1].sys[m - 1].y[z] + h*l2 / 2, k - 1, m - 1);
			l3 = sys[k - 1].sys[m - 1].g(sys[k - 1].sys[m - 1].x[z] + h*k2 / 2);

			k4 = x_corners(sys[k - 2].sys[m - 1].x[z] + h*k3, sys[k - 1].sys[m - 2].x[z] + h*k3, sys[k - 1].sys[m - 1].x[z] + h*k3, sys[k - 1].sys[m - 1].y[z] + h*l3, k - 1, m - 1);
			l4 = sys[k - 1].sys[m - 1].g(sys[k - 1].sys[m - 1].x[z] + h*k3);
		}
		sys[i].sys[j].x[z + 1] = sys[i].sys[j].x[z] + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6.;
		sys[i].sys[j].y[z + 1] = sys[i].sys[j].y[z] + h*(l1 + 2 * l2 + 2 * l3 + l4) / 6.;
	}

	void MethodRunge_Kutta_Edges(int i, int j, int z)
	{
		double k1, k2, k3, k4;
		double l1, l2, l3, l4;
		if (i == 0)
		{
			k1 = x_edges(sys[i].sys[j - 1].x[z], sys[i].sys[j + 1].x[z], sys[i + 1].sys[j].x[z], sys[i].sys[j].x[z], sys[i].sys[j].y[z], i, j);
			l1 = sys[i].sys[j].g(sys[i].sys[j].x[z]);

			k2 = x_edges(sys[i].sys[j - 1].x[z] + h*k1 / 2, sys[i].sys[j + 1].x[z] + h*k1 / 2, sys[i + 1].sys[j].x[z] + h*k1 / 2, sys[i].sys[j].x[z] + h*k1 / 2, sys[i].sys[j].y[z] + h*l1 / 2, i, j);
			l2 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k1 / 2);

			k3 = x_edges(sys[i].sys[j - 1].x[z] + h*k2 / 2, sys[i].sys[j + 1].x[z] + h*k2 / 2, sys[i + 1].sys[j].x[z] + h*k2 / 2, sys[i].sys[j].x[z] + h*k2 / 2, sys[i].sys[j].y[z] + h*l2 / 2, i, j);
			l3 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k2 / 2);

			k4 = x_edges(sys[i].sys[j - 1].x[z] + h*k3, sys[i].sys[j + 1].x[z] + h*k3, sys[i + 1].sys[j].x[z] + h*k3, sys[i].sys[j].x[z] + h*k3, sys[i].sys[j].y[z] + h*l3, i, j);
			l4 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k3);
		}
		if (j == 0)
		{
			k1 = x_edges(sys[i - 1].sys[j].x[z], sys[i].sys[j + 1].x[z], sys[i + 1].sys[j].x[z], sys[i].sys[j].x[z], sys[i].sys[j].y[z], i, j);
			l1 = sys[i].sys[j].g(sys[i].sys[j].x[z]);

			k2 = x_edges(sys[i - 1].sys[j].x[z] + h*k1 / 2, sys[i].sys[j + 1].x[z] + h*k1 / 2, sys[i + 1].sys[j].x[z] + h*k1 / 2, sys[i].sys[j].x[z] + h*k1 / 2, sys[i].sys[j].y[z] + h*l1 / 2, i, j);
			l2 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k1 / 2);

			k3 = x_edges(sys[i - 1].sys[j].x[z] + h*k2 / 2, sys[i].sys[j + 1].x[z] + h*k2 / 2, sys[i + 1].sys[j].x[z] + h*k2 / 2, sys[i].sys[j].x[z] + h*k2 / 2, sys[i].sys[j].y[z] + h*l2 / 2, i, j);
			l3 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k2 / 2);

			k4 = x_edges(sys[i - 1].sys[j].x[z] + h*k3, sys[i].sys[j + 1].x[z] + h*k3, sys[i + 1].sys[j].x[z] + h*k3, sys[i].sys[j].x[z] + h*k3, sys[i].sys[j].y[z] + h*l3, i, j);
			l4 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k3);
		}
		if (i == k - 1)
		{
			k1 = x_edges(sys[i].sys[j - 1].x[z], sys[i].sys[j + 1].x[z], sys[i - 1].sys[j].x[z], sys[i].sys[j].x[z], sys[i].sys[j].y[z], i, j);
			l1 = sys[i].sys[j].g(sys[i].sys[j].x[z]);

			k2 = x_edges(sys[i].sys[j - 1].x[z] + h*k1 / 2, sys[i].sys[j + 1].x[z] + h*k1 / 2, sys[i - 1].sys[j].x[z] + h*k1 / 2, sys[i].sys[j].x[z] + h*k1 / 2, sys[i].sys[j].y[z] + h*l1 / 2, i, j);
			l2 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k1 / 2);

			k3 = x_edges(sys[i].sys[j - 1].x[z] + h*k2 / 2, sys[i].sys[j + 1].x[z] + h*k2 / 2, sys[i - 1].sys[j].x[z] + h*k2 / 2, sys[i].sys[j].x[z] + h*k2 / 2, sys[i].sys[j].y[z] + h*l2 / 2, i, j);
			l3 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k2 / 2);

			k4 = x_edges(sys[i].sys[j - 1].x[z] + h*k3, sys[i].sys[j + 1].x[z] + h*k3, sys[i - 1].sys[j].x[z] + h*k3, sys[i].sys[j].x[z] + h*k3, sys[i].sys[j].y[z] + h*l3, i, j);
			l4 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k3);
		}
		if (j == m - 1)
		{
			k1 = x_edges(sys[i - 1].sys[j].x[z], sys[i].sys[j - 1].x[z], sys[i + 1].sys[j].x[z], sys[i].sys[j].x[z], sys[i].sys[j].y[z], i, j);
			l1 = sys[i].sys[j].g(sys[i].sys[j].x[z]);

			k2 = x_edges(sys[i - 1].sys[j].x[z] + h*k1 / 2, sys[i].sys[j - 1].x[z] + h*k1 / 2, sys[i + 1].sys[j].x[z] + h*k1 / 2, sys[i].sys[j].x[z] + h*k1 / 2, sys[i].sys[j].y[z] + h*l1 / 2, i, j);
			l2 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k1 / 2);

			k3 = x_edges(sys[i - 1].sys[j].x[z] + h*k2 / 2, sys[i].sys[j - 1].x[z] + h*k2 / 2, sys[i + 1].sys[j].x[z] + h*k2 / 2, sys[i].sys[j].x[z] + h*k2 / 2, sys[i].sys[j].y[z] + h*l2 / 2, i, j);
			l3 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k2 / 2);

			k4 = x_edges(sys[i - 1].sys[j].x[z] + h*k3, sys[i].sys[j - 1].x[z] + h*k3, sys[i + 1].sys[j].x[z] + h*k3, sys[i].sys[j].x[z] + h*k3, sys[i].sys[j].y[z] + h*l3, i, j);
			l4 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k3);
		}

		sys[i].sys[j].x[z + 1] = sys[i].sys[j].x[z] + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6.;
		sys[i].sys[j].y[z + 1] = sys[i].sys[j].y[z] + h*(l1 + 2 * l2 + 2 * l3 + l4) / 6.;
	}


	void MethodRunge_Kutta_Inner(int i, int j, int z)
	{
		double k1, k2, k3, k4;
		double l1, l2, l3, l4;

		k1 = x_inner(sys[i].sys[j - 1].x[z], sys[i + 1].sys[j].x[z], sys[i - 1].sys[j].x[z], sys[i].sys[j + 1].x[z], sys[i].sys[j].x[z], sys[i].sys[j].y[z], i, j);
		l1 = sys[i].sys[j].g(sys[i].sys[j].x[z]);

		k2 = x_inner(sys[i].sys[j - 1].x[z] + h*k1 / 2, sys[i + 1].sys[j].x[z] + h*k1 / 2, sys[i - 1].sys[j].x[z] + h*k1 / 2, sys[i].sys[j + 1].x[z] + h*k1 / 2, sys[i].sys[j].x[z] + h*k1 / 2, sys[i].sys[j].y[z] + h*l1 / 2, i, j);
		l2 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k1 / 2);

		k3 = x_inner(sys[i].sys[j - 1].x[z] + h*k2 / 2, sys[i + 1].sys[j].x[z] + h*k2 / 2, sys[i - 1].sys[j].x[z] + h*k2 / 2, sys[i].sys[j + 1].x[z] + h*k2 / 2, sys[i].sys[j].x[z] + h*k2 / 2, sys[i].sys[j].y[z] + h*l2 / 2, i, j);
		l3 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k2 / 2);

		k4 = x_inner(sys[i].sys[j - 1].x[z] + h*k3, sys[i + 1].sys[j].x[z] + h*k3, sys[i - 1].sys[j].x[z] + h*k3, sys[i].sys[j + 1].x[z] + h*k3, sys[i].sys[j].x[z] + h*k3, sys[i].sys[j].y[z] + h*l3, i, j);
		l4 = sys[i].sys[j].g(sys[i].sys[j].x[z] + h*k3);
		sys[i].sys[j].x[z + 1] = sys[i].sys[j].x[z] + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6.;
		sys[i].sys[j].y[z + 1] = sys[i].sys[j].y[z] + h*(l1 + 2 * l2 + 2 * l3 + l4) / 6.;
	}

	void MethodRunge_Kutta(int z)//вычисление методом Цунге-†утта 4-го порІдка
	{

		MethodRunge_Kutta_Corners(0, 0, z);
		MethodRunge_Kutta_Corners(0, m - 1, z);
		MethodRunge_Kutta_Corners(k - 1, 0, z);
		MethodRunge_Kutta_Corners(k - 1, m - 1, z);
		for (int i = 1; i < k-1; i++)
			for (int j = 1; j < m-1; j++)
				MethodRunge_Kutta_Inner(i, j, z);
		for (int j = 1; j < m-1; j++)
			MethodRunge_Kutta_Edges(0, j, z);
		for (int j = 1; j < k-1; j++)
			MethodRunge_Kutta_Edges(j, 0, z);
		for (int j = 1; j < m-1; j++)
			MethodRunge_Kutta_Edges(k - 1, j, z);
		for (int j = 1; j < k-1; j++)
			MethodRunge_Kutta_Edges(j, m - 1, z);
		


	}

	void MethodRunge_Kutta()//вычисление методом ÷унге-Жутта 4-го порІдка
	{
		for (int z = 0; z<n - 1; z++)
		{
			MethodRunge_Kutta(z);
		}

	}
	void CalculationFunctions()//решение системы
	{
		MethodRunge_Kutta();
	}

};
#endif
