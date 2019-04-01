#ifndef SQUAREAP_H
#define SQUAREAP_H
#include <vector>
#include <math.h>
#include "SystemAP.h"
#include "SystemChainAP.h"
#include<iostream>
using namespace std;

class SquareAP {
	int m;
	int k;
	int n;
	float t0, tn, h;
public:
	vector<ChainAP> sys;
	SquareAP() {
		t0 = 0; tn = 100;
		h = 0.1;
		n = (tn - t0) / h;
		m = 100;
		k = 100;
		sys.resize(k);
		for (int i = 0; i < k; i++) {
			sys[i].set(h, t0, tn, m);
		}
	}
	void set(float _h, float _t0, float _tn, int _m, int _k) {
		if (_h != h&& t0 != _t0 && _tn != tn && _m != m && k != _k) {
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
	void setInitialCondition(long double _u, long double _v, double _ta) {
		for (int i = 0; i < k; i++) {
			for (int j = 0; j<m; k++)
				sys[i].sys[j].setInitialCondition(_u, _v,_ta);
		}
	}
	void setParametrs(float _a, float _k, float _kt, float _eps1, float _eps0) {
		for (int i = 0; i < k; i++) {
			for (int j = 0; j<m; j++)
				sys[i].sys[j].setParametrs(_a, _k, _kt, _eps1, _eps0);
		}
	}
	void setD(float _d, int i, int j) {
		if (i < k && j < m)
			sys[i].d[j] = _d;
	}

	 double d_corners( double u1,  double u2,  double u3, int i, int j) {
		 double d = sys[i].d[j] * (u1 + u2 - 2 * u3);
		return d;
	}
	 double d_edges( double u1,  double u2,  double u3,  double u4, int i, int j) {
		 double d = sys[i].d[j] * (u1 + u2 + u3 - 3 * u4);
		return d;
	}
	 double d_inner( double u1,  double u2,  double u3,  double u4,  double u5, int i, int j) {
		 double d = sys[i].d[j] * (u1 + u2 + u3 + u4 - 4 * u5);
		return d;
	}
	 double x_corners( double u1,  double u2,  double u3,  double v, int i, int j,int z) {
		 double res = sys[i].sys[j].f(u3, v,z);
		res += d_corners(u1, u2, u3, i, j);
		return res;
	}
	 double x_edges( double u1,  double u2,  double u3,  double u4,  double v, int i, int j,int z) {
		 double res = sys[i].sys[j].f(u4, v,z);
		res += d_edges(u1, u2, u3, u4, i, j);
		return res;
	}
	 double x_inner( double u1,  double u2,  double u3,  double u4,  double u5,  double v, int i, int j,int z) {
		 double res;
		res = sys[i].sys[j].f(u5, v, z);
		res += d_inner(u1, u2, u3, u4, u5, i, j);
		return res;
	}

	void MethodRunge_Kutta_Corners(int i, int j, int z)
	{

		 double k1, k2, k3, k4;
		 double l1, l2, l3, l4;
		 double m1, m2, m3, m4;
		if (i == k - 1 && j == 0)
		{
			k1 = x_corners(sys[i - 1].sys[j].u[z], sys[i].sys[j + 1].u[z], sys[i].sys[j].u[z], sys[i].sys[j].v[z], i, j,z);
			l1 = sys[i].sys[j].g(sys[i].sys[j].u[z], sys[i].sys[j].v[z]);
			m1 = sys[i].sys[j].th(sys[i].sys[j].u[z], sys[i].sys[j].ta[z]);
			
			k2 = x_corners(sys[i - 1].sys[j].u[z] + h*k1 / 2, sys[i].sys[j + 1].u[z] + h*k1 / 2, sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2, i, j,z);
			l2 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2);
			m2 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].ta[z] + h*m1 / 2);

			k3 = x_corners(sys[i - 1].sys[j].u[z] + h*k2 / 2, sys[i].sys[j + 1].u[z] + h*k2 / 2, sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l2 / 2, i, j,z);
			l3 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l2 / 2);
			m3 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].ta[z] + h*m2 / 2);

			k4 = x_corners(sys[i - 1].sys[j].u[z] + h*k3, sys[i].sys[j + 1].u[z] + h*k3, sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].v[z] + h*l3, i, j,z);
			l4 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].v[z] + h*l3);
			m4 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].ta[z] + h*m3);
		}
		if (i == 0 && j == m - 1)
		{
			k1 = x_corners(sys[i].sys[j - 1].u[z], sys[i + 1].sys[j].u[z], sys[i].sys[j].u[z], sys[i].sys[j].v[z], i, j,z);
			l1 = sys[i].sys[j].g(sys[i].sys[j].u[z], sys[i].sys[j].v[z]);
			m1 = sys[i].sys[j].th(sys[i].sys[j].u[z], sys[i].sys[j].ta[z]);

			k2 = x_corners(sys[i].sys[j - 1].u[z] + h*k1 / 2, sys[i + 1].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2, i, j,z);
			l2 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2);
			m2 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].ta[z] + h*m1 / 2);

			k3 = x_corners(sys[i].sys[j - 1].u[z] + h*k2 / 2, sys[i + 1].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l2 / 2, i, j,z);
			l3 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l2 / 2);
			m3 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].ta[z] + h*m2 / 2);

			k4 = x_corners(sys[i].sys[j - 1].u[z] + h*k3, sys[i + 1].sys[j].u[z] + h*k3, sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].v[z] + h*l3, i, j,z);
			l4 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].v[z] + h*l3);
			m4 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].ta[z] + h*m3);
		}
		if (i == 0 && j == 0)
		{
			k1 = x_corners(sys[0].sys[1].u[z], sys[1].sys[0].u[z], sys[0].sys[0].u[z], sys[0].sys[0].v[z], 0, 0,z);
			l1 = sys[0].sys[0].g(sys[0].sys[0].u[z], sys[0].sys[0].v[z]);
			m1 = sys[0].sys[0].th(sys[0].sys[0].u[z], sys[0].sys[0].ta[z]);

			k2 = x_corners(sys[0].sys[1].u[z] + h*k1 / 2, sys[1].sys[0].u[z] + h*k1 / 2, sys[0].sys[0].u[z] + h*k1 / 2, sys[0].sys[0].v[z] + h*l1 / 2, 0, 0,z);
			l2 = sys[0].sys[0].g(sys[0].sys[0].u[z] + h*k1 / 2, sys[0].sys[0].v[z]+h*l1/2);
			m2 = sys[0].sys[0].th(sys[0].sys[0].u[z] + h*k1 / 2, sys[0].sys[0].ta[z] + h*m1 / 2);

			k3 = x_corners(sys[0].sys[1].u[z] + h*k2 / 2, sys[1].sys[0].u[z] + h*k2 / 2, sys[0].sys[0].u[z] + h*k2 / 2, sys[0].sys[0].v[z] + h*l2 / 2, 0, 0,z);
			l3 = sys[0].sys[0].g(sys[0].sys[0].u[z] + h*k2 / 2, sys[0].sys[0].v[z] + h*l2 / 2);
			m3 = sys[0].sys[0].th(sys[0].sys[0].u[z] + h*k2 / 2, sys[0].sys[0].ta[z] + h*m2 / 2);

			k4 = x_corners(sys[0].sys[1].u[z] + h*k3, sys[1].sys[0].u[z] + h*k3, sys[0].sys[0].u[z] + h*k3, sys[0].sys[0].v[z] + h*l3, 0, 0,z);
			l4 = sys[0].sys[0].g(sys[0].sys[0].u[z] + h*k3, sys[0].sys[0].v[z] + h*l2);
			m4 = sys[0].sys[0].th(sys[0].sys[0].u[z] + h*k3, sys[0].sys[0].ta[z] + h*m3);
		}
		if (i == k - 1 && j == m - 1)
		{
			k1 = x_corners(sys[k - 2].sys[m - 1].u[z], sys[k - 1].sys[m - 2].u[z], sys[k - 1].sys[m - 1].u[z], sys[k - 1].sys[m - 1].v[z], k - 1, m - 1,z);
			l1 = sys[k - 1].sys[m - 1].g(sys[k - 1].sys[m - 1].u[z], sys[k - 1].sys[m - 1].v[z]);
			m1 = sys[k-1].sys[m-1].th(sys[k-1].sys[m-1].u[z], sys[k-1].sys[m-1].ta[z]);

			k2 = x_corners(sys[k - 2].sys[m - 1].u[z] + h*k1 / 2, sys[k - 1].sys[m - 2].u[z] + h*k1 / 2, sys[k - 1].sys[m - 1].u[z] + h*k1 / 2, sys[k - 1].sys[m - 1].v[z] + h*l1 / 2, k - 1, m - 1,z);
			l2 = sys[k - 1].sys[m - 1].g(sys[k - 1].sys[m - 1].u[z] + h*k1 / 2, sys[k - 1].sys[m - 1].v[z]+h*l1/2);
			m2 = sys[k - 1].sys[m - 1].th(sys[k - 1].sys[m - 1].u[z]+h*k1/2, sys[k - 1].sys[m - 1].ta[z]+h*m1/2);

			k3 = x_corners(sys[k - 2].sys[m - 1].u[z] + h*k2 / 2, sys[k - 1].sys[m - 2].u[z] + h*k2 / 2, sys[k - 1].sys[m - 1].u[z] + h*k2 / 2, sys[k - 1].sys[m - 1].v[z] + h*l2 / 2, k - 1, m - 1,z);
			l3 = sys[k - 1].sys[m - 1].g(sys[k - 1].sys[m - 1].u[z] + h*k2 / 2, sys[k - 1].sys[m - 1].v[z] + h*l2 / 2);
			m3 = sys[k - 1].sys[m - 1].th(sys[k - 1].sys[m - 1].u[z] + h*k2 / 2, sys[k - 1].sys[m - 1].ta[z] + h*m2 / 2);

			k4 = x_corners(sys[k - 2].sys[m - 1].u[z] + h*k3, sys[k - 1].sys[m - 2].u[z] + h*k3, sys[k - 1].sys[m - 1].u[z] + h*k3, sys[k - 1].sys[m - 1].v[z] + h*l3, k - 1, m - 1,z);
			l4 = sys[k - 1].sys[m - 1].g(sys[k - 1].sys[m - 1].u[z] + h*k3, sys[k - 1].sys[m - 1].v[z] + h*l3);
			m4 = sys[k - 1].sys[m - 1].th(sys[k - 1].sys[m - 1].u[z] + h*k3, sys[k - 1].sys[m - 1].ta[z] + h*m3);
		}
		sys[i].sys[j].u[z + 1] = sys[i].sys[j].u[z] + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
		sys[i].sys[j].v[z + 1] = sys[i].sys[j].v[z] + h*(l1 + 2 * l2 + 2 * l3 + l4) / 6;
		sys[i].sys[j].ta[z + 1] = sys[i].sys[j].ta[z] + h*(m1 + 2 * m2 + 2 * m3 + m4) / 6;
	}

	void MethodRunge_Kutta_Edges(int i, int j, int z)
	{
		 double k1, k2, k3, k4;
		 double l1, l2, l3, l4;
		 double m1, m2, m3, m4;
		if (i == 0)
		{
			k1 = x_edges(sys[i].sys[j - 1].u[z], sys[i].sys[j + 1].u[z], sys[i + 1].sys[j].u[z], sys[i].sys[j].u[z], sys[i].sys[j].v[z], i, j,z);
			l1 = sys[i].sys[j].g(sys[i].sys[j].u[z],sys[i].sys[j].v[z]);
			m1 = sys[i].sys[j].th(sys[i].sys[j].u[z], sys[i].sys[j].ta[z]);

			k2 = x_edges(sys[i].sys[j - 1].u[z] + h*k1 / 2, sys[i].sys[j + 1].u[z] + h*k1 / 2, sys[i + 1].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2, i, j,z);
			l2 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2);
			m2 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].ta[z] + h*m1 / 2);

			k3 = x_edges(sys[i].sys[j - 1].u[z] + h*k2 / 2, sys[i].sys[j + 1].u[z] + h*k2 / 2, sys[i + 1].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l2 / 2, i, j,z);
			l3 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l1 / 2);
			m3 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].ta[z] + h*m2 / 2);

			k4 = x_edges(sys[i].sys[j - 1].u[z] + h*k3, sys[i].sys[j + 1].u[z] + h*k3, sys[i + 1].sys[j].u[z] + h*k3, sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].u[z] + h*l3, i, j,z);
			l4 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k3,sys[i].sys[j].v[z]+h*l3);
			m4 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].ta[z] + h*m3);
		}
		if (j == 0)
		{
			k1 = x_edges(sys[i - 1].sys[j].u[z], sys[i].sys[j + 1].u[z], sys[i + 1].sys[j].u[z], sys[i].sys[j].u[z], sys[i].sys[j].v[z], i, j,z);
			l1 = sys[i].sys[j].g(sys[i].sys[j].u[z], sys[i].sys[j].v[z]);
			m1 = sys[i].sys[j].th(sys[i].sys[j].u[z], sys[i].sys[j].ta[z]);

			k2 = x_edges(sys[i - 1].sys[j].u[z] + h*k1 / 2, sys[i].sys[j + 1].u[z] + h*k1 / 2, sys[i + 1].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2, i, j,z);
			l2 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2);
			m2 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].ta[z] + h*m1 / 2);

			k3 = x_edges(sys[i - 1].sys[j].u[z] + h*k2 / 2, sys[i].sys[j + 1].u[z] + h*k2 / 2, sys[i + 1].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l2 / 2, i, j,z);
			l3 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l1 / 2);
			m3 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].ta[z] + h*m2 / 2);

			k4 = x_edges(sys[i - 1].sys[j].u[z] + h*k3, sys[i].sys[j + 1].u[z] + h*k3, sys[i + 1].sys[j].u[z] + h*k3, sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].v[z] + h*l3, i, j,z);
			l4 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].v[z] + h*l3);
			m4 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].ta[z] + h*m3);
		}
		if (i == k - 1)
		{
			k1 = x_edges(sys[i].sys[j - 1].u[z], sys[i].sys[j + 1].u[z], sys[i - 1].sys[j].u[z], sys[i].sys[j].u[z], sys[i].sys[j].v[z], i, j,z);
			l1 = sys[i].sys[j].g(sys[i].sys[j].u[z], sys[i].sys[j].v[z]);
			m1 = sys[i].sys[j].th(sys[i].sys[j].u[z], sys[i].sys[j].ta[z]);

			k2 = x_edges(sys[i].sys[j - 1].u[z] + h*k1 / 2, sys[i].sys[j + 1].u[z] + h*k1 / 2, sys[i - 1].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2, i, j,z);
			l2 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2);
			m2 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].ta[z] + h*m1 / 2);

			k3 = x_edges(sys[i].sys[j - 1].u[z] + h*k2 / 2, sys[i].sys[j + 1].u[z] + h*k2 / 2, sys[i - 1].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l2 / 2, i, j,z);
			l3 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l1 / 2);
			m3 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].ta[z] + h*m2 / 2);

			k4 = x_edges(sys[i].sys[j - 1].u[z] + h*k3, sys[i].sys[j + 1].u[z] + h*k3, sys[i - 1].sys[j].u[z] + h*k3, sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].v[z] + h*l3, i, j,z);
			l4 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].v[z] + h*l3);
			m4 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].ta[z] + h*m3);
		}
		if (j == m - 1)
		{
			k1 = x_edges(sys[i - 1].sys[j].u[z], sys[i].sys[j - 1].u[z], sys[i + 1].sys[j].u[z], sys[i].sys[j].u[z], sys[i].sys[j].v[z], i, j,z);
			l1 = sys[i].sys[j].g(sys[i].sys[j].u[z], sys[i].sys[j].v[z]);
			m1 = sys[i].sys[j].th(sys[i].sys[j].u[z], sys[i].sys[j].ta[z]);

			k2 = x_edges(sys[i - 1].sys[j].u[z] + h*k1 / 2, sys[i].sys[j - 1].u[z] + h*k1 / 2, sys[i + 1].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2, i, j,z);
			l2 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2);
			m2 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].ta[z] + h*m1 / 2);

			k3 = x_edges(sys[i - 1].sys[j].u[z] + h*k2 / 2, sys[i].sys[j - 1].u[z] + h*k2 / 2, sys[i + 1].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l2 / 2, i, j,z);
			l3 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l1 / 2);
			m3 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].ta[z] + h*m2 / 2);

			k4 = x_edges(sys[i - 1].sys[j].u[z] + h*k3, sys[i].sys[j - 1].u[z] + h*k3, sys[i + 1].sys[j].u[z] + h*k3, sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].v[z] + h*l3, i, j,z);
			l4 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].v[z] + h*l3);
			m4 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].ta[z] + h*m3);
		}

		sys[i].sys[j].u[z + 1] = sys[i].sys[j].u[z] + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
		sys[i].sys[j].v[z + 1] = sys[i].sys[j].v[z] + h*(l1 + 2 * l2 + 2 * l3 + l4) / 6;
		sys[i].sys[j].ta[z + 1] = sys[i].sys[j].ta[z] + h*(m1 + 2 * m2 + 2 * m3 + m4) / 6;
	}


	void MethodRunge_Kutta_Inner(int i, int j, int z)
	{
		 double k1, k2, k3, k4;
		 double l1, l2, l3, l4;
		 double m1, m2, m3, m4;

		k1 = x_inner(sys[i].sys[j - 1].u[z], sys[i + 1].sys[j].u[z], sys[i - 1].sys[j].u[z], sys[i].sys[j + 1].u[z], sys[i].sys[j].u[z], sys[i].sys[j].v[z], i, j,z);
		l1 = sys[i].sys[j].g(sys[i].sys[j].u[z], sys[i].sys[j].v[z]);
		m1 = sys[i].sys[j].th(sys[i].sys[j].u[z], sys[i].sys[j].ta[z]);

		k2 = x_inner(sys[i].sys[j - 1].u[z] + h*k1 / 2, sys[i + 1].sys[j].u[z] + h*k1 / 2, sys[i - 1].sys[j].u[z] + h*k1 / 2, sys[i].sys[j + 1].u[z] + h*k1 / 2, sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2, i, j,z);
		l2 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].v[z] + h*l1 / 2);
		m2 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k1 / 2, sys[i].sys[j].ta[z] + h*m1 / 2);

		k3 = x_inner(sys[i].sys[j - 1].u[z] + h*k2 / 2, sys[i + 1].sys[j].u[z] + h*k2 / 2, sys[i - 1].sys[j].u[z] + h*k2 / 2, sys[i].sys[j + 1].u[z] + h*k2 / 2, sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l2 / 2, i, j,z);
		l3 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].v[z] + h*l1 / 2);
		m3 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k2 / 2, sys[i].sys[j].ta[z] + h*m2 / 2);

		k4 = x_inner(sys[i].sys[j - 1].u[z] + h*k3, sys[i + 1].sys[j].u[z] + h*k3, sys[i - 1].sys[j].u[z] + h*k3, sys[i].sys[j + 1].u[z] + h*k3, sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].u[z] + h*l3, i, j,z);
		l4 = sys[i].sys[j].g(sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].v[z] + h*l3);
		m4 = sys[i].sys[j].th(sys[i].sys[j].u[z] + h*k3, sys[i].sys[j].ta[z] + h*m3);

		sys[i].sys[j].u[z + 1] = sys[i].sys[j].u[z] + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
		sys[i].sys[j].v[z + 1] = sys[i].sys[j].v[z] + h*(l1 + 2 * l2 + 2 * l3 + l4) / 6;
		sys[i].sys[j].ta[z + 1] = sys[i].sys[j].ta[z] + h*(m1 + 2 * m2 + 2 * m3 + m4) / 6;
	}

	void MethodRunge_Kutta(int z)//вычисление методом Цунге-†утта 4-го порІдка
	{

		MethodRunge_Kutta_Corners(0, 0, z);
		MethodRunge_Kutta_Corners(0, m - 1, z);
		MethodRunge_Kutta_Corners(k - 1, 0, z);
		MethodRunge_Kutta_Corners(k - 1, m - 1, z);
		for (int i = 1; i < k - 1; i++)
			for (int j = 1; j < m - 1; j++)
				MethodRunge_Kutta_Inner(i, j, z);
		for (int j = 1; j < m - 1; j++)
			MethodRunge_Kutta_Edges(0, j, z);
		for (int j = 1; j < k - 1; j++)
			MethodRunge_Kutta_Edges(j, 0, z);
		for (int j = 1; j < m - 1; j++)
			MethodRunge_Kutta_Edges(k - 1, j, z);
		for (int j = 1; j < k - 1; j++)
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
#pragma once
