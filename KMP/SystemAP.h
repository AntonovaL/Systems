#ifndef SYSTEMAP_H
#define SYSTEMAP_H


#include <vector>
#include <math.h>

using namespace std;

class SystemAP {
	int n;
	float h, t0, tn;
	int rezh;
	float IextConst;
	float TextConst;
	vector<float> Iext;
	vector<float> Text;
	float E = 1;
	float Ue = 0.05;
	float Es = 1;
	float Gs = 0;
	bool mex;
	bool elec;
public:
	vector< double> u;
	vector< double> v;
	vector< double> ta;
	float a, eps0, eps1, k, kt;
	SystemAP() {
		a = 0.08;
		eps0 = 1;
		eps1 = 0.1;
		k = 8;
		kt = 1.5;
		rezh = 0; mex = false; elec = false;
		TextConst = 0; IextConst = 0;
		t0 = 0; tn = 100;
		h = 0.01;
		n = (tn - t0) / h;
		u.resize(n);
		v.resize(n);
		ta.resize(n);
		Iext.resize(n);
		Text.resize(n);
		for (int i = 0; i < n; i++) {
			u[i] = 0;
			v[i] = 0;
			ta[i] = 0;
			Iext[i] = 0;
			Text[i] = 0;
		}
	}

	void setParametrs(float _a, float _k, float _kt, float _eps1, float _eps0) {
		a = _a;
		k = _k;
		kt = _kt;
		eps1 = _eps1;
		eps0 = _eps0;
	}

	void set(float _h, float _t0, float _tn) {
		if (_h != h && t0 != _t0 && tn != _tn) {
			h = _h;
			t0 = _t0;
			tn = _tn;
			n = (tn - t0) / h;
			u.resize(n);
			v.resize(n);
			ta.resize(n);
			Iext.resize(n);
			Text.resize(n);
			for (int i = 0; i < n; i++) {
				u[i] = 0;
				v[i] = 0;
				ta[i] = 0;
				Iext[i] = 0;
				Text[i] = 0;
			}
		}
	}

	void setInitialCondition(double _u0, double _v0, double _ta0) {
		u[0] = _u0;
		v[0] = _v0;
		ta[0] = _ta0;
	}
	void setElec(float _Iext,int i) {
		Iext[i] = _Iext;
	}
	void setMex(float _Text,int i) {
		Text[i] = _Text;
	}
	void setElecConst(float _Iext) {
		elec = true;
		IextConst = _Iext;
	}
	void setMexConst(float _Text) {
		mex = true;
		TextConst = _Text;
	}
	int getN() {
		return n;
	}
	void SetRezh(int _rezh) {
		rezh = _rezh;
	}

	 double FuncHeviside( double u) {
		if (u >= 0) {
			return 1;
		}
		else {
			return 0;
		}
	}
	float Eps(long double u) {
		return eps0 + FuncHeviside(u - Ue)*(eps1 - eps0);
	}

	 double calcD(int i) {
		if (rezh == 0) {
			return 0;
		}
		else {
			if (mex) {
				return (TextConst - ta[i]) / E;
			}
			else {
				return (Text[i] - ta[i]) / E;
			}
		}
	}

	 double CalcIs(int i) {
		if (rezh == 0) {
			return 0;
		}
		else {
			return Gs*FuncHeviside(calcD(i))*calcD(i)*(u[i] - Es);
		}
	}

	 double calcIext(int i) {
		if (elec) {
			return IextConst;
		}
		else {
			return Iext[i];
		}
	}
	 double f( double ui,  double vi, int i) {
		double tmp = -k*ui*(ui - a)*(ui - 1) - ui * vi-CalcIs(i)+calcIext(i);
		 //double tmp = -k*ui*(ui - a)*(ui - 1) - ui * vi;
		return tmp;
	}

	 double g( double ui, double vi) {
		return Eps(ui)*(k*ui - vi);
	}

	 double th( double ui,  double tai) {
		return Eps(ui)*(kt*ui - tai);
	}


	void methodRungeKutta() {
		 double k1, k2, k3, k4;
		 double l1, l2, l3, l4;
		 double m1, m2, m3, m4;
		for (int i = 0; i < n - 1; i++) {
			k1 = f(u[i], v[i],i);
			l1 = g(u[i], v[i]);
			m1 = th(u[i], ta[i]);

			k2 = f(u[i] + h*k1 / 2, v[i] + h*l1 / 2,i);
			l2 = g(u[i] + h*k1 / 2, v[i] + h*l1 / 2);
			m2 = th(u[i]+h*k1/2, ta[i]+h*m1/2);

			k3 = f(u[i] + h*k2 / 2, v[i] + h*l2 / 2, i);
			l3 = g(u[i] + h*k2 / 2, v[i] + h*l2 / 2);
			m3 = th(u[i] + h*k2 / 2, ta[i] + h*m2 / 2);

			k4 = f(u[i] + h*k3, v[i] + h*l3, i);
			l4 = g(u[i] + h*k3, v[i] + h*l3);
			m4 = th(u[i] + h*k3, ta[i] + h*m3);

			u[i + 1] = u[i] + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
			v[i + 1] = v[i] + h*(l1 + 2 * l2 + 2 * l3 + l4) / 6;
			ta[i + 1] = ta[i] + h*(m1 + 2 * m2 + 2 * m3 + m4) / 6;
		}
	}
};
#endif
