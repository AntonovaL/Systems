#include "SystemFHN.h"
#include "ChainFHN.h"
#include "SquareFHN.h"
#include "SystemAP.h"
#include "SystemChainAP.h"
#include "SystemSquareAP.h"
#include <iostream>
#include <fstream>
#include "iomanip" 
using namespace std;

int main() {
	/*ChainFHN system;
	ofstream fout;
	double x0 = 0, y0 = 0;
	x0 = -1.1; y0 = -1.1 - pow(-1.1, 3) / 3;
	std::ofstream result_stream("testing.bin", ios::out | ios::binary);
	for (int i = 0; i < 100; i++)
	{
		system.sys[i].setParametrs(-1.1, 0.05);
		if (i<10)
			system.sys[i].setInitialCondition(2, -1.5);
		else
			system.sys[i].setInitialCondition(x0, y0);
	}

	system.Calculate();
	for (int i = 0; i < 100; i++)
	{

		for (int j = 0; j < 10000; j = j + 100)
		{
			result_stream << system.sys[i].x[j] << "  ";
		}
		result_stream << endl;
	}
	result_stream.close();*/
	/*SquareFHN system;
	ofstream fout;
	double x0 = 0, y0 = 0;
	x0 = -1.1; y0 = -1.1 - pow(-1.1, 3) / 3;
	std::ofstream result_stream1("test1.bin", ios::out | ios::binary);
	std::ofstream result_stream2("test2.bin", ios::out | ios::binary);
	std::ofstream result_stream3("test3.bin", ios::out | ios::binary);
	std::ofstream result_stream4("test4.bin", ios::out | ios::binary);
	std::ofstream result_stream5("test5.bin", ios::out | ios::binary);
	std::ofstream result_stream6("test6.bin", ios::out | ios::binary);
	std::ofstream result_stream7("test7.bin", ios::out | ios::binary);
	std::ofstream result_stream8("test8.bin", ios::out | ios::binary);
	for (int i = 0; i < 100; i++)
	{
		for (int j = 0; j < 100; j++)
		{
			system.sys[i].sys[j].setParametrs(-1.1, 0.05);
			if (i < 3 && j < 3)
				system.sys[i].sys[j].setInitialCondition(2, -1.5);
			else
				system.sys[i].sys[j].setInitialCondition(x0, y0);

		}
	}

	system.CalculationFunctions();
	for (int i = 0; i < 100; i++)
	{
		for (int j = 0; j < 100; j = j++)
		{
			result_stream1 <<system.sys[i].sys[j].x[0] << " ";
			result_stream2 <<  system.sys[i].sys[j].x[100] << " ";
			result_stream3 <<  system.sys[i].sys[j].x[300] << " ";
			result_stream4 <<  system.sys[i].sys[j].x[500] << " ";
			result_stream5 <<  system.sys[i].sys[j].x[600] << " ";
			result_stream6 << system.sys[i].sys[j].x[700] << " ";
			result_stream7 << system.sys[i].sys[j].x[800] << " ";
			result_stream8 << system.sys[i].sys[j].x[900] << " ";
		}
		result_stream1 << endl;
		result_stream2 << endl;
		result_stream3 << endl;
		result_stream4 << endl;
		result_stream5 << endl;
		result_stream6 << endl;
		result_stream7 << endl;
		result_stream8 << endl;

	}
	result_stream1.close();
	result_stream2.close();
	result_stream3.close();
	result_stream4.close();
	result_stream5.close();
	result_stream6.close();
	result_stream7.close();
	result_stream8.close();*/
	/*ChainAP system;
	ofstream fout;
	std::ofstream result_stream("testing.bin", ios::out | ios::binary);
	for (int i = 0; i < 100; i++)
	{
		system.sys[i].setParametrs(0.08, 8, 1.5, 0.1, 1);
		if (i < 6)
			system.sys[i].setElecConst(0.3); //system.sys[i].setInitialCondition(0.8,0,0);
		system.sys[i].setInitialCondition(0, 0, 0);
		//else
			//system.sys[i].setInitialCondition(0, 0,0);
	}

	system.Calculate();
	for (int i = 0; i < 100; i++)
	{

		for (int j = 0; j < 10000; j = j + 100)
		{
			result_stream << system.sys[i].u[j] << "  ";
		}
		result_stream << endl;
	}
	result_stream.close();*/
	SquareAP system;
	ofstream fout;
	std::ofstream result_stream1("test1.bin", ios::out | ios::binary);
	std::ofstream result_stream2("test2.bin", ios::out | ios::binary);
	std::ofstream result_stream3("test3.bin", ios::out | ios::binary);
	std::ofstream result_stream4("test4.bin", ios::out | ios::binary);
	std::ofstream result_stream5("test5.bin", ios::out | ios::binary);
	std::ofstream result_stream6("test6.bin", ios::out | ios::binary);
	std::ofstream result_stream7("test7.bin", ios::out | ios::binary);
	std::ofstream result_stream8("test8.bin", ios::out | ios::binary);
	for (int i = 0; i < 100; i++)
	{
	for (int j = 0; j < 100; j++)
	{
	system.sys[i].sys[j].setParametrs(0.08, 8, 1.5, 0.1, 1);
	if (i < 3 && j < 3)
	system.sys[i].sys[j].setInitialCondition(0.8, 0, 0);
	else
	system.sys[i].sys[j].setInitialCondition(0,0,0);

	}
	}

	system.CalculationFunctions();
	for (int i = 0; i < 100; i++)
	{
	for (int j = 0; j < 100; j = j++)
	{
	result_stream1 <<system.sys[i].sys[j].u[0] << " ";
	result_stream2 <<  system.sys[i].sys[j].u[100] << " ";
	result_stream3 <<  system.sys[i].sys[j].u[300] << " ";
	result_stream4 <<  system.sys[i].sys[j].u[500] << " ";
	result_stream5 <<  system.sys[i].sys[j].u[600] << " ";
	result_stream6 << system.sys[i].sys[j].u[700] << " ";
	result_stream7 << system.sys[i].sys[j].u[800] << " ";
	result_stream8 << system.sys[i].sys[j].u[900] << " ";
	}
	result_stream1 << endl;
	result_stream2 << endl;
	result_stream3 << endl;
	result_stream4 << endl;
	result_stream5 << endl;
	result_stream6 << endl;
	result_stream7 << endl;
	result_stream8 << endl;

	}
	result_stream1.close();
	result_stream2.close();
	result_stream3.close();
	result_stream4.close();
	result_stream5.close();
	result_stream6.close();
	result_stream7.close();
	result_stream8.close();
}