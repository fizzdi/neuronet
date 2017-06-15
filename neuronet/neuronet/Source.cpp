#define _CRT_SECURE_NO_WARNINGS
#include "NeuralNetwork.h"
#include <algorithm>
#include <ctime>
#include <iomanip>
using namespace std;
const int num_check = 15;
const int HiddenNeuron = 15; //sin - 11 15
const int SampleCount = 10;
const int Epoh = 10000;

double testfun(double x)
{
	return sin(x);
}

const double xmin = -1;
const double xmax = 1;

int main()
{
	std::ios::sync_with_stdio(false);
	freopen("debug.txt", "w+", stdout);
	NeuroNet::NeuralNetwork net;
	net.Xmin = -2;
	net.Xmax = 2;
	net.Ymin = -1;
	net.Ymax = 1;
	//srand(1);
	srand(time(NULL));

	bool failed = false;
	net.Init(1, 1, HiddenNeuron, NeuroNet::SIGM, NeuroNet::RPROP);
	for (int i = 1; i <= SampleCount; ++i)
	{
		double x = xmin + rand()*(xmax - xmin) / RAND_MAX;
		//double x =rand()%10;
		//double x = i;
		double y = testfun(x);
		net.TrainingSet.push_back(NeuroNet::Problem({ x }, { y}));
	}
	cout.flush();
	double lstError = 0.0;
	double Error = 0.0;
	int curEpoh = 1;
	do
	{
		lstError = Error;
		bool fl = curEpoh % 500 == 0;
		//Error = net.RunTrainingSet();
		Error = net.RunTrainingSetOffline(0);
		if (fl)
		{
			cout << endl << "================================================================================================" << endl;
			cout << setw(15) << curEpoh << ") " << setw(15) << Error << " (" << setw(15) << abs(Error - lstError) << ")" << endl;
			cout.flush();
		}
		if (Error < 1e-6 || abs(Error - lstError) < 1e-12)
		{
			cout << "STOP Epoh: " << setw(15) << curEpoh << " " << setw(15) << Error << " (" << setw(15) << abs(Error - lstError) << ")" << endl;
			break;
		}

		auto y = -nan("");
		if (Error > 1e12 /*|| !(!(Error < -nan("")) && !(Error < -nan("")))*/)
		{
			cout << "Failed" << endl;
			failed = true;
			break;
		}
		curEpoh++;
	} while (curEpoh < Epoh);

	cout << "End Epoh: " << setw(15) << curEpoh << " Error:" << setw(15) << Error << " (" << setw(15) << abs(Error - lstError) << ")" << endl;


	if (!failed)
	{
		cout << endl << "================================================================================================" << endl;
		cout << "Tests:" << endl;
		vector<double> inputs(1);
		{
			double mse = 0.0;
			double max_mse = 0.0;
			double max_error = 0.0;
			for (int i = 1; i <= num_check; ++i)
			{
				//inputs[0] = i;
				inputs[0] = xmin + rand()*(xmax - xmin) / RAND_MAX;
				//inputs[0] = rand() % 100;

				double ideal = testfun(inputs[0]);
				net.Run(NeuroNet::Problem({ inputs[0] }, { testfun(inputs[0]) }));
				auto outputs = net.GetOut();
				double curerror = abs(outputs[0][0] - ideal) * abs(outputs[0][0] - ideal) / 2;
				mse += curerror;
				max_error = max(max_error, abs(outputs[0][0] - ideal));
				cout << setw(15) << inputs[0] << " " << setw(15) << outputs[0][0] << " (" << setw(15) << ideal << ") error: " << setw(15) << abs(outputs[0][0] - ideal)  << " mse: " << setw(15) << curerror << endl;
				max_mse = max(max_mse, curerror);
			}
			cout << endl << "Sum mse: " << mse << endl;
			cout << "MAX mse: " << max_mse << endl;
			cout << "Avg mse: " << mse / num_check << endl;
			cout << "MAX error: " << max_error << endl;


		}
	}

	//system("pause");
	return 0;
}