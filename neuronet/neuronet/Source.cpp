#include "NeuralNetwork.h"
#include <algorithm>
#include <ctime>
using namespace std;
const int num_check = 150;
const int HiddenNeuron = 5; //sin - 11 15
const int SampleCount = 15;
const int Epoh = 200;

double testfun(double x)
{
	return sin(x);
}
int main()
{
	std::ios::sync_with_stdio(false);
	NeuroNet::NeuralNetwork net;
	//srand(1);
	srand(time(NULL));

	bool failed = false;
	net.Init(1, 1, HiddenNeuron, NeuroNet::TANH, NeuroNet::RPROP);
	for (int i = 0; i < SampleCount; ++i)
	{
		double x = -1+ rand()*2.0/RAND_MAX;
		double y = testfun(x);
		net.TrainingSet.push_back(NeuroNet::Problem({ x }, { y}));
	}
	double lstError = 0.0;
	double Error = 0.0;
	int curEpoh = 1;
	do
	{
		lstError = Error;
		bool fl = curEpoh % 10 == 0;
		Error = net.RunTrainingSet(true);
		if (fl)
		{
			cout << endl << "================================================================================================" << endl;
			cout << curEpoh << ") " << Error << " (" << abs(Error - lstError) << ")" << endl;
			cout.flush();
		}
		if (Error < 1e-6 || abs(Error - lstError) < 1e-12)
		{
			cout << "STOP Epoh: " << curEpoh << " " << Error << " (" << abs(Error - lstError) << ")" << endl;
			break;
		}

		auto y = -nan("");
		if (Error > 1e9 /*|| !(!(Error < -nan("")) && !(Error < -nan("")))*/)
		{
			cout << "Failed" << endl;
			failed = true;
			break;
		}
		curEpoh++;
	} while (curEpoh < Epoh);

	cout << "End Epoh: " << curEpoh << " Error:" << Error << " (" << abs(Error - lstError) << ")" << endl;


	if (!failed)
	{
		vector<double> inputs(1);
		{
			double error = 0.0;
			double max_error = 0.0;
			for (int i = 0; i < num_check; ++i)
			{
				inputs[0] = -1 + rand()*2.0 / RAND_MAX;
				double ideal = testfun(inputs[0]);
				NeuroNet::Matrix2d inpm;
				inpm = inputs;
				net.Run(inpm);
				auto outputs = net.GetOut();
				double curerror = abs(outputs[0][0] - ideal) * abs(outputs[0][0] - ideal) / 2;
				error += curerror;
				cout << inputs[0] << " " << outputs[0][0] << " (" << ideal << ") error: " << abs(outputs[0][0] - ideal)  << " mse: "<< curerror << endl;
				max_error = max(max_error, curerror);
			}
			cout << endl << "Sum Error: " << error << endl;
			cout << "MAX Error: " << max_error << endl;
			cout << "Avg Error: " << error / num_check << endl;

		}
	}

	system("pause");
	return 0;
}