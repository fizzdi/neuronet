#include "NeuralNetwork.h"
#include <algorithm>
#include <ctime>
using namespace std;
const int num_check = 30;
const int HiddenNeuron =20; //sin - 30
const int SampleCount = 5;
const int Epoh = 100000;
#define TESTFUNC sin
int main()
{
	std::ios::sync_with_stdio(false);
	NeuroNet::NeuralNetwork net;
	//srand(1);
	srand(time(NULL));

	net.Init(1, 1, HiddenNeuron, NeuroNet::TANH, NeuroNet::BACKPROP);
	for (int i = 0; i < SampleCount; ++i)
	{
		double x = -1+rand()*2.0 / RAND_MAX;
		double y = TESTFUNC(x);
		net.TrainingSet.push_back(NeuroNet::Problem({ x }, { y}));
	}
	double lstError = 0.0;
	double maxError = 0.0;
	int curEpoh = 1;
	do
	{
		lstError = maxError;
		maxError = net.RunTrainingSet();
		bool fl = curEpoh % 1000 == 0;
		if (fl)
		{
			cout << endl << "================================================================================================" << endl;
			cout << curEpoh << ") " << maxError << " (" << abs(maxError - lstError) << ")" << endl;
		}
		if (maxError < 1e-6 || abs(maxError - lstError) < 1e-10)
		{
			cout << "STOP Epoh: " << curEpoh << " " << maxError << " (" << abs(maxError - lstError) << ")" << endl;
			break;
		}
		curEpoh++;
	} while (curEpoh < Epoh);

	vector<double> inputs(1);
	{
		double error = 0.0;
		double max_error = 0.0;
		for (int i = 0; i < num_check; ++i)
		{
			inputs[0] = -1+rand()*2.0 / RAND_MAX;
			double ideal = TESTFUNC(inputs[0]);
			NeuroNet::Matrix2d inpm;
			inpm = inputs;
			net.Run(inpm);
			auto outputs = net.GetOut();
			double curerror = abs(outputs[0][0] - ideal);
			error += curerror;
			cout << inputs[0] << " " << outputs[0][0] << " (" << ideal << ") error: " << curerror << endl;
			max_error = max(max_error, curerror);
		}
		cout << endl << "Sum Error: " << error << endl;
		cout << "MAX Error: " << max_error << endl;
		cout << "Avg Error: " << error / num_check << endl;

	}

	system("pause");
	return 0;
}