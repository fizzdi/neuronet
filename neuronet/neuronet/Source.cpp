#include "NeuralNetwork.h"
#include "Elman.h"
#include <algorithm>
#include <ctime>
using namespace std;
const int num_check = 30;
const int HiddenNeuron =20; //sin - 30
const int SampleCount = 20;
const int Epoh = 10000;
#define TESTFUNC sin
int main()
{
	std::ios::sync_with_stdio(false);
	NeuroNet::NeuralNetwork net;
	srand(1);	//srand(time(NULL));

	net.Init(1, 1, HiddenNeuron, NeuroNet::TANH);
	for (int i = 0; i < SampleCount; ++i)
	{
		double x = rand()*1.0 / RAND_MAX;
		double y = TESTFUNC(x);
		net.TrainingSet.push_back(NeuroNet::Problem({ x }, { y}));
	}
	net.RunTrainingSet();
	double lstError = 0.0;
	double maxError = 0.0;
	for (int i = 0; i < Epoh; ++i)
	{
		lstError = maxError;
		net.CorrectWeights();
		//cout << endl << endl << "CORRECT" << endl << endl;
		maxError = net.RunTrainingSet();
		if (i % 1000 == 0)
			cout << i << ") " << maxError << endl;
		if (maxError < 0.000001 || abs(maxError - lstError) < 0.00000001)
		{
			cout << i << endl;
			break;
		}
	}

	vector<double> inputs(1);
	{
		double error = 0.0;
		double max_error = 0.0;
		for (int i = 0; i < num_check; ++i)
		{
			inputs[0] = rand()*1.0 / RAND_MAX;
			double ideal = TESTFUNC(inputs[0]);
			NeuroNet::Matrix2d inpm;
			inpm.operator=(inputs);
			auto outputs = net.Run(inpm);
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