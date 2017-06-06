#include "NeuralNetwork.h"
#include "Elman.h"
#include <ctime>
using namespace std;
const int num_check = 10;
const int HiddenNeuron = 20; //sin - 30
const int SampleCount = 20;
const int Epoh = 1000;
#define TESTFUNC sin
int main()
{
	std::ios::sync_with_stdio(false);
	NeuroNet::Elman net;
	srand(time(NULL));
	/*net.Init(2, 1, 2);
	net.TrainingSet.push_back(NeuroNet::Problem({ 1,1 }, { 0 }));
	net.TrainingSet.push_back(NeuroNet::Problem({ 0,0 }, { 0 }));
	net.TrainingSet.push_back(NeuroNet::Problem({ 0,1 }, { 1 }));
	net.TrainingSet.push_back(NeuroNet::Problem({ 1,0 }, { 1 }));
	*/

	net.Init(1, 1, HiddenNeuron);
	for (int i = 0; i < SampleCount; ++i)
	{
		double x = rand()*1.0 / RAND_MAX;
		double y = TESTFUNC(x);
		net.TrainingSet.push_back(NeuroNet::Problem({ x }, { y }));
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
		cout << i << ") " << maxError << endl;
		if (maxError < 0.000001 || abs(maxError - lstError) < 0.00000001)
		{
			cout << i << endl;
			break;
		}
	}

	vector<double> inputs(1);
	vector<double> outputs;
	double error = 0.0;
	for (int i = 0; i < num_check; ++i)
	{
		inputs[0] = rand()*1.0 / RAND_MAX;
		double ideal = TESTFUNC(inputs[0]);
		NeuroNet::Matrix inpm;
		inpm.operator=(inputs);
		outputs = net.Run(inpm);
		double curerror = abs(outputs[0] - ideal);
		error += curerror;
		cout << inputs[0] << " " << outputs[0] << " (" << ideal << ") error: " << curerror << endl;
	}
	cout << endl << "Error: " << error << " sr error: " << error / num_check << endl;


	system("pause");
	return 0;
}