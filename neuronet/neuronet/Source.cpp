#include "NeuralNetwork.h"
using namespace std;
const int num_check = 10;
const int HiddenNeuron = 30; //sin - 30
const int SampleCount = 10;
const int Epoh = 100;
#define TESTFUNC cos
int main()
{
	std::ios::sync_with_stdio(false);
	NeuroNet::NeuralNetwork net;
	srand(17);
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
	//	cout << i << ") " << maxError << endl;
		if (maxError < 0.001 || abs(maxError - lstError) < 0.0001)
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
		outputs = net.Run(inputs);
		double curerror = abs(outputs[0] - ideal);
		error += curerror;
		cout << inputs[0] << " " << outputs[0] << " (" << ideal << ") error: " << curerror << endl;
	}
	cout << endl << "Error: " << error << " sr error: " << error / num_check << endl;


	system("pause");
	return 0;
}