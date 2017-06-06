#include "NeuralNetwork.h"
using namespace std;
const int num_check = 10;
const int HiddenNeuron = 7;
const int SampleCount = 10;
const int Epoh = 10;

int main()
{
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
		double y = sin(x);
		net.TrainingSet.push_back(NeuroNet::Problem({ x }, { y }));
	}
	net.RunTrainingSet(true);
	for (int i = 0; i < Epoh; ++i)
	{
		net.CorrectWeights();
		//cout << endl << endl << "CORRECT" << endl << endl;
		net.RunTrainingSet(true);
	}

	vector<double> inputs(1);
	vector<double> outputs;
	double error = 0.0;
	for (int i = 0; i < num_check; ++i)
	{
		inputs[0] = rand()*1.0 / RAND_MAX;
		double ideal = sin(inputs[0]);
		outputs = net.Run(inputs);
		error += abs(outputs[0] - ideal);
		cout << inputs[0] << " " << outputs[0] << " (" << ideal << ") " << endl;
	}
	error /= num_check;
	cout << endl << "Error: " << error << endl;


	system("pause");
	return 0;
}