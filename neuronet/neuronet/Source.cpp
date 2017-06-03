#include "NeuralNetwork.h"
using namespace std;

int main()
{
	NeuroNet::NeuralNetwork net;
	srand(17);
	/*net.Init(2, 1, 2);
	net.ProblemSet.push_back(NeuroNet::Problem({ 1,1 }, { 0 }));
	net.ProblemSet.push_back(NeuroNet::Problem({ 0,0 }, { 0 }));
	net.ProblemSet.push_back(NeuroNet::Problem({ 0,1 }, { 1 }));
	net.ProblemSet.push_back(NeuroNet::Problem({ 1,0 }, { 1 }));
	*/

	net.Init(1, 1, 7);
	for (int i = 0; i < 10; ++i)
	{
		double x = rand()*1.0 / RAND_MAX;
		double y = sin(x);
		net.ProblemSet.push_back(NeuroNet::Problem({ x }, { y }));
	}
	net.RunProblemSet();
	int epoh = 10;
	while (epoh--)
	{
		net.CorrectWeights();
		cout << endl << endl << "CORRECT" << endl << endl;
		net.RunProblemSet();
	}

	vector<double> inputs(1);
	vector<double> outputs;
	for (int i = 0; i < 150; ++i)
	{
		inputs[0] = rand()*1.0 / RAND_MAX;
		outputs = net.Run(inputs);
		cout << inputs[0] << " " << outputs[0] << " (" << sin(inputs[0]) << ") " << outputs[0] - sin(inputs[0]) << endl;
	}
	


	system("pause");
	return 0;
}