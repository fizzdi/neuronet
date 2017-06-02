#include "NeuralNetwork.h"
using namespace std;

int main()
{
	NeuroNet::NeuralNetwork net;
	srand(1);
	net.Init(2, 1, 10);
	net.ProblemSet.push_back(NeuroNet::Problem({ 1,1 }, { 0 }));
	net.ProblemSet.push_back(NeuroNet::Problem({ 0,1 }, { 1 }));
	net.ProblemSet.push_back(NeuroNet::Problem({ 1,0 }, { 1 }));
	net.ProblemSet.push_back(NeuroNet::Problem({ 0,0 }, { 0 }));
	net.RunProblemSet();
	net.CorrectWeights();
	cout << endl << endl << "CORRECT" << endl << endl;
	net.RunProblemSet();
	net.CorrectWeights();
	cout << endl << endl << "CORRECT" << endl << endl;
	net.RunProblemSet();
	net.CorrectWeights();
	cout << endl << endl << "CORRECT" << endl << endl;
	net.RunProblemSet();
	net.CorrectWeights();
	cout << endl << endl << "CORRECT" << endl << endl;
	net.RunProblemSet();

	system("pause");
	return 0;
}