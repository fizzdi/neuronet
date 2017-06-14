#include "NeuralNetwork.h"
#include "Elman.h"
#include <algorithm>
#include <ctime>
const int num_check = 50;
const int HiddenNeuron = 7; //sin - 30
const int SampleCount = 1;
const int Epoh = 20;

using namespace std;
double test_fun(double x)
{
	return sin(x);
}

int main()
{
	std::ios::sync_with_stdio(false);
	NeuroNet::Elman net;
	//srand(1);	
	srand(time(NULL));

	net.Init(6, 6, HiddenNeuron, NeuroNet::SIGM);
	/*for (int i = 0; i < SampleCount; ++i)
	{
		double x = -1.0+(rand()*1.0 / (RAND_MAX/2));
		double y = test_fun(x);*/
		//0 3 5 2 0
	net.TrainingSet.push_back(NeuroNet::Problem({ 1,0,0,0,0,0 }, { 0,0,0,1,0,0 }));
	net.TrainingSet.push_back(NeuroNet::Problem({ 0,0,0,1,0,0 }, { 0,0,0,0,0,1 }));
	net.TrainingSet.push_back(NeuroNet::Problem({ 0,0,0,0,0,1 }, { 0,0,1,0,0,0 }));
	net.TrainingSet.push_back(NeuroNet::Problem({ 0,0,1,0,0,0 }, { 1,0,0,0,0,0 }));
	//}
	cout << "0) " << net.RunTrainingSet(true) << endl;
	double lstError = 0.0;
	double maxError = 0.0;
	for (int i = 1; i <= Epoh; ++i)
	{
		cout << "=======================================================================" << endl;
		lstError = maxError;
		//net.CorrectWeights();
		//cout << endl << endl << "CORRECT" << endl << endl;
		maxError = net.RunTrainingSet(true);
		//if (i % 1000 == 0)
		cout << i << ") " << maxError << " (" << abs(maxError - lstError) << ")" << endl;
		cout.flush();
		int t = 1e6;
		if (maxError < 1e-3 || abs(maxError - lstError) < 1e-8 || maxError > 1e6)
		{
			cout << "STOP " << i << ") " << maxError << " (" << abs(maxError - lstError) << ")" << endl;
			break;
		}
	}


	cout << endl << "------------------------------------------" << endl;

	{
		double error = 0.0;
		double max_error = 0.0;
		NeuroNet::Matrix2d inp(1, 6);
		for (int t = 0; t < num_check; ++t)
		{
			int num = rand() % 6;
			auto outputs = net.Run(net.TrainingSet[0].inputs);
			cout << "(0)";
			int p;
			int index = 0;
			bool flag = false;
			bool stop = false;
			bool sucs = false;
			do
			{
				for (int j = 0; j < 6; ++j)
				{
					cout << outputs[0][j] << " ";
					if (outputs[0][j] > 0.4)
						p = j;
				}
				cout << endl;
				num = rand() % 6; //Enter a random letter.

				index += 1;              //Increment to the next position.
				if (index == 5) {
					flag = true;
				}
				else {
					cout << "(" << num << ") ";
				}

				for (int i = 0; i <= 5; i++) {
					if (i == num) {
						inp[0][i] = 1.0;
						if (i == p) {
							sucs = true;
						}
						else {
							//Failure. Stop this sample and try a new sample.
							stop = true;
						}
					}
					else {
						inp[0][i] = 0.0;
					}
				} // i

				outputs = net.Run(inp);

			} while (stop == false);
			if ((index > 4) && (sucs == true)) {
				//If the random sequence happens to be in the correct order,
				//the network reports success.
				cout << "Success." << endl;
				cout << "Completed " << t << " tests." << endl;
				stop = true;
				break;
			}
			else {
				cout << "Failed." << endl;
			}
		}
		cout << endl << "Sum Error: " << error << endl;
		cout << "MAX Error: " << max_error << endl;
		cout << "Avg Error: " << error / num_check << endl;

	}

	system("pause");
	return 0;
}