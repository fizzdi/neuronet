#define _CRT_SECURE_NO_WARNINGS
#include "NeuralNetwork.h"
#include "ElmanNetwork.h"
#include "Common.h"
#include <algorithm>
#include <ctime>
#include <iomanip>
using namespace std;
const int TestCount = 10000;
const int HiddenNeuron = 500; //sin - 11 15
const int SampleCount = 10;
const int MaxEpoch = 10000;
const double MinMSE = 1e-7;
const double MaxMSE = 1e12;
const double MinMSEIncrement = 1e-9;

double testfun(double x)
{
	return x + 1;
}

const double xmin = -1;
const double xmax = 1;

int main()
{
	std::ios::sync_with_stdio(false);
	freopen("debug.txt", "w+", stdout);
	NeuroNet::ElmanNetwork net;
	//srand(1);
	srand((unsigned)time(NULL));

	bool failed = false;
	net.Init(10, 10, HiddenNeuron, NeuroNet::SIGM, NeuroNet::RPROP);
	/*for (int i = 1; i <= SampleCount; ++i)
	{
		//double x = NeuroNet::Common::getRand(0, 15, true);
		//double x =rand()%10;
		double x = i;
		double y = testfun(x);
		net.TrainingSet.push_back(NeuroNet::Problem({ x }, { y }));
	}*/

	cout.flush();
	double lstError = 0.0;
	double Error = 0.0;
	int Epoch = 1;
	do
	{
		lstError = Error;
		bool fl = Epoch % 1000 == 0;
		//Error = net.RunTrainingSet();
		Error = net.RunTrainingSetOffline(0);
		if (fl)
		{
			cout << endl << "================================================================================================" << endl;
			cout << setw(15) << Epoch << ") " << setw(15) << Error << " (" << setw(15) << abs(Error - lstError) << ")" << endl;
			cout.flush();
		}
		if (Error < MinMSE || abs(Error - lstError) < MinMSEIncrement)
		{
			cout << "STOP Epoch: " << setw(15) << Epoch << " " << setw(15) << Error << " (" << setw(15) << abs(Error - lstError) << ")" << endl;
			break;
		}

		if (Error > MaxMSE)
		{
			cout << "Failed" << endl;
			failed = true;
			break;
		}
		Epoch++;
	} while (Epoch < MaxEpoch);



	if (!failed)
	{
		cout << "End Epoch: " << setw(15) << Epoch << " Error:" << setw(15) << Error << " (" << setw(15) << abs(Error - lstError) << ")" << endl;
		cout << endl << "================================================================================================" << endl;
		cout << "Tests:" << endl;
		cout << setw(15) << "Input" << " |" << setw(15) << "Output" << " |" << setw(15) << "Answer" << " |" << setw(15) << "Different" << " |" << setw(15) << "MSE" << endl;
		/*vector<double> inputs(1);
		{
			double mse = 0.0;
			double max_mse = 0.0;
			double max_error = 0.0;
			for (int i = 0; i <= TestCount; ++i)
			{
				inputs[0] = i;
				//inputs[0] = NeuroNet::Common::getRand(0, 15, true);
				//inputs[0] = rand() % 100;

				double ideal = testfun(inputs[0]);
				net.Run(NeuroNet::Problem({ inputs[0] }, { testfun(inputs[0]) }));
				auto outputs = net.GetOut();
				double curerror = abs(outputs(0,0) - ideal) * abs(outputs(0, 0) - ideal) / 2;
				mse += curerror;
				max_error = max(max_error, abs(outputs(0, 0) - ideal));
				cout << setw(15) << inputs[0] << " |" << setw(15) << outputs(0, 0) << " |" << setw(15) << ideal << " |" << setw(15) << abs(outputs(0, 0) - ideal) << " |" << setw(15) << curerror << endl;
				max_mse = max(max_mse, curerror);


			}
			cout << endl << "Sum mse: " << mse << endl;
			cout << "MAX mse: " << max_mse << endl;
			cout << "Avg mse: " << mse / TestCount << endl;
			cout << "MAX error: " << max_error << endl;
		}
*/

		cout << fixed << setprecision(6);
		vector<double> inputs(10), outputs(10);
		for (int test = 0; test < TestCount; ++test)
		{
			int x = 0;
			cout << "(0) ";
			net.Run(net.TrainingSet[0]);
			bool successful = false;
			bool stop = false;
			int predicted = -1;
			int index = 0;
			do
			{
				auto res = net.GetOut();
				for (int i = 0; i < 10; ++i)
				{
					cout << setw(10) << res(0, i) << " ";
					if (res(0, i) >= 0.4)
						predicted = i;
				}
				cout << " <- " << predicted << endl;

				x = rand() % 10;
				index++;
				if (index == 10)
					stop = true;
				else cout << "(" << x << ") ";

				for (int i = 0; i < 10; ++i)
				{
					if (i == x) {
						inputs[i] = 1.0;
						if (i == predicted) {
							successful = true;
						}
						else {
							//Failure. Stop this sample and try a new sample.
							stop = true;
						}
					}
					else {
						inputs[i] = 0.0;
					}
				}
				net.Run(NeuroNet::Problem(inputs, outputs));
			} while (!stop);

			if ((index > 9) && (successful == true)) {
				cout << "Success." << endl;
				cout << "Completed " << test << " tests." << endl;
				stop = true;
				break;
			}
			else {
				cout << "Failed." << endl << endl;
				if (test > TestCount) {
					stop = true;
					cout << "Completed " << test << " tests with no success." << endl;
					break;
				}
			}
		}

	}

	//system("pause");
	return 0;
}