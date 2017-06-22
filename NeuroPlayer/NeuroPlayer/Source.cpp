#define _CRT_SECURE_NO_WARNINGS
#include "NeuralNetwork.h"
#include "ElmanNetwork.h"
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <fstream>
using namespace std;
const int TestCount = 10000;
const int HiddenNeuron = 500; //sin - 11 15
const int SampleCount = 100;
const int MaxEpoch = 10000;
const double MinMSE = 1e-7;
const double MaxMSE = 1e12;
const double MinMSEIncrement = 1e-9;

double testfun(double x)
{
	return x + 1;
#include "ElmanNetwork.h"
}

const double xmin = -1;
const double xmax = 1;

	ofstream debug("tdebut.txt");

deque<NeuroNet::Problem> TrainingSet;

int main()
{
	NeuroNet::Matrix2d a(3, 3);
	NeuroNet::Matrix2d b(3, 3);
	a.at(0, 0) = 1;
	a.at(0, 1) = 2;
	a.at(0, 2) = 3;
	a.at(1, 0) = 4;
	a.at(1, 1) = 5;
	a.at(1, 2) = 6;
	a.at(2, 0) = 7;
	a.at(2, 1) = 8;
	a.at(2, 2) = 9;
	b.at(0, 0) = 1;
	b.at(0, 1) = 2;
	b.at(0, 2) = 3;
	b.at(1, 0) = 4;
	b.at(1, 1) = 5;
	b.at(1, 2) = 6;
	b.at(2, 0) = 7;
	b.at(2, 1) = 8;
	b.at(2, 2) = 9;
	cout << a << endl << b << endl << endl;
	cout << a << endl << b << endl << 10 / a;
	cout << a << endl << b << endl << 10 / a;
	system("pause");
	return 0;


	std::ios::sync_with_stdio(false);
	freopen("debug.txt", "w+", stdout);
	NeuroNet::NeuralNetwork net(1,1,2,NeuroNet::AFType::TANH);
	//srand(1);
	srand((unsigned)time(NULL));

	bool failed = false;
	for (int i = 1; i <= SampleCount; ++i)
	{
	//double x = NeuroNet::Common::getRand(0, 15, true);
	//double x =rand()%10;
	double x = i;
	double y = testfun(x);
	TrainingSet.emplace_back(NeuroNet::Problem({ x }, { y }));
	}

	cout.flush();
	double lstError = 0.0;
	double Error = 0.0;
	int Epoch = 1;
	do
	{
		lstError = Error;
		bool fl = Epoch % 1 == 0;
		//Error = net.RunTrainingSet();
		Error = net.RMSTraining(TrainingSet);
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
		vector<double> inputs(1);
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
		net.Run(inputs);
		auto outputs = net.GetOut();
		double curerror = abs(outputs.at(0,0) - ideal) * abs(outputs.at(0, 0) - ideal) / 2;
		mse += curerror;
		max_error = max(max_error, abs(outputs.at(0, 0) - ideal));
		cout << setw(15) << inputs[0] << " |" << setw(15) << outputs.at(0, 0) << " |" << setw(15) << ideal << " |" << setw(15) << abs(outputs.at(0, 0) - ideal) << " |" << setw(15) << curerror << endl;
		max_mse = max(max_mse, curerror);


		}
		cout << endl << "Sum mse: " << mse << endl;
		cout << "MAX mse: " << max_mse << endl;
		cout << "Avg mse: " << mse / TestCount << endl;
		cout << "MAX error: " << max_error << endl;
		}
		
/*
		cout << fixed << setprecision(6);
		vector<double> inputs(10), outputs(10);
		for (int test = 0; test < TestCount; ++test)
		{
			int x = 0;
			cout << "(0) ";
			net.Run(TrainingSet[0].inputs);
			bool successful = false;
			bool stop = false;
			int predicted = -1;
			int index = 0;
			do
			{
				auto res = net.GetOut();
				for (int i = 0; i < 10; ++i)
				{
					cout << setw(10) << res.at(0, i) << " ";
					if (res.at(0, i) >= 0.4)
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
		*/
	}

	//system("pause");
	return 0;
}