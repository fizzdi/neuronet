#include "NeuralNetwork.h"
#include <algorithm>
#include <string>
using namespace NeuroNet;

NeuralNetwork::NeuralNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
{
	Layers.clear();
	Layers.emplace_back(Layer(InputCount, 0, LINE));
	Layers.emplace_back(Layer(NeuronCount, InputCount, HiddenLayerFunction));
	Layers.back().NguenWidrow(-2, 2, -1, 1);
	/*Layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
	Layers.back().NguenWidrow(-2, 2, -1, 1);*/
	/*Layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
	Layers.back().NguenWidrow(-2, 2, -1, 1);
	*/
	
	Layers.emplace_back(Layer(OutputCount, NeuronCount, LINE));
	//Layers.back().NguenWidrow(-1, 1, -1, 1);
	_countlayers = Layers.size();

	eligibility.resize(_countlayers);
	for (int i = 0; i < _countlayers; ++i)
	{
		eligibility[i] = Matrix2d(Layers[i].Grad.GetVerticalSize(), Layers[i].Grad.GetHorizontalSize());
		eligibility[i].Fill(0.0);
	}
}

double NeuralNetwork::RunTrainingSetOffline(std::deque<Problem> &TrainingSet, bool print)
{
	for (int i = 0; i < _countlayers; ++i)
	{
		Layers[i].GradSum.Fill(0.0);
		Layers[i].DeltaSum.Fill(0.0);
	}

	for (int t = 0; t < TrainingSet.size(); ++t)
	{
		Run(TrainingSet[t].inputs);
		CalcGradDelta(TrainingSet[t].outputs);
		for (int i = 1; i < _countlayers; ++i)
		{
			Layers[i].GradSum += Layers[i].Grad;
			Layers[i].DeltaSum += Layers[i].Delta;
		}
		if (print) PrintProblemResult(TrainingSet[t]);
	}

	if (print) std::cout << std::endl << "=======CORRECT==========" << std::endl;

	double normGrad = 0.0;
	for (int i = 0; i < _countlayers; ++i)
		normGrad += sqrt(Layers[i].GradSum * !Layers[i].GradSum).sum();

	if (normGrad < 1e-6)
		return 0.0;
	ResilientPropagationOffline();
	for (int i = 0; i < _countlayers; ++i)
	{
		Layers[i].LastGradSum = Layers[i].GradSum;
		Layers[i].LastDeltaSum = Layers[i].DeltaSum;
	}


	double error = 0.0;
	for each (Problem test in TrainingSet)
	{
		Run(test.inputs);
		if (print) PrintProblemResult(test);
		error += CalculateError(test, print);
	}
	return error;
}

void NeuralNetwork::PrintProblemResult(Problem & test)
{
	std::cout << *this;
	std::cout << "Expected results:" << std::endl;
	for (int i = 0; i < test.outputs.GetHorizontalSize(); ++i)
		std::cout << test.outputs.at(0, i) << " ";
	std::cout << std::endl;
}

void NeuralNetwork::ResilientPropagation()
{
	const double EttaPlus = 1.2, EttaMinus = 0.5;
	for (int i = 1; i < _countlayers; ++i)
	{
		for (int j = 0; j < Layers[i].Grad.GetVerticalSize(); ++j)
		{
			for (int k = 0; k < Layers[i].Grad.GetHorizontalSize(); ++k)
			{
				double cur_correct = 0.0;
				double curmatrixult = Layers[i].Grad.at(j, k) * Layers[i].LastGrad.at(j, k);
				if (curmatrixult == 0.0)
					cur_correct = getRand(0, 1);
				else if (curmatrixult > 0.0)
					cur_correct = std::min(EttaPlus * Layers[i].CorrectVal.at(j, k), 50.0);
				else if (curmatrixult < 0.0)
					cur_correct = std::max(EttaMinus * Layers[i].CorrectVal.at(j, k), 1e-6);

				Layers[i].CorrectVal.at(j, k) = cur_correct;

				if (Layers[i].Grad.at(j, k) == 0.0) continue;

				if (Layers[i].Grad.at(j, k) > 0)
					Layers[i].Weights.at(j, k) += -cur_correct;
				else
					Layers[i].Weights.at(j, k) += cur_correct;
			}
		}

		for (int j = 0; j < Layers[i].Delta.GetHorizontalSize(); ++j)
		{
			double cur_correct = 0.0;
			double curmatrixult = Layers[i].Delta.at(0, j) * Layers[i].LastDelta.at(0, j);
			if (curmatrixult == 0.0)
				cur_correct = getRand(0, 1);
			else if (curmatrixult > 0.0)
				cur_correct = std::min(EttaPlus * Layers[i].BiasCorrectVal.at(0, j), 50.0);
			else if (curmatrixult < 0.0)
				cur_correct = std::max(EttaMinus * Layers[i].BiasCorrectVal.at(0, j), 1e-6);

			Layers[i].BiasCorrectVal.at(0, j) = cur_correct;

			if (Layers[i].Delta.at(0, j) == 0.0) continue;

			if (Layers[i].Delta.at(0, j) > 0)
				Layers[i].Bias.at(0, j) += -cur_correct;
			else
				Layers[i].Bias.at(0, j) += cur_correct;
		}
	}
}

double NeuralNetwork::CalculateError(Problem & test, bool print)
{
	//MSE
	double error = (test.outputs - Layers.back().Axons).multiplication(test.outputs - Layers.back().Axons).sum() / 2.0;
	//RootMSE
	//double error = (test.outputs - Layers.back().Axons).abs().sqrt().sum() / 2.0;

	if (print)
	{
		std::cout << "MSE: " << error << std::endl;
		std::cout << "---------------------------------------------------------------------------------" << std::endl;
	}
	return error;
}

void NeuralNetwork::ResilientPropagationOffline()
{
	const double EttaPlus = 1.2, EttaMinus = 0.5;
	for (int i = 1; i < _countlayers; ++i)
	{
		for (int j = 0; j < Layers[i].GradSum.GetVerticalSize(); ++j)
		{
			for (int k = 0; k < Layers[i].GradSum.GetHorizontalSize(); ++k)
			{
				double cur_correct = 0.0;
				double curmatrixult = Layers[i].GradSum.at(j, k) * Layers[i].LastGradSum.at(j, k);
				if (curmatrixult == 0.0)
					cur_correct = getRand(0, 1);
				else if (curmatrixult > 0.0)
					cur_correct = std::min(EttaPlus * Layers[i].CorrectVal.at(j, k), 50.0);
				else if (curmatrixult < 0.0)
					cur_correct = std::max(EttaMinus * Layers[i].CorrectVal.at(j, k), 1e-6);

				Layers[i].CorrectVal.at(j, k) = cur_correct;

				if (Layers[i].GradSum.at(j, k) == 0.0) continue;

				if (Layers[i].GradSum.at(j, k) > 0)
					Layers[i].Weights.at(j, k) += -cur_correct;
				else
					Layers[i].Weights.at(j, k) += cur_correct;
			}
		}

		for (int j = 0; j < Layers[i].DeltaSum.GetHorizontalSize(); ++j)
		{
			double cur_correct = 0.0;
			double curmatrixult = Layers[i].DeltaSum.at(0, j) * Layers[i].LastDeltaSum.at(0, j);
			if (curmatrixult == 0.0)
				cur_correct = getRand(0, 1);
			else if (curmatrixult > 0.0)
				cur_correct = std::min(EttaPlus * Layers[i].BiasCorrectVal.at(0, j), 50.0);
			else if (curmatrixult < 0.0)
				cur_correct = std::max(EttaMinus * Layers[i].BiasCorrectVal.at(0, j), 1e-6);

			Layers[i].BiasCorrectVal.at(0, j) = cur_correct;

			if (Layers[i].DeltaSum.at(0, j) == 0.0) continue;

			if (Layers[i].DeltaSum.at(0, j) > 0)
				Layers[i].Bias.at(0, j) += -cur_correct;
			else
				Layers[i].Bias.at(0, j) += cur_correct;
		}
	}
}

double NeuroNet::NeuralNetwork::RMSTraining(std::deque<Problem> &TrainingSet)
{
	for (int i = 1; i < _countlayers; ++i)
	{
		Layers[i].GradSum.Fill(0.0);
		Layers[i].DeltaSum.Fill(0.0);
	}

	for (int t = 0; t < TEST_COUNT; ++t)
	{
		int test = myrand() % TrainingSet.size();
		Run(TrainingSet[test].inputs);
		//debug << "test.inputs: " << std::endl<< test.inputs << std::endl;
		CalcGradDelta(TrainingSet[test].outputs);
		for (int i = 1; i < _countlayers; ++i)
		{
			Layers[i].GradSum += Layers[i].Grad;
			Layers[i].DeltaSum += Layers[i].Delta;
		}
	}

	for (int i = 1; i < _countlayers; ++i)
	{
		//Layers[i].RMS.Fill(0.0);
	//	//debug << "RMS1: " << std::endl<<  Layers[i].RMS << std::endl;
		//debug.flush();
		Layers[i].RMS *= RMS_GAMMA;
		//	//debug << "RMS2: " << std::endl << Layers[i].RMS << std::endl;
		//debug.flush();
		Layers[i].RMS += Layers[i].GradSum.multiplication(Layers[i].GradSum) * (1.0 - RMS_GAMMA);
		//	//debug << "RMS3: " << std::endl << Layers[i].RMS << std::endl << std::endl;
		//debug.flush();
		Layers[i].RMSBias *= RMS_GAMMA;
		Layers[i].RMSBias += Layers[i].DeltaSum.multiplication(Layers[i].DeltaSum) * (1.0 - RMS_GAMMA);

		Layers[i].RMSN *= RMS_GAMMA;
		Layers[i].RMSN += Layers[i].GradSum * (1.0 - RMS_GAMMA);
		Layers[i].RMSNBias *= RMS_GAMMA;
		Layers[i].RMSNBias += Layers[i].DeltaSum * (1.0 - RMS_GAMMA);
	}


	double normGrad = 0.0;
	for (int i = 1; i < _countlayers; ++i)
		normGrad += sqrt(Layers[i].GradSum * !Layers[i].GradSum).sum();

	if (normGrad < 1e-6)
		return 0.0;

	RMSPropagation();

	double error = 0.0;
	for each (Problem test in TrainingSet)
	{
		Run(test.inputs);
		error += CalculateError(test);
	}
	return error;
}

void NeuroNet::NeuralNetwork::RMSPropagation()
{
	for (int i = 1; i < _countlayers; ++i)
	{
		Layers[i].Weights -= (RMS_LEARNRATE / sqrt(Layers[i].RMS - Layers[i].RMSN.multiplication(Layers[i].RMSN) + RMS_EPSILON)).multiplication(Layers[i].GradSum);
		Layers[i].Bias -= (RMS_LEARNRATE / sqrt(Layers[i].RMSBias - Layers[i].RMSNBias.multiplication(Layers[i].RMSNBias) + RMS_EPSILON)).multiplication(Layers[i].DeltaSum);
		//Layers[i].Weights -= (RMS_LEARNRATE / sqrt(Layers[i].RMS - Layers[i].RMSN.multiplication(Layers[i].RMSN) + RMS_EPSILON)).multiplication(Layers[i].GradSum);
		//Layers[i].Bias -= (RMS_LEARNRATE / sqrt(Layers[i].RMSBias - Layers[i].RMSNBias.multiplication(Layers[i].RMSNBias) + RMS_EPSILON)).multiplication(Layers[i].DeltaSum);
	}
}

void NeuralNetwork::CalcGradDelta(const double output)
{
	CalcGradDelta(std::vector<double>(1, output));
}

void NeuralNetwork::CalcGradDelta(const std::vector<double>& outputs)
{
	CalcGradDelta(Matrix2d(outputs));
}
void NeuralNetwork::CalcGradDelta(Matrix2d &outputs)
{
	for (int i = 1; i < _countlayers; ++i)
	{
		Layers[i].LastDelta = Layers[i].Delta;
		Layers[i].LastGrad = Layers[i].Grad;
	}

	Layers.back().Delta = (Layers.back().Axons - outputs).multiplication(Layers.back().GetDiff());
	////debug << "Layers.back().Delta: " << std::endl << Layers.back().Delta << std::endl;
	////debug << "Layers.back().Axons: " << std::endl << Layers.back().Axons << std::endl;
	////debug << "outputs: " << std::endl << outputs << std::endl;
	////debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;
	////debug << "Layers.back().GetDiff(): " << std::endl << Layers.back().GetDiff() << std::endl;
	////debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;

	////debug.flush();
	////debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;
	for (int i = _countlayers - 2; i >= 0; --i)
	{
		////debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;
		Layers[i].Delta = (Layers[i + 1].Delta * Layers[i + 1].Weights).multiplication(Layers[i].GetDiff());
		////debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;
		////debug << "Layers["<<i<<"].Delta: " << std::endl << Layers[i].Delta << std::endl;
		////debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;
	}
	for (int i = 1; i < _countlayers; ++i)
	{
		Layers[i].Grad = !Layers[i].Delta * Layers[i - 1].Axons;
		//debug << "Layers["<<i<<"].Grad: " << std::endl << Layers[i].Grad << std::endl;
	}
}

Matrix2d& NeuralNetwork::GetOutputGrad()
{
	return Layers.back().Grad;
}

void NeuralNetwork::MCQLCorrect()
{
	/*for (int i = 1; i < _countlayers; ++i)
	Layers[i].Weights += eligibility[i] * ALPHA*(r + GAMMA*Q - lastQ);

	std::vector<double> out;
	out.push_back(Q);
	CalcGradDelta(out);

	for (int i = 1; i < _countlayers; ++i)
	eligibility[i] = Layers[i].Grad + eligibility[i] * GAMMA*LAMBDA;
	}
	catch (std::exception ex)
	{
	//debug << "EXCEPTION: " << std::endl << ex.what() << std::endl;
	throw new std::exception("Exit");
	}*/
}

void NeuroNet::NeuralNetwork::debuginfo(std::ostream & debug) const
{
	//input neurons
	//hiddens neurons
	//outputs
	for (int i = 1; i < _countlayers; ++i)
	{
		//debug << Layers[i].Weights << std::endl;
		//debug << Layers[i].Bias << std::endl;
	}
	//debug << std::endl;
}


Matrix2d NeuralNetwork::GetOut() const
{
	return Layers.back().Axons;
}


void NeuralNetwork::Run()
{
	//init hidden and output layers
	for (int i = 1; i < (int)Layers.size(); ++i)
	{
		Layers[i].CalculateStates(Layers[i - 1]);
		Layers[i].CalculateAxons();
		//debug << "Layers[" << i << "].Axons: " << std::endl << Layers[i].Axons << std::endl;
	}
}

void NeuralNetwork::Run(const std::vector<double>& input)
{
	Run(Matrix2d(input));
}

void NeuralNetwork::Run(Matrix2d& input)
{
	//init input layer
	for (int i = 0; i < (int)input.GetHorizontalSize(); ++i)
		Layers[0].States.at(0, i) = input.at(0, i);
	//debug << "Layers[0].States: " << std::endl<< Layers[0].States << std::endl;
	Layers[0].CalculateAxons();
	//debug << "Layers[0].Axons: " << std::endl << Layers[0].Axons << std::endl;

	Run();
}
void NeuralNetwork::AddTest(std::deque<Problem> &TrainingSet, const std::vector<double> &input, const std::vector<double> &ideal) const
{
	if (TrainingSet.size() + 1 == MAX_TEST_COUNT)
		TrainingSet.pop_front();
	auto pr = Problem();
	pr.inputs = input;
	pr.outputs = ideal;
	TrainingSet.emplace_back(pr);
}
void NeuralNetwork::AddTest(std::deque<Problem> &TrainingSet, const Matrix2d& input, const Matrix2d& ideal) const
{
	if (TrainingSet.size() + 1 == MAX_TEST_COUNT)
		TrainingSet.pop_front();
	auto pr = Problem();
	pr.inputs = input;
	pr.outputs = ideal;
	TrainingSet.emplace_back(pr);
}

std::ostream & NeuroNet::operator<<(std::ostream & os, NeuralNetwork & net)
{
	for (int i = 1; i < (int)net.Layers.size(); ++i)
	{
		os << "Weights " << i << " -> " << i - 1 << ": " << std::endl;
		os << net.Layers[i].Weights << std::endl;
	}

	for (int i = 1; i < (int)net.Layers.size(); ++i)
	{
		os << "Bias " << i << " -> " << i - 1 << ": " << std::endl;
		os << net.Layers[i].Bias << std::endl;
	}

	os << std::endl << "Input neurons:" << std::endl;
	for (int i = 0; i < net.Layers[0].Axons.GetHorizontalSize(); ++i)
		os << net.Layers[0].Axons.at(0, i) << " ";

	os << std::endl;

	os << std::endl << "Output neurons:" << std::endl;
	for (int i = 0; i < net.Layers.back().Axons.GetHorizontalSize(); ++i)
		os << net.Layers.back().Axons.at(0, i) << " ";
	os << std::endl;
	os << "===============" << std::endl;
	return os;
}
