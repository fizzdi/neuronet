#include "NeuralNetwork.h"
#include <algorithm>
#include <string>
using namespace NeuroNet;

NeuralNetwork::NeuralNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
{
	layers.clear();
	layers.emplace_back(Layer(InputCount, 0, LINE));
	layers.emplace_back(Layer(NeuronCount, InputCount, HiddenLayerFunction));
	layers.back().NguenWidrow(-2, 2, -1, 1);
	/*layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
	layers.back().NguenWidrow(-2, 2, -1, 1);*/
	/*layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
	layers.back().NguenWidrow(-2, 2, -1, 1);
	*/
	
	layers.emplace_back(Layer(OutputCount, NeuronCount, LINE));
	//layers.back().NguenWidrow(-1, 1, -1, 1);
	countLayers = (int)layers.size();
}

double NeuralNetwork::RunTrainingSetOffline(training_set &TrainingSet)
{
	for (int i = 0; i < countLayers; ++i)
	{
		layers[i].GradSum.Fill(0.0);
		layers[i].DeltaSum.Fill(0.0);
	}

	for (int t = 0; t < TrainingSet.size(); ++t)
	{
		Run(TrainingSet[t].inputs);
		CalcGradDelta(TrainingSet[t].outputs);
		for (int i = 1; i < countLayers; ++i)
		{
			layers[i].GradSum += layers[i].Grad;
			layers[i].DeltaSum += layers[i].Delta;
		}
	}

	double normGrad = 0.0;
	for (int i = 0; i < countLayers; ++i)
		normGrad += sqrt(layers[i].GradSum * !layers[i].GradSum).sum();

	if (normGrad < 1e-6)
		return 0.0;
	ResilientPropagationOffline();
	for (int i = 0; i < countLayers; ++i)
	{
		layers[i].LastGradSum = layers[i].GradSum;
		layers[i].LastDeltaSum = layers[i].DeltaSum;
	}


	double error = 0.0;
	for each (Problem test in TrainingSet)
	{
		Run(test.inputs);
		error += CalculateError(test);
	}
	return error;
}

void NeuralNetwork::ResilientPropagation()
{
	const double EttaPlus = 1.2, EttaMinus = 0.5;
	for (int i = 1; i < countLayers; ++i)
	{
		for (int j = 0; j < layers[i].Grad.GetVerticalSize(); ++j)
		{
			for (int k = 0; k < layers[i].Grad.GetHorizontalSize(); ++k)
			{
				double cur_correct = 0.0;
				double curmatrixult = layers[i].Grad.at(j, k) * layers[i].LastGrad.at(j, k);
				if (curmatrixult == 0.0)
					cur_correct = getRand(0, 1);
				else if (curmatrixult > 0.0)
					cur_correct = std::min(EttaPlus * layers[i].CorrectVal.at(j, k), 50.0);
				else if (curmatrixult < 0.0)
					cur_correct = std::max(EttaMinus * layers[i].CorrectVal.at(j, k), 1e-6);

				layers[i].CorrectVal.at(j, k) = cur_correct;

				if (layers[i].Grad.at(j, k) == 0.0) continue;

				if (layers[i].Grad.at(j, k) > 0)
					layers[i].Weights.at(j, k) += -cur_correct;
				else
					layers[i].Weights.at(j, k) += cur_correct;
			}
		}

		for (int j = 0; j < layers[i].Delta.GetHorizontalSize(); ++j)
		{
			double cur_correct = 0.0;
			double curmatrixult = layers[i].Delta.at(0, j) * layers[i].LastDelta.at(0, j);
			if (curmatrixult == 0.0)
				cur_correct = getRand(0, 1);
			else if (curmatrixult > 0.0)
				cur_correct = std::min(EttaPlus * layers[i].BiasCorrectVal.at(0, j), 50.0);
			else if (curmatrixult < 0.0)
				cur_correct = std::max(EttaMinus * layers[i].BiasCorrectVal.at(0, j), 1e-6);

			layers[i].BiasCorrectVal.at(0, j) = cur_correct;

			if (layers[i].Delta.at(0, j) == 0.0) continue;

			if (layers[i].Delta.at(0, j) > 0)
				layers[i].Bias.at(0, j) += -cur_correct;
			else
				layers[i].Bias.at(0, j) += cur_correct;
		}
	}
}

double NeuralNetwork::CalculateError(Problem & test, bool print)
{
	//MSE
	double error = (test.outputs - layers.back().Axons).multiplication(test.outputs - layers.back().Axons).sum() / 2.0;
	//RootMSE
	//double error = (test.outputs - layers.back().Axons).abs().sqrt().sum() / 2.0;

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
	for (int i = 1; i < countLayers; ++i)
	{
		for (int j = 0; j < layers[i].GradSum.GetVerticalSize(); ++j)
		{
			for (int k = 0; k < layers[i].GradSum.GetHorizontalSize(); ++k)
			{
				double cur_correct = 0.0;
				double curmatrixult = layers[i].GradSum.at(j, k) * layers[i].LastGradSum.at(j, k);
				if (curmatrixult == 0.0)
					cur_correct = getRand(0, 1);
				else if (curmatrixult > 0.0)
					cur_correct = std::min(EttaPlus * layers[i].CorrectVal.at(j, k), 50.0);
				else if (curmatrixult < 0.0)
					cur_correct = std::max(EttaMinus * layers[i].CorrectVal.at(j, k), 1e-6);

				layers[i].CorrectVal.at(j, k) = cur_correct;

				if (layers[i].GradSum.at(j, k) == 0.0) continue;

				if (layers[i].GradSum.at(j, k) > 0)
					layers[i].Weights.at(j, k) += -cur_correct;
				else
					layers[i].Weights.at(j, k) += cur_correct;
			}
		}

		for (int j = 0; j < layers[i].DeltaSum.GetHorizontalSize(); ++j)
		{
			double cur_correct = 0.0;
			double curmatrixult = layers[i].DeltaSum.at(0, j) * layers[i].LastDeltaSum.at(0, j);
			if (curmatrixult == 0.0)
				cur_correct = getRand(0, 1);
			else if (curmatrixult > 0.0)
				cur_correct = std::min(EttaPlus * layers[i].BiasCorrectVal.at(0, j), 50.0);
			else if (curmatrixult < 0.0)
				cur_correct = std::max(EttaMinus * layers[i].BiasCorrectVal.at(0, j), 1e-6);

			layers[i].BiasCorrectVal.at(0, j) = cur_correct;

			if (layers[i].DeltaSum.at(0, j) == 0.0) continue;

			if (layers[i].DeltaSum.at(0, j) > 0)
				layers[i].Bias.at(0, j) += -cur_correct;
			else
				layers[i].Bias.at(0, j) += cur_correct;
		}
	}
}

double NeuroNet::NeuralNetwork::RMSTraining(training_set &TrainingSet)
{
	for (int i = 1; i < countLayers; ++i)
	{
		layers[i].GradSum.Fill(0.0);
		layers[i].DeltaSum.Fill(0.0);
	}

	for (int t = 0; t < TEST_COUNT; ++t)
	{
		int test = myrand() % TrainingSet.size();
		Run(TrainingSet[test].inputs);
		CalcGradDelta(TrainingSet[test].outputs);
		for (int i = 1; i < countLayers; ++i)
		{
			layers[i].GradSum += layers[i].Grad;
			layers[i].DeltaSum += layers[i].Delta;
		}
	}

	for (int i = 1; i < countLayers; ++i)
	{
		layers[i].RMS *= RMS_GAMMA;
		layers[i].RMS += layers[i].GradSum.multiplication(layers[i].GradSum) * (1.0 - RMS_GAMMA);
		layers[i].RMSBias *= RMS_GAMMA;
		layers[i].RMSBias += layers[i].DeltaSum.multiplication(layers[i].DeltaSum) * (1.0 - RMS_GAMMA);

		layers[i].RMSN *= RMS_GAMMA;
		layers[i].RMSN += layers[i].GradSum * (1.0 - RMS_GAMMA);
		layers[i].RMSNBias *= RMS_GAMMA;
		layers[i].RMSNBias += layers[i].DeltaSum * (1.0 - RMS_GAMMA);
	}


	double normGrad = 0.0;
	for (int i = 1; i < countLayers; ++i)
		normGrad += sqrt(layers[i].GradSum * !layers[i].GradSum).sum();

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
	for (int i = 1; i < countLayers; ++i)
	{
		layers[i].Weights -= (RMS_LEARNRATE / sqrt(layers[i].RMS - layers[i].RMSN.multiplication(layers[i].RMSN) + RMS_EPSILON)).multiplication(layers[i].GradSum);
		layers[i].Bias -= (RMS_LEARNRATE / sqrt(layers[i].RMSBias - layers[i].RMSNBias.multiplication(layers[i].RMSNBias) + RMS_EPSILON)).multiplication(layers[i].DeltaSum);
		//layers[i].Weights -= (RMS_LEARNRATE / sqrt(layers[i].RMS - layers[i].RMSN.multiplication(layers[i].RMSN) + RMS_EPSILON)).multiplication(layers[i].GradSum);
		//layers[i].Bias -= (RMS_LEARNRATE / sqrt(layers[i].RMSBias - layers[i].RMSNBias.multiplication(layers[i].RMSNBias) + RMS_EPSILON)).multiplication(layers[i].DeltaSum);
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
void NeuralNetwork::CalcGradDelta(const Matrix2d &outputs)
{
	for (int i = 1; i < countLayers; ++i)
	{
		layers[i].LastDelta = layers[i].Delta;
		layers[i].LastGrad = layers[i].Grad;
	}

	layers.back().Delta = (layers.back().Axons - outputs).multiplication(layers.back().GetDiff());
	for (int i = countLayers - 2; i >= 0; --i)
		layers[i].Delta = (layers[i + 1].Delta * layers[i + 1].Weights).multiplication(layers[i].GetDiff());
	for (int i = 1; i < countLayers; ++i)
		layers[i].Grad = !layers[i].Delta * layers[i - 1].Axons;
}

void NeuroNet::NeuralNetwork::PrintFullInfo(std::ostream & debug) const
{
	debug << typeid(*this).name() << std::endl;
	//Neural config
	debug << SENSOR_COUNT << std::endl;
	debug << INPUT_NEURON_COUNT << std::endl;
	debug << HIDDEN_NEURON_COUNT << std::endl;
	debug << OUTPUT_NEURON_COUNT << std::endl;
	debug << TEST_COUNT << std::endl;
	debug << MAX_TEST_COUNT << std::endl;
	debug << TRAIN_EPOCH << std::endl;
	debug << TRAIN_PERIOD << std::endl;
	debug << TRAIN_EPS << std::endl;

	//RMS config
	debug << RMS_GAMMA << std::endl;
	debug << RMS_LEARNRATE << std::endl;
	debug << RMS_EPSILON << std::endl;

	for (int i = 1; i < countLayers; ++i)
		debug << layers[i].Weights << std::endl;

	for (int i = 1; i < countLayers; ++i)
		debug << layers[i].Bias << std::endl;
}


Matrix2d NeuralNetwork::GetOut() const
{
	return layers.back().Axons;
}

void NeuroNet::NeuralNetwork::SetInput(const Matrix2d & inputs)
{
	for (int i = 0; i < inputNeuron; ++i)
		layers[0].States.at(0, i);
}


void NeuralNetwork::Run()
{
	//init hidden and output layers
	for (int i = 1; i < countLayers; ++i)
	{
		layers[i].CalculateStates(layers[i - 1]);
		layers[i].CalculateAxons();
	}
}

void NeuralNetwork::Run(const std::vector<double>& input)
{
	Run(Matrix2d(input));
}

void NeuralNetwork::Run(const Matrix2d& input)
{
	//init input layer
	SetInput(input);
	layers[0].CalculateAxons();

	Run();
}

void NeuroNet::AddTest(training_set & TrainingSet, const std::vector<double>& input, const std::vector<double>& ideal)
{
	if (TrainingSet.size() + 1 == MAX_TEST_COUNT)
		TrainingSet.pop_front();
	auto pr = Problem();
	pr.inputs = input;
	pr.outputs = ideal;
	TrainingSet.emplace_back(pr);
}

void NeuroNet::AddTest(training_set & TrainingSet, const Matrix2d & input, const Matrix2d & ideal)
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
	for (int i = 1; i < (int)net.countLayers; ++i)
	{
		os << "Weights " << i << " -> " << i - 1 << ": " << std::endl;
		os << net.layers[i].Weights << std::endl;
	}

	for (int i = 1; i < (int)net.countLayers; ++i)
	{
		os << "Bias " << i << " -> " << i - 1 << ": " << std::endl;
		os << net.layers[i].Bias << std::endl;
	}

	os << std::endl << "Input neurons:" << std::endl;
	for (int i = 0; i < net.layers[0].Axons.GetHorizontalSize(); ++i)
		os << net.layers[0].Axons.at(0, i) << " ";

	os << std::endl;

	os << std::endl << "Output neurons:" << std::endl;
	for (int i = 0; i < net.layers.back().Axons.GetHorizontalSize(); ++i)
		os << net.layers.back().Axons.at(0, i) << " ";
	os << std::endl;
	os << "===============" << std::endl;
	return os;
}
