#include "MLP.h"
#include <queue>
using namespace NeuroNet;

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

#pragma region Protected Methods
void NeuroNet::MLP::Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
{
	layers.clear();
	layers.emplace_back(Layer(InputCount, 0, LINE));
	for (int i = 0; i < HIDDEN_LAYER_COUNT; ++i)
	{
		layers.emplace_back(Layer(NeuronCount, InputCount, HiddenLayerFunction));
		layers.back().NguenWidrow(-2, 2, -1, 1);
	}
	layers.emplace_back(Layer(OutputCount, NeuronCount, LINE));
	countLayers = (int)layers.size();
}

void MLP::Run()
{
	//init hidden and output layers
	for (int i = 1; i < countLayers; ++i)
	{
		layers[i].CalculateStates(layers[i - 1]);
		layers[i].CalculateAxons();
	}
}

void MLP::ResilientPropagation()
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

double NeuroNet::MLP::RMSTraining(training_set &TrainingSet)
{
	for (int i = 1; i < countLayers; ++i)
	{
		layers[i].GradSum.Fill(0.0);
		layers[i].DeltaSum.Fill(0.0);
	}

	std::queue<int> tests;
	if (TrainingSet.size() < TEST_COUNT)
	{
		for (int test = 0; test < TrainingSet.size(); ++test)
		{
			tests.push(test);
			Run(TrainingSet[test].inputs);
			CalcGradDelta(TrainingSet[test].outputs);
			for (int i = 1; i < countLayers; ++i)
			{
				layers[i].GradSum += layers[i].Grad;
				layers[i].DeltaSum += layers[i].Delta;
			}
		}
	}
	else
	{
		for (int t = 0; t < TEST_COUNT; ++t)
		{
			int test = myrand() % TrainingSet.size();
			tests.push(test);
			Run(TrainingSet[test].inputs);
			CalcGradDelta(TrainingSet[test].outputs);
			for (int i = 1; i < countLayers; ++i)
			{
				layers[i].GradSum += layers[i].Grad;
				layers[i].DeltaSum += layers[i].Delta;
			}
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
	while (!tests.empty())
	{
		Run(TrainingSet[tests.front()].inputs);
		error += CalculateError(TrainingSet[tests.front()]);
		tests.pop();
	}
	return error;
}

double MLP::CalculateError(Problem & test, bool print)
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

void MLP::CalcGradDelta(const double output)
{
	CalcGradDelta(std::vector<double>(1, output));
}

void MLP::CalcGradDelta(const std::vector<double>& outputs)
{
	CalcGradDelta(Matrix2d(outputs));
}

void MLP::CalcGradDelta(const Matrix2d &outputs)
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
#pragma endregion

#pragma region Public methods
MLP::MLP(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
{
	Init(InputCount, OutputCount, NeuronCount, HiddenLayerFunction);
}

NeuroNet::MLP::MLP(const std::string& FileName)
{
	std::ifstream fi(FileName);
	std::string str;
	std::getline(fi, str);
	int fun_act;
	//Neural config
	fi >> SENSOR_COUNT;
	fi >> INPUT_NEURON_COUNT;
	fi >> HIDDEN_NEURON_COUNT;
	fi >> HIDDEN_LAYER_COUNT;
	fi >> OUTPUT_NEURON_COUNT;
	fi >> TEST_COUNT;
	fi >> MAX_TEST_COUNT;
	fi >> TRAIN_EPOCH;
	fi >> TRAIN_PERIOD;
	fi >> TRAIN_EPS;
	fi >> END_TRAIN_TICK;
	fi >> fun_act;
	FUN_ACT = (AFType)fun_act;

	//RMS config
	fi >> RMS_GAMMA;
	fi >> RMS_LEARNRATE;
	fi >> RMS_EPSILON;

	for (int i = 1; i < countLayers; ++i)
		fi >> layers[i].Weights;

	for (int i = 1; i < countLayers; ++i)
		fi >> layers[i].Bias;
	Init(INPUT_NEURON_COUNT, OUTPUT_NEURON_COUNT, HIDDEN_NEURON_COUNT, NeuroNet::AFType::TANH);
}

void MLP::Run(const std::vector<double>& input)
{
	Run(Matrix2d(input));
}

void MLP::Run(const Matrix2d& input)
{
	//init input layer
	SetInput(input);
	layers[0].CalculateAxons();
	Run();
}

double MLP::RPROPTraining(training_set &TrainingSet)
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
	ResilientPropagation();
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

void NeuroNet::MLP::RMSPropagation()
{
	for (int i = 1; i < countLayers; ++i)
	{
		layers[i].Weights -= (RMS_LEARNRATE / sqrt(layers[i].RMS - layers[i].RMSN.multiplication(layers[i].RMSN) + RMS_EPSILON)).multiplication(layers[i].GradSum);
		layers[i].Bias -= (RMS_LEARNRATE / sqrt(layers[i].RMSBias - layers[i].RMSNBias.multiplication(layers[i].RMSNBias) + RMS_EPSILON)).multiplication(layers[i].DeltaSum);
	}
}

void NeuroNet::MLP::PrintFullInfo(std::ostream & outstream) const
{
	outstream << typeid(*this).name() << std::endl;
	//Neural config
	outstream << SENSOR_COUNT << std::endl;
	outstream << INPUT_NEURON_COUNT << std::endl;
	outstream << HIDDEN_NEURON_COUNT << std::endl;
	outstream << HIDDEN_LAYER_COUNT << std::endl;
	outstream << OUTPUT_NEURON_COUNT << std::endl;
	outstream << TEST_COUNT << std::endl;
	outstream << MAX_TEST_COUNT << std::endl;
	outstream << TRAIN_EPOCH << std::endl;
	outstream << TRAIN_PERIOD << std::endl;
	outstream << TRAIN_EPS << std::endl;
	outstream << END_TRAIN_TICK << std::endl;

	//RMS config
	outstream << RMS_GAMMA << std::endl;
	outstream << RMS_LEARNRATE << std::endl;
	outstream << RMS_EPSILON << std::endl;

	for (int i = 1; i < countLayers; ++i)
		outstream << layers[i].Weights << std::endl;

	for (int i = 1; i < countLayers; ++i)
		outstream << layers[i].Bias << std::endl;
	outstream.flush();
}

Matrix2d MLP::GetOut() const
{
	return layers.back().Axons;
}

std::ostream & NeuroNet::operator<<(std::ostream & os, MLP & net)
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

void NeuroNet::MLP::SetInput(const Matrix2d & inputs)
{
	layers[0].States.copy(0, 0, inputNeuron, inputs);
}
#pragma endregion