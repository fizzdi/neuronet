#include "ElmanNetwork.h"
#include <algorithm>
#include <string>
#include <queue>

using namespace NeuroNet;

ElmanNetwork::ElmanNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
{
	Init(InputCount, OutputCount, NeuronCount, HiddenLayerFunction);
}

NeuroNet::ElmanNetwork::ElmanNetwork(const std::string & FileName)
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

Matrix2d NeuroNet::ElmanNetwork::GetContext()
{
	Matrix2d res(1, context_neuron);
	res.copy(inputNeuron, 0, context_neuron, layers[0].Axons);
	return std::move(res);
}

void NeuroNet::ElmanNetwork::SetContext(const Matrix2d & context)
{
	layers[0].States.copy(0, inputNeuron, context_neuron, context);
}

double NeuroNet::ElmanNetwork::RMSTraining(training_set& TrainingSet)
{
	for (int i = 1; i < countLayers; ++i)
	{
		layers[i].GradSum.Fill(0.0);
		layers[i].DeltaSum.Fill(0.0);
	}

	std::queue<int> tests;
	Matrix2d Context = GetContext();
	if (TrainingSet.size() < TEST_COUNT)
	{
		for (int test = 0; test < TrainingSet.size(); ++test)
		{
			tests.push(test);
			SetContext(Context);
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
			SetContext(Context);
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
		normGrad += layers[i].GradSum.multiplication(layers[i].GradSum).sum();
	normGrad = std::sqrt(normGrad);
	if (normGrad < 1e-2)
		return 0.0;

	RMSPropagation();

	double error = 0.0;
	while (!tests.empty())
	{
		SetContext(Context);
		Run(TrainingSet[tests.front()].inputs);
		error += CalculateError(TrainingSet[tests.front()]);
		tests.pop();
	}
	SetContext(Context);
	return error;
}

void NeuroNet::ElmanNetwork::Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
{
	layers.clear();
	context_neuron = NeuronCount;
	inputNeuron = InputCount;
	outputNeuron = OutputCount;
	layers.emplace_back(Layer(InputCount + context_neuron, 0, LINE));
	layers.emplace_back(Layer(NeuronCount, InputCount + context_neuron, HiddenLayerFunction));
	layers.back().NguenWidrow(-2, 2, -1, 1);
	for (int i = 1; i < HIDDEN_LAYER_COUNT; ++i)
	{
		layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
		layers.back().NguenWidrow(-2, 2, -1, 1);
	}
	layers.emplace_back(Layer(OutputCount, NeuronCount, LINE));
	countLayers = (int)layers.size();
}

void ElmanNetwork::Run()
{
	//init hidden and output layers

	for (int i = 1; i < countLayers; ++i)
	{
		layers[i].CalculateStates(layers[i - 1]);
		layers[i].CalculateAxons();
	}

	SetContext(layers[countLayers - 2].Axons);
}