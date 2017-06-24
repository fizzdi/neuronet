#include "ElmanNetwork.h"
#include <algorithm>
#include <string>

using namespace NeuroNet;

ElmanNetwork::ElmanNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
{
	layers.clear();
	context_neuron = NeuronCount;
	layers.emplace_back(Layer(InputCount + context_neuron, 0, LINE));
	layers.emplace_back(Layer(NeuronCount, InputCount + context_neuron, HiddenLayerFunction));
	layers.back().NguenWidrow(-2, 2, -1, 1);
	layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
	layers.back().NguenWidrow(-2, 2, -1, 1);
	/*layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
	layers.back().NguenWidrow(-2, 2, -1, 1);*/
	layers.emplace_back(Layer(OutputCount, NeuronCount, LINE));
	countLayers = (int)layers.size();
}

Matrix2d NeuroNet::ElmanNetwork::GetContext()
{
	Matrix2d res(1, context_neuron);
	for (int i = 1; i <= context_neuron; ++i)
		res.at(0, context_neuron - i) = layers[0].States.at(0, layers[0].States.GetHorizontalSize() - i);
	return std::move(res);
}

void NeuroNet::ElmanNetwork::SetContext(const Matrix2d & context)
{
	for (int i = 1; i <= context_neuron; ++i)
		layers[0].States.at(0, layers[0].States.GetHorizontalSize() - i) = context.at(0, context_neuron - i);
}

double NeuroNet::ElmanNetwork::RMSTraining(training_set& TrainingSet)
{
	{
		for (int i = 1; i < countLayers; ++i)
		{
			layers[i].GradSum.Fill(0.0);
			layers[i].DeltaSum.Fill(0.0);
		}

		Matrix2d Context = GetContext();
		for (int t = 0; t < TEST_COUNT; ++t)
		{
			int test = myrand() % TrainingSet.size();
			SetContext(Context);
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
			SetContext(Context);
			Run(test.inputs);
			error += CalculateError(test);
		}
		SetContext(Context);
		return error;
	}
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