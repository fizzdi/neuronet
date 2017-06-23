#include "ElmanNetwork.h"
#include <algorithm>

using namespace NeuroNet;

ElmanNetwork::ElmanNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
{
	Layers.clear();
	Layers.emplace_back(Layer(InputCount + NeuronCount, 0, LINE));
	Layers.emplace_back(Layer(NeuronCount, InputCount + NeuronCount, HiddenLayerFunction));
	Layers.back().NguenWidrow(-2, 2, -1, 1);
	Layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
	Layers.back().NguenWidrow(-2, 2, -1, 1);
	/*Layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
	Layers.back().NguenWidrow(-2, 2, -1, 1);*/
	Layers.emplace_back(Layer(OutputCount, NeuronCount, LINE));
	_countlayers = Layers.size();

	eligibility.resize(_countlayers);
	for (int i = 0; i < _countlayers; ++i)
	{
		eligibility[i] = Matrix2d(Layers[i].Grad.GetVerticalSize(), Layers[i].Grad.GetHorizontalSize());
		eligibility[i].Fill(0.0);
	}
}

void ElmanNetwork::Run()
{
	for (int i = 0; i < Layers.size(); ++i)
	{
		//debug << "Layers[" << i << "].Axons: " << std::endl << Layers[i].Axons.GetHorizontalSize() << " " << Layers[i].Axons.GetVerticalSize() << std::endl;
	}

	//init hidden and output layers
	for (int i = 1; i < (int)Layers.size(); ++i)
	{
		Layers[i].CalculateStates(Layers[i - 1]);
		Layers[i].CalculateAxons();
		//debug << "Layers["<<i<<"].Axons: " << std::endl << Layers[i].Axons << std::endl;
	}

	//copy hidden into input (context)
	for (int i = 1; i <= Layers[1].Axons.GetHorizontalSize(); ++i)
		Layers[0].States.at(0, Layers[0].States.GetHorizontalSize() - i) = Layers[1].Axons.at(0, Layers[1].Axons.GetHorizontalSize() - i);
	//debug << "Layers[0].States: " << std::endl << Layers[0].States << std::endl;
}