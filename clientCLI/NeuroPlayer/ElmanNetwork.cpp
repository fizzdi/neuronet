#include "ElmanNetwork.h"
#include <algorithm>

using namespace NeuroNet;

ElmanNetwork::ElmanNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
{
	_layers.clear();
	_layers.emplace_back(Layer(InputCount + NeuronCount, 0, LINE));
	_layers.emplace_back(Layer(NeuronCount, InputCount + NeuronCount, HiddenLayerFunction));
	_layers.back().NguenWidrow(-2, 2, -1, 1);
	_layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
	_layers.back().NguenWidrow(-2, 2, -1, 1);
	/*_layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
	_layers.back().NguenWidrow(-2, 2, -1, 1);*/
	_layers.emplace_back(Layer(OutputCount, NeuronCount, LINE));
	//_layers.back().NguenWidrow(-1, 1, -1, 1);
	_countlayers = _layers.size();

	eligibility.resize(_countlayers);
	for (int i = 0; i < _countlayers; ++i)
	{
		eligibility[i] = Matrix2d(_layers[i].Grad.GetVerticalSize(), _layers[i].Grad.GetHorizontalSize());
		eligibility[i].Fill(0.0);
	}
}

double ElmanNetwork::RunTrainingSetOffline(std::deque<Problem> &TrainingSet, bool print)
{
	for (int i = 0; i < _countlayers; ++i)
	{
		_layers[i].GradSum.Fill(0.0);
		_layers[i].DeltaSum.Fill(0.0);
	}

	for (int t = 0; t < TrainingSet.size(); ++t)
	{
		Run(TrainingSet[t].inputs);
		CalcGradDelta(TrainingSet[t].outputs);
		for (int i = 1; i < _countlayers; ++i)
		{
			_layers[i].GradSum += _layers[i].Grad;
			_layers[i].DeltaSum += _layers[i].Delta;
		}
		if (print) PrintProblemResult(TrainingSet[t]);
	}

	if (print) std::cout << std::endl << "=======CORRECT==========" << std::endl;

	double normGrad = 0.0;
	for (int i = 0; i < _countlayers; ++i)
		normGrad += sqrt(_layers[i].GradSum * !_layers[i].GradSum).sum();

	if (normGrad < 1e-6)
		return 0.0;
	ResilientPropagationOffline();
	for (int i = 0; i < _countlayers; ++i)
	{
		_layers[i].LastGradSum = _layers[i].GradSum;
		_layers[i].LastDeltaSum = _layers[i].DeltaSum;
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

void ElmanNetwork::PrintProblemResult(Problem & test)
{
	std::cout << *this;
	std::cout << "Expected results:" << std::endl;
	for (int i = 0; i < test.outputs.GetHorizontalSize(); ++i)
		std::cout << test.outputs.at(0, i) << " ";
	std::cout << std::endl;
}

void ElmanNetwork::ResilientPropagation()
{
	const double EttaPlus = 1.2, EttaMinus = 0.5;
	for (int i = 1; i < _countlayers; ++i)
	{
		for (int j = 0; j < _layers[i].Grad.GetVerticalSize(); ++j)
		{
			for (int k = 0; k < _layers[i].Grad.GetHorizontalSize(); ++k)
			{
				double cur_correct = 0.0;
				double cur_mult = _layers[i].Grad.at(j, k) * _layers[i].LastGrad.at(j, k);
				if (cur_mult == 0.0)
					cur_correct = getRand(0, 1);
				else if (cur_mult > 0.0)
					cur_correct = std::min(EttaPlus * _layers[i].CorrectVal.at(j, k), 50.0);
				else if (cur_mult < 0.0)
					cur_correct = std::max(EttaMinus * _layers[i].CorrectVal.at(j, k), 1e-6);

				_layers[i].CorrectVal.at(j, k) = cur_correct;

				if (_layers[i].Grad.at(j, k) == 0.0) continue;

				if (_layers[i].Grad.at(j, k) > 0)
					_layers[i].Weights.at(j, k) += -cur_correct;
				else
					_layers[i].Weights.at(j, k) += cur_correct;
			}
		}

		for (int j = 0; j < _layers[i].Delta.GetHorizontalSize(); ++j)
		{
			double cur_correct = 0.0;
			double cur_mult = _layers[i].Delta.at(0, j) * _layers[i].LastDelta.at(0, j);
			if (cur_mult == 0.0)
				cur_correct = getRand(0, 1);
			else if (cur_mult > 0.0)
				cur_correct = std::min(EttaPlus * _layers[i].BiasCorrectVal.at(0, j), 50.0);
			else if (cur_mult < 0.0)
				cur_correct = std::max(EttaMinus * _layers[i].BiasCorrectVal.at(0, j), 1e-6);

			_layers[i].BiasCorrectVal.at(0, j) = cur_correct;

			if (_layers[i].Delta.at(0, j) == 0.0) continue;

			if (_layers[i].Delta.at(0, j) > 0)
				_layers[i].Bias.at(0, j) += -cur_correct;
			else
				_layers[i].Bias.at(0, j) += cur_correct;
		}
	}
}

double ElmanNetwork::CalculateError(Problem & test, bool print)
{
	//MSE
	double error = (test.outputs - _layers.back().Axons).multiplication(test.outputs - _layers.back().Axons).sum() / 2.0;
	//RootMSE
	//double error = (test.outputs - _layers.back().Axons).abs().sqrt().sum() / 2.0;

	if (print)
	{
		std::cout << "MSE: " << error << std::endl;
		std::cout << "---------------------------------------------------------------------------------" << std::endl;
	}
	return error;
}

void ElmanNetwork::ResilientPropagationOffline()
{
	const double EttaPlus = 1.2, EttaMinus = 0.5;
	for (int i = 1; i < _countlayers; ++i)
	{
		for (int j = 0; j < _layers[i].GradSum.GetVerticalSize(); ++j)
		{
			for (int k = 0; k < _layers[i].GradSum.GetHorizontalSize(); ++k)
			{
				double cur_correct = 0.0;
				double cur_mult = _layers[i].GradSum.at(j, k) * _layers[i].LastGradSum.at(j, k);
				if (cur_mult == 0.0)
					cur_correct = getRand(0, 1);
				else if (cur_mult > 0.0)
					cur_correct = std::min(EttaPlus * _layers[i].CorrectVal.at(j, k), 50.0);
				else if (cur_mult < 0.0)
					cur_correct = std::max(EttaMinus * _layers[i].CorrectVal.at(j, k), 1e-6);

				_layers[i].CorrectVal.at(j, k) = cur_correct;

				if (_layers[i].GradSum.at(j, k) == 0.0) continue;

				if (_layers[i].GradSum.at(j, k) > 0)
					_layers[i].Weights.at(j, k) += -cur_correct;
				else
					_layers[i].Weights.at(j, k) += cur_correct;
			}
		}

		for (int j = 0; j < _layers[i].DeltaSum.GetHorizontalSize(); ++j)
		{
			double cur_correct = 0.0;
			double cur_mult = _layers[i].DeltaSum.at(0, j) * _layers[i].LastDeltaSum.at(0, j);
			if (cur_mult == 0.0)
				cur_correct = getRand(0, 1);
			else if (cur_mult > 0.0)
				cur_correct = std::min(EttaPlus * _layers[i].BiasCorrectVal.at(0, j), 50.0);
			else if (cur_mult < 0.0)
				cur_correct = std::max(EttaMinus * _layers[i].BiasCorrectVal.at(0, j), 1e-6);

			_layers[i].BiasCorrectVal.at(0, j) = cur_correct;

			if (_layers[i].DeltaSum.at(0, j) == 0.0) continue;

			if (_layers[i].DeltaSum.at(0, j) > 0)
				_layers[i].Bias.at(0, j) += -cur_correct;
			else
				_layers[i].Bias.at(0, j) += cur_correct;
		}
	}
}

void ElmanNetwork::CalcGradDelta(const double output)
{
	CalcGradDelta(std::vector<double>(1, output));
}

void ElmanNetwork::CalcGradDelta(const std::vector<double>& outputs)
{
	CalcGradDelta(Matrix2d(outputs));
}
void ElmanNetwork::CalcGradDelta(Matrix2d &outputs)
{
	for (int i = 0; i < _countlayers; ++i)
	{
		_layers[i].LastDelta = _layers[i].Delta;
		_layers[i].LastGrad = _layers[i].Grad;
	}

	_layers.back().Delta = (_layers.back().Axons - outputs).multiplication(_layers.back().GetDiff());

	for (int i = _countlayers - 2; i >= 0; --i)
		_layers[i].Delta = (_layers[i + 1].Delta * _layers[i + 1].Weights).multiplication(_layers[i].GetDiff());

	for (int i = 1; i < _countlayers; ++i)
		_layers[i].Grad = !_layers[i].Delta * _layers[i - 1].Axons;
}

Matrix2d& ElmanNetwork::GetOutputGrad()
{
	return _layers.back().Grad;
}

void ElmanNetwork::MCQLCorrect()
{
	/*try {
	for (int i = 1; i < _countlayers; ++i)
	_layers[i].Weights += eligibility[i] * ALPHA*(r + GAMMA*Q - lastQ);

	std::vector<double> out;
	out.push_back(Q);
	CalcGradDelta(out);

	for (int i = 1; i < _countlayers; ++i)
	eligibility[i] = _layers[i].Grad + eligibility[i] * GAMMA*LAMBDA;
	}
	catch (std::exception ex)
	{
	debug << "EXCEPTION: " << endl << ex.what() << endl;
	throw new std::exception("Exit");
	}*/
}

Matrix2d ElmanNetwork::GetOut() const
{
	return _layers.back().Axons;
}


void ElmanNetwork::Run()
{
	//init hidden and output layers
	for (int i = 1; i < (int)_layers.size(); ++i)
	{
		_layers[i].CalculateStates(_layers[i - 1]);
		_layers[i].CalculateAxons();
	}

	//copy hidden into input (context)
	for (int i = 1; i <= _layers[1].Axons.GetHorizontalSize(); ++i)
		_layers[0].States.at(0, _layers[0].States.GetHorizontalSize() - i) = _layers[1].Axons.at(0, _layers[1].Axons.GetHorizontalSize() - i);
}

void ElmanNetwork::Run(const std::vector<double>& input)
{
	Run(Matrix2d(input));
}

void ElmanNetwork::Run(Matrix2d& input)
{
	//init input layer
	for (int i = 0; i < (int)input.GetHorizontalSize(); ++i)
		_layers[0].States.at(0, i) = input.at(0, i);
	_layers[0].CalculateAxons();

	Run();
}
void ElmanNetwork::AddTest(std::deque<Problem> &TrainingSet, const std::vector<double> &ideal) const
{
	if (TrainingSet.size() + 1 == TEST_COUNT)
		TrainingSet.pop_front();
	auto pr = Problem();
	pr.inputs = _layers[0].States;
	pr.outputs = ideal;
	TrainingSet.emplace_back(pr);
}
void ElmanNetwork::AddTest(std::deque<Problem> &TrainingSet, Matrix2d& ideal) const
{
	if (TrainingSet.size() + 1 == TEST_COUNT)
		TrainingSet.pop_front();
	auto pr = Problem();
	pr.inputs = _layers[0].States;
	pr.outputs = ideal;
	TrainingSet.emplace_back(pr);
}

std::ostream & NeuroNet::operator<<(std::ostream & os, ElmanNetwork & net)
{
	for (int i = 1; i < (int)net._layers.size(); ++i)
	{
		os << "Weights " << i << " -> " << i - 1 << ": " << std::endl;
		os << net._layers[i].Weights << std::endl;
	}

	for (int i = 1; i < (int)net._layers.size(); ++i)
	{
		os << "Bias " << i << " -> " << i - 1 << ": " << std::endl;
		os << net._layers[i].Bias << std::endl;
	}

	os << std::endl << "Input neurons:" << std::endl;
	for (int i = 0; i < net._layers[0].Axons.GetHorizontalSize(); ++i)
		os << net._layers[0].Axons.at(0, i) << " ";

	os << std::endl;

	os << std::endl << "Output neurons:" << std::endl;
	for (int i = 0; i < net._layers.back().Axons.GetHorizontalSize(); ++i)
		os << net._layers.back().Axons.at(0, i) << " ";
	os << std::endl;
	os << "===============" << std::endl;
	return os;
}
