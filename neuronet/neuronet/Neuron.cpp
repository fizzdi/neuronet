#include "Neuron.h"
#include <cmath>
#include <fstream>
#include <iostream>

NeuroNet::Neuron::Neuron()
{
	_state = 0;
	_axon = 0;
	_aftype = SIGM;
}

NeuroNet::Neuron::Neuron(AFType aftype)
{
	_state = 0;
	_axon = 0;
	_aftype = aftype;
}

void NeuroNet::Neuron::SetState(double state)
{
	this->_state = state;
}

double NeuroNet::Neuron::GetAxon()
{
	return _axon;
}

void NeuroNet::Neuron::CalculateAxon()
{
	switch (_aftype)
	{
	case SIGM:
		if (_axon <= -35) _axon = 10e-15;
		else _axon = 1.0 / (1.0 + exp(-_state));
		break;
	case LINE:
		_axon = _state;
		break;
	case TANH:
		_axon = (exp(2*_state)-1.0)/(exp(2*_state)+1.0);
		break;
	default:
		if (_axon <= -35) _axon = 10e-15;
		else _axon = 1.0 / (1.0 + exp(-_state));
		break;
	}

	std::ofstream f("debug.txt", std::ios::app);
	f << (_aftype == SIGM ? "SIGM" : "LINE") << " f(" << _state << ") = " << _axon << std::endl;
	f.close();
}

double NeuroNet::Neuron::Delta()
{
	return _delta;
}

double NeuroNet::Neuron::Delta(double delta)
{
	return _delta = delta;
}

void NeuroNet::Neuron::SetLastdw(double ldw)
{
}

double NeuroNet::Neuron::GetLastdw()
{
	return 0.0;
}

double NeuroNet::Neuron::GetDiff()
{
	switch (_aftype)
	{
	case SIGM:
		return (1.0 - _axon)*_axon;
	case TANH:
		return 1.0 - _axon*_axon;
	default:
		return 1;
	}
}
