#pragma once

//activation function type
enum AFType {LINE, SIGM, TANH};

namespace NeuroNet
{
	class Neuron {
	private:
		double _state;
		double _axon;
		AFType _aftype;
		double _delta;
		double _last_dw;
	public:
		Neuron();
		Neuron(AFType aftype);
		void SetState(double state);
		double GetAxon();
		void CalculateAxon();
		double Delta();
		double Delta(double delta);
		void SetLastdw(double ldw);
		double	GetLastdw();
		double GetDiff();
	};
}