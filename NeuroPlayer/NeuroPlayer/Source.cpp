#include <vector>
#include <iostream>
#include <algorithm>
#include "World.h"
#include "MyPlayer.h"
#include <fstream>
#include <ctime>
#include <exception>
#include <climits>
#include <omp.h>

using namespace std;
ofstream debug("neurodebug.txt");

const double ALPHA = 2.0;
const double GAMMA = 0.99;
const double LAMBDA = 0.5;
double r;
double Q;
double lastQ = 0.0;
int lastAction = 0;

//World config
const int PARAMS_COUNT = 3; //HEALTH, FULLNESS, ANGLE
const int FOOD_COUNT = 100;
const int POISON_COUNT = 0;
const int TRAP_COUNT = 0;
const int CORNUCOPIA_COUNT = 0;
const int BLOCK_COUNT = 0;
const int PLAYER_COUNT = 1;

//Neural config
const int SENSOR_COUNT = 8;
const int INPUT_NEURON_COUNT = SENSOR_COUNT * 2;
const int HIDDEN_NEURON_COUNT = INPUT_NEURON_COUNT * 1.5;
const int OUTPUT_NEURON_COUNT = 1;



namespace NeuroNet
{
	//activation function type
	enum AFType { LINE, SIGM, TANH };

	double getRand(double vmin, double vmax, bool integer = false);

	class Matrix2d
	{
	private:
		std::vector<std::vector<double>> _m;
	public:
		Matrix2d() {};
		Matrix2d(int n, int m, double val = 0.0);
		Matrix2d(const std::vector<double> &rhs);
		void Init(int n, int m, double val = 0.0);
		void InitRandom(int n, int m, double minv = -1.0, double maxv = 1.0);

		void Clear();
		int GetHorizontalSize() const;
		int GetVerticalSize() const;

		Matrix2d operator! () const;
		Matrix2d operator- () const;

		Matrix2d operator+= (const Matrix2d &rhs);
		Matrix2d operator+ (const Matrix2d &rhs) const;
		Matrix2d operator+= (const double rhs);
		Matrix2d operator+ (const double rhs) const;

		Matrix2d operator-= (const Matrix2d &rhs);
		Matrix2d operator- (const Matrix2d &rhs) const;
		Matrix2d operator-= (const double rhs);
		Matrix2d operator- (const double rhs) const;

		Matrix2d operator* (const Matrix2d &rhs) const;
		Matrix2d operator* (const double &rhs);

		Matrix2d operator= (const std::vector<std::vector<double>> &rhs);
		Matrix2d operator= (const Matrix2d &rhs);
		Matrix2d operator= (const std::vector<double> &rhs);

		Matrix2d abs() const;
		Matrix2d multiplication(const Matrix2d &rhs) const;
		const double sum() const;

		double& operator() (const int i, const int j);

		friend std::ostream& operator<< (std::ostream &os, const Matrix2d &m);
		friend Matrix2d sqrt(const Matrix2d&rhs);
	};

	struct Problem
	{
		Matrix2d inputs, outputs;
		Problem(std::vector<double> input, std::vector<double> output)
		{
			inputs = input;
			outputs = output;
		}
		Problem()
		{};
	};

	class Layer
	{
	private:
		AFType _aftype;

		Matrix2d sigm_function(Matrix2d m);
		Matrix2d tanh_function(Matrix2d m);
		Matrix2d diff_tanh_function(Matrix2d m);
		Matrix2d diff_sigm_function(Matrix2d m);
	public:
		Matrix2d Axons;
		Matrix2d States;
		Matrix2d Delta;
		Matrix2d LastDelta;
		Matrix2d Grad;
		Matrix2d LastGrad;
		Matrix2d CorrectVal;
		Matrix2d BiasCorrectVal;
		Matrix2d GetDiff();
		Matrix2d Weights;
		Matrix2d Bias;


		//----------------------------------
		//Functions
		Layer(int neuronCount, int prevNeuronCount, AFType activationFunction, bool bias = false);
		void CalculateStates(Layer &prevLayer);
		void CalculateAxons();
		void NguenWidrow(double Xmin, double Xmax, double Ymin, double Ymax);
	};

	class ElmanNetwork
	{
	protected:
		std::vector<Layer> _layers;
		std::vector<Matrix2d> LastDeltaSum;
		std::vector<Matrix2d> LastGradSum;
		int _countlayers;
	public:
		ElmanNetwork() {};
		ElmanNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction);
		std::vector<Problem> TrainingSet;
		std::vector<Matrix2d> eligibility;
		void Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction);
		//void Run();
		void Run();
		void Run(std::vector<double> input);
		void Run(Matrix2d &input);
		void AddTest(vector<double> ideal);
		double RunTrainingSetOffline(bool print = false);
		double CalculateError(Problem& test, bool print = false);
		void PrintProblemResult(Problem& test);
		void ResilientPropagation();
		void ResilientPropagation(std::vector<Matrix2d> PrevSumGrad, std::vector<Matrix2d> SumGrad, std::vector<Matrix2d> PrevSumDelta, std::vector<Matrix2d> SumDelta);
		void CalcGradDelta(double output);
		void CalcGradDelta(std::vector<double> output);
		void CalcGradDelta(Matrix2d &output);
		Matrix2d& GetOutputGrad();
		void MCQLCorrect();

		Matrix2d GetOut() const;
		friend std::ostream& operator<< (std::ostream &os, ElmanNetwork &net);
	};

	double getRand(double vmin, double vmax, bool integer)
	{
		if (integer)
			return (int)(vmin + rand() * (vmax - vmin) / RAND_MAX);
		return vmin + rand() * (vmax - vmin) / RAND_MAX;//
	}

	Matrix2d::Matrix2d(int n, int m, double val)
	{
		_m.resize(0);
		_m.resize(n, std::vector<double>(m, val));
	}

	Matrix2d::Matrix2d(const std::vector<double>& rhs)
	{
		Init(1, rhs.size());
		for (int i = 0; i < rhs.size(); ++i)
			_m[0][i] = rhs[i];
	}

	void Matrix2d::Init(int n, int m, double val)
	{
		_m.resize(0);
		_m.resize(n, std::vector<double>(m, val));
	}

	void Matrix2d::InitRandom(int n, int m, double minv, double maxv)
	{
		Init(n, m);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				_m[i][j] = getRand(minv, maxv);
			}
		}
	}

	void Matrix2d::Clear()
	{
		Init(GetVerticalSize(), GetHorizontalSize());
	}

	int Matrix2d::GetHorizontalSize() const
	{
		return (GetVerticalSize() > 0 ? (int)_m[0].size() : 0);
	}

	int Matrix2d::GetVerticalSize() const
	{
		return (int)_m.size();
	}

	Matrix2d Matrix2d::operator!() const
	{
		Matrix2d res;
		int m = GetVerticalSize();
		int n = GetHorizontalSize();
		res.Init(n, m);
		for (int i = 0; i < m; ++i)
			for (int j = 0; j < n; ++j)
				res._m[j][i] = _m[i][j];
		return res;
	}

	Matrix2d Matrix2d::operator-() const
	{
		Matrix2d res = *this;
		int n = GetVerticalSize();
		int m = GetHorizontalSize();
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				res._m[i][j] = -_m[i][j];
		return res;
	}

	Matrix2d Matrix2d::operator+=(const Matrix2d & rhs)
	{
		if (GetHorizontalSize() != rhs.GetHorizontalSize() || GetVerticalSize() != rhs.GetVerticalSize())
			throw std::logic_error("Wrong sizes in addition operation");

		for (int i = 0; i < (int)this->_m.size(); ++i)
			for (int j = 0; j < (int)this->_m[i].size(); ++j)
				_m[i][j] += rhs._m[i][j];
		return *this;
	}

	Matrix2d Matrix2d::operator+(const Matrix2d & rhs) const
	{
		Matrix2d res = *this;
		return res += rhs;
	}

	Matrix2d Matrix2d::operator+=(const double rhs)
	{
		for (int i = 0; i < (int)this->_m.size(); ++i)
			for (int j = 0; j < (int)this->_m[i].size(); ++j)
				_m[i][j] += rhs;
		return *this;
	}

	Matrix2d Matrix2d::operator+(const double rhs) const
	{
		Matrix2d res = *this;
		return res += rhs;
	}

	Matrix2d Matrix2d::operator-=(const Matrix2d & rhs)
	{
		if (GetHorizontalSize() != rhs.GetHorizontalSize() || GetVerticalSize() != rhs.GetVerticalSize())
			throw std::logic_error("Wrong sizes in subtraction operation");

		for (int i = 0; i < (int)this->_m.size(); ++i)
			for (int j = 0; j < (int)this->_m[i].size(); ++j)
				_m[i][j] -= rhs._m[i][j];
		return *this;
	}

	Matrix2d Matrix2d::operator-(const Matrix2d & rhs) const
	{
		Matrix2d res = *this;
		return res -= rhs;
	}

	Matrix2d Matrix2d::operator-=(const double rhs)
	{
		return *this += -rhs;
	}

	Matrix2d Matrix2d::operator-(const double rhs) const
	{
		Matrix2d res = *this;
		return res += -rhs;
	}

	Matrix2d Matrix2d::operator*(const Matrix2d & rhs) const
	{
		if (GetHorizontalSize() != rhs.GetVerticalSize())
			throw std::logic_error("Wrong sizes in multiplication operation");

		int n = GetVerticalSize();
		int m = rhs.GetHorizontalSize();
		int nm = GetHorizontalSize();
		Matrix2d res(n, m);

#pragma omp parallel for private(j,k,res._m[i][j])
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				for (int k = 0; k < nm; ++k)
					res._m[i][j] += _m[i][k] * rhs._m[k][j];
		return res;
	}

	Matrix2d Matrix2d::operator*(const double & rhs)
	{
		int n = GetVerticalSize();
		int m = GetHorizontalSize();
		Matrix2d res(n, m);

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				res._m[i][j] = _m[i][j] * rhs;
		return res;
	}

	Matrix2d Matrix2d::operator=(const std::vector<std::vector<double>>& rhs)
	{
		_m.resize(rhs.size());
		for (int i = 0; i < (int)rhs.size(); ++i)
		{
			_m[i].resize(rhs[i].size());
			for (int j = 0; j < (int)rhs[i].size(); ++j)
				_m[i][j] = rhs[i][j];
		}
		return *this;
	}

	Matrix2d Matrix2d::operator=(const Matrix2d & rhs)
	{
		_m.resize(rhs._m.size());
		for (int i = 0; i < (int)this->_m.size(); ++i)
		{
			_m[i].resize(rhs._m[i].size());
			for (int j = 0; j < (int)this->_m[i].size(); ++j)
			{
				this->_m[i][j] = rhs._m[i][j];
			}
		}
		return *this;
	}

	Matrix2d Matrix2d::operator=(const std::vector<double>& rhs)
	{
		this->_m.resize(1);
		_m[0].resize(rhs.size());
		std::copy(rhs.begin(), rhs.end(), _m[0].begin());
		return *this;
	}

	Matrix2d Matrix2d::abs() const
	{
		Matrix2d old;
		int n = GetHorizontalSize();
		int m = GetVerticalSize();
		old.Init(n, m);
		for (int i = 0; i < m; ++i)
			for (int j = 0; j < n; ++j)
				old._m[j][i] = std::abs(_m[i][j]);
		return old;
	}

	Matrix2d Matrix2d::multiplication(const Matrix2d & rhs) const
	{
		if (GetHorizontalSize() != rhs.GetHorizontalSize() || GetVerticalSize() != rhs.GetVerticalSize())
			throw std::logic_error("Wrong sizes in multiplication operation");

		Matrix2d res((int)rhs._m.size(), (int)rhs._m[0].size());
		for (int i = 0; i < (int)rhs._m.size(); ++i)
			for (int j = 0; j < (int)rhs._m[i].size(); ++j)
				res._m[i][j] = this->_m[i][j] * rhs._m[i][j];
		return res;
	}

	const double Matrix2d::sum() const
	{
		double ans = 0;
		for (int i = 0; i < (int)_m.size(); ++i)
			for (int j = 0; j < (int)_m[i].size(); ++j)
				ans += _m[i][j];
		return ans;
	}

	double & Matrix2d::operator()(const int i, const int j)
	{
		return _m[i][j];
	}

	std::ostream & operator<<(std::ostream & os, const Matrix2d & m)
	{
		for (int j = 0; j < m.GetVerticalSize(); ++j)
		{
			for (int k = 0; k < m.GetHorizontalSize(); ++k)
				os << m._m[j][k] << " ";
			os << std::endl;
		}

		return os;
	}

	Matrix2d sqrt(const Matrix2d & rhs)
	{
		Matrix2d res;
		int n = rhs.GetHorizontalSize();
		int m = rhs.GetVerticalSize();
		res.Init(n, m);
		for (int i = 0; i < m; ++i)
			for (int j = 0; j < n; ++j)
				res._m[j][i] = std::abs(rhs._m[i][j]);
		return res;
	}

	Matrix2d Layer::sigm_function(Matrix2d x)
	{
		int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
		Matrix2d res(n, m);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				res(i, j) = x(i, j) <= -35 ? x(i, j) = 10e-15 : 1.0 / (1.0 + exp(-x(i, j)));
		return res;
	}

	Matrix2d Layer::tanh_function(Matrix2d x)
	{
		int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
		Matrix2d res(n, m);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
			{
				//res(i,j) = (exp(2 * x(i,j)) - 1.0) / (exp(2 * x(i,j)) + 1.0);
				res(i, j) = std::tanh(x(i, j));
			}
		return res;
	}

	Matrix2d Layer::diff_tanh_function(Matrix2d x)
	{
		int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
		Matrix2d res(n, m);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				res(i, j) = 1.0 - x(i, j) * x(i, j);
		return res;
	}

	Matrix2d Layer::diff_sigm_function(Matrix2d x)
	{
		int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
		Matrix2d res(n, m);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				res(i, j) = (1.0 - x(i, j)) * x(i, j);
		return res;
	}

	Layer::Layer(int neuronCount, int prevNeuronCount, AFType activationFunction, bool bias)
	{
		_aftype = activationFunction;
		Weights.InitRandom(neuronCount, prevNeuronCount, -1, 1);
		States.Init(1, neuronCount);
		Axons.Init(1, neuronCount);
		Delta.Init(1, neuronCount);
		LastDelta.Init(1, neuronCount);
		Grad.Init(neuronCount, prevNeuronCount);
		LastGrad.Init(neuronCount, prevNeuronCount);
		CorrectVal.InitRandom(neuronCount, prevNeuronCount, 0, 1);
		if (bias)
			Bias.InitRandom(1, neuronCount);
		else Bias.Init(1, neuronCount);
		//BiasCorrectVal.Init(1, neuronCount, 0.1);
		BiasCorrectVal.InitRandom(1, neuronCount);
	}

	//TODO add const operation
	void Layer::CalculateStates(Layer & prevLayer)
	{
		States = prevLayer.Axons * !Weights + Bias;
	}

	void Layer::CalculateAxons()
	{
		switch (_aftype)
		{
		case SIGM:
			Axons = sigm_function(States);
			break;
		case LINE:
			Axons = States;
			break;
		case TANH:
			Axons = tanh_function(States);
			break;
		default:
			break;
		}
	}

	Matrix2d Layer::GetDiff()
	{
		switch (_aftype)
		{
		case SIGM:
			return diff_sigm_function(Axons);
		case LINE:
			return Matrix2d(Axons.GetVerticalSize(), Axons.GetHorizontalSize(), 1.0);
		case TANH:
			return diff_tanh_function(Axons);
		default:
			return Matrix2d(Axons.GetVerticalSize(), Axons.GetHorizontalSize(), 1.0);
			break;
		}
	}

	void Layer::NguenWidrow(double Xmin, double Xmax, double Ymin, double Ymax)
	{
		double beta = 0.7*std::pow(Weights.GetVerticalSize(), 1.0 / Weights.GetHorizontalSize());
		for (int i = 0; i < Weights.GetVerticalSize(); ++i)
		{
			double mvij = 0.0;
			for (int j = 0; j < Weights.GetHorizontalSize(); ++j)
				mvij += Weights(i, j) * Weights(i, j);
			mvij = std::sqrt(mvij);
			if (mvij == 0) mvij = 1;

			for (int j = 0; j < Weights.GetHorizontalSize(); ++j)
			{
				Weights(i, j) *= beta / mvij;
			}
			Bias(0, i) = getRand(-beta, beta);
		}

		double x = 0.5*(Xmax - Xmin), y = 0.5*(Xmax + Xmin);
		Weights = Weights*x; Bias = Bias*x + y;
		Matrix2d a(1, Weights.GetHorizontalSize()), c(1, Weights.GetHorizontalSize());
		for (int j = 0; j < Weights.GetHorizontalSize(); j++) {
			a(0, j) = 2.0 / (Ymax - Ymin);
			c(0, j) = 1.0 - Ymax*a(0, j);
		}
		Bias = !(Weights * !c + !Bias);

		for (int j = 0; j < Weights.GetVerticalSize(); j++)
			for (int k = 0; k < Weights.GetHorizontalSize(); k++)
				Weights(j, k) *= a(0, k);
	}

	double operator* (const std::vector<double> &lhs, const std::vector<double> &rhs)
	{
		double res = 0.0;
		for (int i = 0; i < (int)lhs.size(); ++i)
		{
			res += lhs[i] * rhs[i];
		}
		return res;
	}

	double ElmanNetwork::RunTrainingSetOffline(bool print)
	{
		std::vector<Matrix2d> GradSum(3);
		std::vector<Matrix2d> DeltaSum(3);

		for (int i = 0; i < 3; ++i)
		{
			GradSum[i] = Matrix2d(_layers[i].Grad.GetVerticalSize(), _layers[i].Grad.GetHorizontalSize());
			DeltaSum[i] = Matrix2d(_layers[i].Delta.GetVerticalSize(), _layers[i].Delta.GetHorizontalSize());
		}
		int TestCount = 30;

		for (int i = min((int)TrainingSet.size() - 1, TestCount); i >= 0; --i)
		{
			Run(TrainingSet[TrainingSet.size() - i - 1].inputs);
			CalcGradDelta(TrainingSet[TrainingSet.size() - i - 1].outputs);
			for (int i = 1; i < 3; ++i)
			{
				GradSum[i] += _layers[i].Grad;
				DeltaSum[i] += _layers[i].Delta;
			}
			if (print) PrintProblemResult(TrainingSet[TrainingSet.size() - i - 1]);
		}
		if (print) std::cout << std::endl << "=======CORRECT==========" << std::endl;

		double normGrad = 0.0;
		for (int i = 0; i < _countlayers; ++i)
			normGrad += sqrt(GradSum[i] * !GradSum[i]).sum();

		if (normGrad < 1e-6)
			return 0.0;
		ResilientPropagation(LastGradSum, GradSum, LastDeltaSum, DeltaSum);
		LastGradSum = GradSum;
		LastDeltaSum = DeltaSum;


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
			std::cout << test.outputs(0, i) << " ";
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
					double cur_mult = _layers[i].Grad(j, k) * _layers[i].LastGrad(j, k);
					if (cur_mult == 0.0)
						cur_correct = getRand(0, 1);
					else if (cur_mult > 0.0)
						cur_correct = min(EttaPlus * _layers[i].CorrectVal(j, k), 50.0);
					else if (cur_mult < 0.0)
						cur_correct = max(EttaMinus * _layers[i].CorrectVal(j, k), 1e-6);

					_layers[i].CorrectVal(j, k) = cur_correct;

					if (_layers[i].Grad(j, k) == 0.0) continue;

					if (_layers[i].Grad(j, k) > 0)
						_layers[i].Weights(j, k) += -cur_correct;
					else
						_layers[i].Weights(j, k) += cur_correct;
				}
			}

			for (int j = 0; j < _layers[i].Delta.GetHorizontalSize(); ++j)
			{
				double cur_correct = 0.0;
				double cur_mult = _layers[i].Delta(0, j) * _layers[i].LastDelta(0, j);
				if (cur_mult == 0.0)
					cur_correct = getRand(0, 1);
				else if (cur_mult > 0.0)
					cur_correct = min(EttaPlus * _layers[i].BiasCorrectVal(0, j), 50.0);
				else if (cur_mult < 0.0)
					cur_correct = max(EttaMinus * _layers[i].BiasCorrectVal(0, j), 1e-6);

				_layers[i].BiasCorrectVal(0, j) = cur_correct;

				if (_layers[i].Delta(0, j) == 0.0) continue;

				if (_layers[i].Delta(0, j) > 0)
					_layers[i].Bias(0, j) += -cur_correct;
				else
					_layers[i].Bias(0, j) += cur_correct;
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

	void ElmanNetwork::ResilientPropagation(std::vector<Matrix2d> PrevSumGrad, std::vector<Matrix2d> SumGrad, std::vector<Matrix2d> PrevSumDelta, std::vector<Matrix2d> SumDelta)
	{
		const double EttaPlus = 1.2, EttaMinus = 0.5;
		for (int i = 1; i < _countlayers; ++i)
		{
			for (int j = 0; j < SumGrad[i].GetVerticalSize(); ++j)
			{
				for (int k = 0; k < SumGrad[i].GetHorizontalSize(); ++k)
				{
					double cur_correct = 0.0;
					double cur_mult = SumGrad[i](j, k) * PrevSumGrad[i](j, k);
					if (cur_mult == 0.0)
						cur_correct = getRand(0, 1);
					else if (cur_mult > 0.0)
						cur_correct = min(EttaPlus * _layers[i].CorrectVal(j, k), 50.0);
					else if (cur_mult < 0.0)
						cur_correct = max(EttaMinus * _layers[i].CorrectVal(j, k), 1e-6);

					_layers[i].CorrectVal(j, k) = cur_correct;

					if (SumGrad[i](j, k) == 0.0) continue;

					if (SumGrad[i](j, k) > 0)
						_layers[i].Weights(j, k) += -cur_correct;
					else
						_layers[i].Weights(j, k) += cur_correct;
				}
			}

			for (int j = 0; j < SumDelta[i].GetHorizontalSize(); ++j)
			{
				double cur_correct = 0.0;
				double cur_mult = SumDelta[i](0, j) * PrevSumDelta[i](0, j);
				if (cur_mult == 0.0)
					cur_correct = getRand(0, 1);
				else if (cur_mult > 0.0)
					cur_correct = min(EttaPlus * _layers[i].BiasCorrectVal(0, j), 50.0);
				else if (cur_mult < 0.0)
					cur_correct = max(EttaMinus * _layers[i].BiasCorrectVal(0, j), 1e-6);

				_layers[i].BiasCorrectVal(0, j) = cur_correct;

				if (SumDelta[i](0, j) == 0.0) continue;

				if (SumDelta[i](0, j) > 0)
					_layers[i].Bias(0, j) += -cur_correct;
				else
					_layers[i].Bias(0, j) += cur_correct;
			}
		}
	}

	void ElmanNetwork::CalcGradDelta(double output)
	{
		CalcGradDelta(std::vector<double>(1, output));
	}

	void ElmanNetwork::CalcGradDelta(std::vector<double> outputs)
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
		try {
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
		}
	}

	Matrix2d ElmanNetwork::GetOut() const
	{
		return _layers.back().Axons;
	}

	std::ostream & operator<<(std::ostream & os, ElmanNetwork & net)
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
			os << net._layers[0].Axons(0, i) << " ";

		os << std::endl;

		os << std::endl << "Output neurons:" << std::endl;
		for (int i = 0; i < net._layers.back().Axons.GetHorizontalSize(); ++i)
			os << net._layers.back().Axons(0, i) << " ";
		os << std::endl;
		os << "===============" << std::endl;
		return os;
	}

	ElmanNetwork::ElmanNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
	{
		Init(InputCount, OutputCount, NeuronCount, HiddenLayerFunction);
	}

	void ElmanNetwork::Init(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
	{
		_layers.clear();
		_layers.push_back(Layer(InputCount + NeuronCount, 0, LINE));
		_layers.push_back(Layer(NeuronCount, InputCount + NeuronCount, HiddenLayerFunction, true));
		_layers.back().NguenWidrow(-2, 2, -1, 1);
		_layers.push_back(Layer(OutputCount, NeuronCount, LINE, true));
		//_layers.back().NguenWidrow(-1, 1, -1, 1);
		_countlayers = _layers.size();

		LastDeltaSum.resize(3);
		LastGradSum.resize(3);
		eligibility.resize(3);
		for (int i = 0; i < 3; ++i)
		{
			LastGradSum[i].Init(_layers[i].Grad.GetVerticalSize(), _layers[i].Grad.GetHorizontalSize());
			eligibility[i].Init(_layers[i].Grad.GetVerticalSize(), _layers[i].Grad.GetHorizontalSize());
			LastDeltaSum[i].Init(_layers[i].Delta.GetVerticalSize(), _layers[i].Delta.GetHorizontalSize());
		}

		debug << _layers[1].Weights(0, 0) << endl << endl;;
	}

	void ElmanNetwork::Run()
	{
		try {
			//init hidden and output layers
			for (int i = 1; i < (int)_layers.size(); ++i)
			{
				_layers[i].CalculateStates(_layers[i - 1]);
				_layers[i].CalculateAxons();
			}

			//copy hidden into input (context)
			for (int i = 1; i <= _layers[1].Axons.GetHorizontalSize(); ++i)
				_layers[0].States(0, _layers[0].States.GetHorizontalSize() - i) = _layers[1].Axons(0, _layers[1].Axons.GetHorizontalSize() - i);
		}
		catch (std::exception ex)
		{
			debug << "EXCEPTION: " << endl << ex.what() << endl;
			throw new std::exception("Exit");
		}
	}

	void ElmanNetwork::Run(std::vector<double> input)
	{
		Run(Matrix2d(input));
	}

	void ElmanNetwork::Run(Matrix2d& input)
	{
		//init input layer
		for (int i = 0; i < (int)input.GetHorizontalSize(); ++i)
			_layers[0].States(0, i) = input(0, i);
		_layers[0].CalculateAxons();

		Run();
	}
	void ElmanNetwork::AddTest(vector<double> ideal)
	{
		TrainingSet.resize(TrainingSet.size() + 1);
		TrainingSet.back().inputs = _layers[0].States;
		TrainingSet.back().outputs = ideal;
	}
}

enum Actions { FORWARD, BACKWARD, LEFTSTEP, RIGHTSTEP, COUNT };

void DoAction(MyPlayer* me, Actions action)
{
	switch (action)
	{
	case Actions::FORWARD:
		me->StepForward();
		break;

	case Actions::BACKWARD:
		me->StepBackward();
		break;

	case Actions::LEFTSTEP:
		me->StepLeft();
		break;

	case Actions::RIGHTSTEP:
		me->StepRight();
		break;

	default:
		break;
	}
}

bool FirstStep = true;
vector<NeuroNet::ElmanNetwork> nets;
//NeuroNet::ElmanNetwork net;
void MyPlayer::Init()
{
	SetName(L"NeuroPlayer");

	nets.resize(Actions::COUNT);
	for (int i = 0; i < nets.size(); ++i)
	{
		debug << "net " << i << endl;
		nets[i].Init(INPUT_NEURON_COUNT, OUTPUT_NEURON_COUNT, HIDDEN_NEURON_COUNT, NeuroNet::AFType::SIGM);
		debug << "net " << i << "inited" << endl << endl;
	}
	//net.Init(INPUT_NEURON_COUNT, OUTPUT_NEURON_COUNT, HIDDEN_NEURON_COUNT, NeuroNet::AFType::TANH);
}

bool check(Element *player, Element *elem, double angle, double step_angle)
{
	double b = player->GetY() - elem->GetY();
	double k = tan(angle);
	double r = elem->GetR();
	double x0 = elem->GetX();
	double angleB = atan((player->GetY() - elem->GetY())*1.0 / (elem->GetX() - player->GetX()));
	return (k*b - x0)*(k*b - x0) - (1.0 + k*k)*(x0*x0 + b*b - r*r) >= 0.0 && abs(angleB - angle) <= step_angle;
}

enum ElementType { TFOOD, TENEMY, TBLOCK, TCOUNT };

double getDistanceOnWall(double x, double y, Player *me, double angle)
{
	double x0 = me->GetX();
	double y0 = me->GetY();
	double vx = x - x0;
	double vy = y - y0;
	double beamx = cos(angle) - x0;
	double beamy = sin(angle) - y0;
	//scalar product
	if (vx*beamx + vy*beamy > 0.0)
		return me->GetDistanceTo(x, y);
	else return DBL_MAX;

}

void setInput(NeuroNet::Matrix2d &inputs, Player *me, const int eyes, World *w)
{
	const double step_angle = M_PI*2.0 / eyes;
	double angle = 0.0;

	vector<pair<double, ElementType>> sensors(eyes, make_pair(DBL_MAX, ElementType::TCOUNT));
	double dist;
	int eye = 0;
	while (angle < 2 * M_PI)
	{
		double min_dist = DBL_MAX;
		ElementType cur_type = ElementType::TBLOCK;
		//food
		for each (auto cur in w->GetFood())
		{
			dist = me->GetDistanceTo(cur);
			if (check(me, cur, angle, step_angle) && dist < min_dist)
			{
				min_dist = dist;
				cur_type = ElementType::TFOOD;
			}
		}

		//block
		for each (auto cur in w->GetBlocks())
		{
			dist = me->GetDistanceTo(cur);
			if (check(me, cur, angle, step_angle) && dist < min_dist)
			{
				min_dist = dist;
				cur_type = ElementType::TBLOCK;
			}
		}

		double x, y;
		double x0 = me->GetX();
		double y0 = me->GetY();
		double k = tan(angle);
		double angleB;

		dist = DBL_MAX;
		//x=0
		dist = min(getDistanceOnWall(0, y0, me, angle), dist);
		//x = getw
		dist = min(getDistanceOnWall(w->GetWidth(), k*w->GetWidth() + y0, me, angle), dist);
		//y = 0;
		dist = min(getDistanceOnWall(-y0 / k, 0, me, angle), dist);
		//y = geth;
		dist = min(getDistanceOnWall((w->GetHeight() - y0) / k, w->GetHeight(), me, angle), dist);
		if (dist < min_dist)
		{
			min_dist = dist;
			cur_type = ElementType::TBLOCK;
		}

		inputs(0, eye) = min_dist;
		inputs(0, eye + 1) = cur_type;

		eye += 2;
		angle += step_angle;
	}


}

void MyPlayer::Move()
{
	//Input order: My coordinates, Health, Fullness, Angle, Food coordinates, Enemy coordinates, Trap coordinates, Poison coordinates, Cornucopia coordinates, Block coordinates
	NeuroNet::Matrix2d inputs(1, INPUT_NEURON_COUNT, -1.0);

	/////////////////////////////////////////////////////////////////////
	///////////////////////// Init input vector /////////////////////////
	setInput(inputs, this, SENSOR_COUNT, GetWorld());
	/////////////////////////////////////////////////////////////////////

	r = GetHealth() * GetFullness() * 1.0 / 300000;

	int action = -1;
	int rnd = rand();
	debug << "RAND: " << rnd;
	if (!FirstStep)
	{
		//nets[lastAction].MCQLCorrect();
		//nets[lastAction].CalcGradDelta(r);
		//nets[lastAction].AddTest(vector<double>(1, r));
		int Epoch = 10;
		for (int i = 0; i < Actions::COUNT; ++i)
		{
			nets[i].AddTest(vector<double>(1, r));
			/*while (Epoch--)
			{
			if (nets[i].RunTrainingSetOffline() < 1e-2)
			break;
			}*/
		}
		while (Epoch--)
		{
			if (nets[lastAction].RunTrainingSetOffline() < 1e-2)
				break;
		}
	}

	debug << "R = " << r << endl;
	if (rnd % 2 != 0)
	{

		Q = -DBL_MAX;
		for (int i = 0; i < nets.size(); ++i)
		{
			nets[i].Run(inputs);
			double curQ = nets[i].GetOut().sum();
			debug << curQ << " " << i << endl;

			if (Q < curQ)
			{
				Q = curQ;
				action = i;
			}
		}
	}
	else
	{
		action = rnd % Actions::COUNT;
	}
	debug << "SELECT " << Q << " " << action << endl;


	DoAction(this, (Actions)action);
	FirstStep = false;

	lastAction = action;
	lastQ = Q;
}

#ifdef _DEBUG
int main()
{
	MyPlayer *a = new MyPlayer();
	a->Init();
	a->Move();
	return 0;
}
#endif