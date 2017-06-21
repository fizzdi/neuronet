#define _CRT_SECURE_NO_WARNINGS
#include <vector>
#include <iostream>
#include <algorithm>
#include "World.h"
#include "MyPlayer.h"
#include <fstream>
#include <ctime>
#include <climits>
#include <omp.h>
#include <random>
#include <memory>
#include <deque>

using namespace std;
#define EPS (1e-4)
#define ERRORDEF "###ERROR"
ofstream debug("neurodebug.txt");
std::mt19937 eng;
std::uniform_int_distribution<> dist(1, 50000);

const double ALPHA = 2.0;
const double GAMMA = 0.99;
const double LAMBDA = 0.5;
//double r;
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
const int HIDDEN_NEURON_COUNT = 40;
const int OUTPUT_NEURON_COUNT = 1;
const int TEST_COUNT = 300;
const int TRAIN_EPOCH = 10;
const int TRAIN_PERIOD = 5;
const double TRAIN_EPS = 1e-1;

namespace NeuroNet
{
	class Matrix2d
	{
	private:
		double* _m;
		int n, m;
	public:
		Matrix2d() { n = m = 0; _m = nullptr; };
		~Matrix2d();
		Matrix2d(const Matrix2d& rhs);
		Matrix2d(Matrix2d&& rhs);
		Matrix2d(int n, int m);
		Matrix2d(const std::vector<double> &rhs);
		void Fill(double val = 0.0);
		void InitRandom(double minv = -1.0, double maxv = 1.0);

		int GetHorizontalSize() const;
		int GetVerticalSize() const;

		Matrix2d operator! () const;
		Matrix2d operator- () const;

		Matrix2d& operator+= (const Matrix2d &rhs);
		Matrix2d operator+ (const Matrix2d &rhs) const;
		Matrix2d& operator+= (const double rhs);
		Matrix2d operator+ (const double rhs) const;

		Matrix2d& operator-= (const Matrix2d &rhs);
		Matrix2d operator- (const Matrix2d &rhs) const;
		Matrix2d& operator-= (const double rhs);
		Matrix2d operator- (const double rhs) const;

		Matrix2d operator* (const Matrix2d &rhs) const;
		Matrix2d operator* (const double &rhs);

		Matrix2d& operator= (const std::vector<std::vector<double>> &rhs);
		Matrix2d& operator= (const Matrix2d &rhs);
		Matrix2d& operator= (Matrix2d &&rhs);
		Matrix2d& operator= (const std::vector<double> &rhs);

		Matrix2d abs() const;
		Matrix2d multiplication(const Matrix2d &rhs) const;
		const double sum() const;

		double& at(const int i, const int j);
		double at(const int i, const int j) const;

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
}

deque<NeuroNet::Problem> TrainingSet;

namespace NeuroNet
{
	//activation function type
	enum AFType { LINE, SIGM, TANH };

	double getRand(double vmin, double vmax, bool integer = false);

	class Layer
	{
	private:
		AFType _aftype;

		Matrix2d sigm_function(Matrix2d& m);
		Matrix2d tanh_function(Matrix2d& m);
		Matrix2d diff_tanh_function(Matrix2d& m);
		Matrix2d diff_sigm_function(Matrix2d& m);
	public:
		Matrix2d Axons;
		Matrix2d States;
		Matrix2d Delta;
		Matrix2d LastDelta;
		Matrix2d Grad;
		Matrix2d GradSum;
		Matrix2d DeltaSum;
		Matrix2d LastGrad;
		Matrix2d CorrectVal;
		Matrix2d BiasCorrectVal;
		Matrix2d Weights;
		Matrix2d Bias;
		Matrix2d LastDeltaSum;
		Matrix2d LastGradSum;

		//----------------------------------
		//Functions
		Layer(int neuronCount, int prevNeuronCount, AFType activationFunction);
		void CalculateStates(Layer &prevLayer);
		void CalculateAxons();
		void NguenWidrow(double Xmin, double Xmax, double Ymin, double Ymax);
		Matrix2d GetDiff();
	};

	class ElmanNetwork
	{
	protected:
		int _countlayers;
	public:
		std::vector<Layer> _layers;
		ElmanNetwork() {};
		ElmanNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction);
		std::vector<Matrix2d> eligibility;
		void Run();
		void Run(const std::vector<double>& input);
		void Run(Matrix2d &input);
		void AddTest(const vector<double>& ideal) const;
		void AddTest(Matrix2d& ideal) const;
		double RunTrainingSetOffline(bool print = false);
		double CalculateError(Problem& test, bool print = false);
		void PrintProblemResult(Problem& test);
		void ResilientPropagation();
		void ResilientPropagationOffline();
		void CalcGradDelta(const double output);
		void CalcGradDelta(const std::vector<double>& output);
		void CalcGradDelta(Matrix2d &output);
		Matrix2d& GetOutputGrad();
		void MCQLCorrect();
		void debuginfo()
		{
			for (int i = 1; i < _countlayers; ++i)
				debug << _layers[i].Weights << endl;
			debug << endl;
		}
		Matrix2d GetOut() const;
		friend std::ostream& operator<< (std::ostream &os, ElmanNetwork &net);
	};

	double getRand(double vmin, double vmax, bool integer)
	{
		if (integer)
			return (int)(vmin + rand() * (vmax - vmin) / RAND_MAX);
		return vmin + rand() * (vmax - vmin) / RAND_MAX;//
	}

	Matrix2d::~Matrix2d()
	{
		if (_m == nullptr) return;
		delete[](_m);
		_m = nullptr;
	}

	Matrix2d::Matrix2d(const Matrix2d & rhs)
	{
		if (_m == rhs._m) return;
		this->n = rhs.n;
		this->m = rhs.m;
		_m = new double[n*m];
		memcpy_s(_m, n*m * sizeof(*_m), rhs._m, n*m * sizeof(*rhs._m));
	}

	Matrix2d::Matrix2d(Matrix2d && rhs)
	{
		_m = move(rhs._m);
		n = rhs.n;
		m = rhs.m;
		rhs._m = nullptr;
	}

	Matrix2d::Matrix2d(int n, int m)
	{
		this->n = n;
		this->m = m;
		_m = new double[n*m];
	}

	Matrix2d::Matrix2d(const std::vector<double>& rhs)
	{
		this->n = 1;
		this->m = rhs.size();
		_m = new double[m];
		memcpy_s(_m, n*m * sizeof(*_m), rhs.data(), n*m * sizeof(*rhs.data()));
	}

	void Matrix2d::Fill(double val)
	{
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				at(i, j) = val;
	}

	void Matrix2d::InitRandom(double minv, double maxv)
	{
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				at(i, j) = getRand(minv, maxv);
	}

	int Matrix2d::GetHorizontalSize() const
	{
		return m;
	}

	int Matrix2d::GetVerticalSize() const
	{
		return n;
	}

	Matrix2d Matrix2d::operator!() const
	{
		int m = GetVerticalSize();
		int n = GetHorizontalSize();
		Matrix2d res(n, m);
		for (int i = 0; i < m; ++i)
			for (int j = 0; j < n; ++j)
				res.at(j, i) = this->at(i, j);
		return res;
	}

	Matrix2d Matrix2d::operator-() const
	{
		Matrix2d res = *this;
		int n = GetVerticalSize();
		int m = GetHorizontalSize();
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				res.at(i, j) = -at(i, j);
		return res;
	}

	Matrix2d& Matrix2d::operator+=(const Matrix2d & rhs)
	{
		if (GetHorizontalSize() != rhs.GetHorizontalSize() || GetVerticalSize() != rhs.GetVerticalSize())
		{
			debug << ERRORDEF << " " << string(__FILE__) << "(" << __LINE__ << "):" << string(__FUNCTION__) << endl;
			throw std::logic_error("Wrong sizes in addition operation");
		}

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				at(i, j) += rhs.at(i, j);
		return *this;
	}

	Matrix2d Matrix2d::operator+(const Matrix2d & rhs) const
	{
		Matrix2d res = *this;
		return res += rhs;
	}

	Matrix2d& Matrix2d::operator+=(const double rhs)
	{
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				at(i, j) += rhs;
		return *this;
	}

	Matrix2d Matrix2d::operator+(const double rhs) const
	{
		Matrix2d res = *this;
		return res += rhs;
	}

	Matrix2d& Matrix2d::operator-=(const Matrix2d & rhs)
	{
		if (GetHorizontalSize() != rhs.GetHorizontalSize() || GetVerticalSize() != rhs.GetVerticalSize())
		{
			debug << ERRORDEF << " " << string(__FILE__) << "(" << __LINE__ << "):" << string(__FUNCTION__) << endl;
			throw std::logic_error("Wrong sizes in subtraction operation");
		}

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				at(i, j) -= rhs.at(i, j);
		return *this;
	}

	Matrix2d Matrix2d::operator-(const Matrix2d & rhs) const
	{
		Matrix2d res = *this;
		return res -= rhs;
	}

	Matrix2d& Matrix2d::operator-=(const double rhs)
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
		{
			debug << ERRORDEF << " " << string(__FILE__) << "(" << __LINE__ << "):" << string(__FUNCTION__) << endl;
			throw std::logic_error("Wrong sizes in multiplication operation");
		}

		int n = GetVerticalSize();
		int m = rhs.GetHorizontalSize();
		int nm = GetHorizontalSize();
		Matrix2d res(n, m);
		res.Fill(0.0);

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				for (int k = 0; k < nm; ++k)
					res.at(i, j) += this->at(i, k) * rhs.at(k, j);
		return res;
	}

	Matrix2d Matrix2d::operator*(const double & rhs)
	{
		int n = GetVerticalSize();
		int m = GetHorizontalSize();
		Matrix2d res(n, m);

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				res.at(i, j) = at(i, j) * rhs;
		return res;
	}

	Matrix2d& Matrix2d::operator=(const std::vector<std::vector<double>>& rhs)
	{

		if (_m != nullptr)
			delete[](_m);
		this->n = rhs.size();
		this->m = n > 0 ? rhs[0].size() : 0;
		_m = new double[n*m];

		for (int i = 0; i < (int)rhs.size(); ++i)
			memcpy_s(_m + (i*m) * sizeof(*_m), m * sizeof(*_m), rhs.data(), m * sizeof(*rhs.data()));
		return *this;
	}

	Matrix2d& Matrix2d::operator=(const Matrix2d & rhs)
	{
		if (rhs._m == _m)
		{
			return *this;
		}

		if (_m != nullptr)
			delete[](_m);


		this->n = rhs.n;
		this->m = rhs.m;
		_m = new double[n*m];

		memcpy_s(_m, n*m * sizeof(*_m), rhs._m, n*m * sizeof(*rhs._m));
		return *this;
	}

	Matrix2d & Matrix2d::operator=(Matrix2d && rhs)
	{
		if (this->_m == rhs._m) return *this;
		_m = move(rhs._m);
		n = rhs.n;
		m = rhs.m;

		rhs._m = nullptr;
		return *this;
	}

	Matrix2d& Matrix2d::operator=(const std::vector<double>& rhs)
	{
		if (_m != nullptr)
			delete[](_m);
		this->n = 1;
		this->m = rhs.size();
		_m = new double[n*m];
		memcpy_s(_m, n*m * sizeof(*_m), rhs.data(), n*m * sizeof(*rhs.data()));
		return *this;
	}

	Matrix2d Matrix2d::abs() const
	{
		int n = GetVerticalSize();
		int m = GetHorizontalSize();
		Matrix2d old(n, m);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				old.at(i, i) = std::abs(at(i, j));
		return old;
	}

	Matrix2d Matrix2d::multiplication(const Matrix2d & rhs) const
	{
		if (GetHorizontalSize() != rhs.GetHorizontalSize() || GetVerticalSize() != rhs.GetVerticalSize())
		{
			debug << ERRORDEF << " " << string(__FILE__) << "(" << __LINE__ << "):" << string(__FUNCTION__) << endl;
			throw std::logic_error("Wrong sizes in multiplication operation");
		}

		Matrix2d res(rhs.n, rhs.m);
		for (int i = 0; i < rhs.n; ++i)
			for (int j = 0; j < rhs.m; ++j)
				res.at(i, j) = this->at(i, j) * rhs.at(i, j);
		return res;
	}

	const double Matrix2d::sum() const
	{
		double ans = 0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				ans += at(i, j);
		return ans;
	}

	double & Matrix2d::at(const int i, const int j)
	{
		if (i >= n || j >= m)
		{
			debug << ERRORDEF << " " << string(__FILE__) << "(" << __LINE__ << "):" << string(__FUNCTION__) << endl;
			throw std::logic_error("Wrong sizes in at operation");
		}
		return _m[i*m + j];
	}

	double Matrix2d::at(const int i, const int j) const
	{
		if (i >= n || j >= m)
		{
			debug << ERRORDEF << " " << string(__FILE__) << "(" << __LINE__ << "):" << string(__FUNCTION__) << endl;
			throw std::logic_error("Wrong sizes in at operation");
		}
		return _m[i*m + j];
	}

	std::ostream & operator<<(std::ostream & os, const Matrix2d & m)
	{
		for (int j = 0; j < m.GetVerticalSize(); ++j)
		{
			for (int k = 0; k < m.GetHorizontalSize(); ++k)
				os << m.at(j, k) << " ";
			os << std::endl;
		}

		return os;
	}

	Matrix2d sqrt(const Matrix2d & rhs)
	{
		int n = rhs.GetVerticalSize();
		int m = rhs.GetHorizontalSize();
		Matrix2d res(n, m);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				res.at(i, j) = std::sqrt(rhs.at(i, j));
		return res;
	}

	Matrix2d Layer::sigm_function(Matrix2d& x)
	{
		int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
		Matrix2d res(n, m);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				res.at(i, j) = x.at(i, j) <= -35 ? x.at(i, j) = 10e-15 : 1.0 / (1.0 + exp(-x.at(i, j)));
		return res;
	}

	Matrix2d Layer::tanh_function(Matrix2d& x)
	{
		int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
		Matrix2d res(n, m);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
			{
				//res(i,j) = (exp(2 * x(i,j)) - 1.0) / (exp(2 * x(i,j)) + 1.0);
				res.at(i, j) = std::tanh(x.at(i, j));
			}
		return res;
	}

	Matrix2d Layer::diff_tanh_function(Matrix2d& x)
	{
		int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
		Matrix2d res(n, m);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				res.at(i, j) = 1.0 - x.at(i, j) * x.at(i, j);
		return res;
	}

	Matrix2d Layer::diff_sigm_function(Matrix2d& x)
	{
		int n = x.GetVerticalSize(), m = x.GetHorizontalSize();
		Matrix2d res(n, m);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				res.at(i, j) = (1.0 - x.at(i, j)) * x.at(i, j);
		return res;
	}

	Layer::Layer(int neuronCount, int prevNeuronCount, AFType activationFunction)
	{
		_aftype = activationFunction;
		Weights = Matrix2d(neuronCount, prevNeuronCount);
		States = Matrix2d(1, neuronCount);
		Axons = Matrix2d(1, neuronCount);
		Delta = Matrix2d(1, neuronCount);
		LastDelta = Matrix2d(1, neuronCount);
		Grad = Matrix2d(neuronCount, prevNeuronCount);
		LastGrad = Matrix2d(neuronCount, prevNeuronCount);
		CorrectVal = Matrix2d(neuronCount, prevNeuronCount);
		DeltaSum = Matrix2d(1, neuronCount);
		GradSum = Matrix2d(neuronCount, prevNeuronCount);
		LastDeltaSum = Matrix2d(1, neuronCount);
		LastGradSum = Matrix2d(neuronCount, prevNeuronCount);
		Bias = Matrix2d(1, neuronCount);
		BiasCorrectVal = Matrix2d(1, neuronCount);

		Weights.InitRandom(0, 1);
		CorrectVal.InitRandom(0, 1);
		DeltaSum.Fill(0.0);
		GradSum.Fill(0.0);
		LastDeltaSum.Fill(0.0);
		LastGradSum.Fill(0.0);
		BiasCorrectVal.InitRandom(-1.0, 1.0);
		Bias.InitRandom(-1.0, 1.0);
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
		{
			Matrix2d res(Axons.GetVerticalSize(), Axons.GetHorizontalSize());
			res.Fill(1.0);
			return res;
		}
		case TANH:
			return diff_tanh_function(Axons);
		default:
			Matrix2d res(Axons.GetVerticalSize(), Axons.GetHorizontalSize());
			res.Fill(1.0);
			return res;
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
				mvij += Weights.at(i, j) * Weights.at(i, j);
			mvij = std::sqrt(mvij);
			if (mvij == 0) mvij = 1;

			for (int j = 0; j < Weights.GetHorizontalSize(); ++j)
			{
				Weights.at(i, j) *= beta / mvij;
			}
			Bias.at(0, i) = getRand(-beta, beta);
		}

		double x = 0.5*(Xmax - Xmin), y = 0.5*(Xmax + Xmin);
		Weights = Weights*x; Bias = Bias*x + y;
		Matrix2d a(1, Weights.GetHorizontalSize()), c(1, Weights.GetHorizontalSize());
		for (int j = 0; j < Weights.GetHorizontalSize(); j++) {
			a.at(0, j) = 2.0 / (Ymax - Ymin);
			c.at(0, j) = 1.0 - Ymax*a.at(0, j);
		}
		Bias = !(Weights * !c + !Bias);

		for (int j = 0; j < Weights.GetVerticalSize(); j++)
			for (int k = 0; k < Weights.GetHorizontalSize(); k++)
				Weights.at(j, k) *= a.at(0, k);
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
						cur_correct = min(EttaPlus * _layers[i].CorrectVal.at(j, k), 50.0);
					else if (cur_mult < 0.0)
						cur_correct = max(EttaMinus * _layers[i].CorrectVal.at(j, k), 1e-6);

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
					cur_correct = min(EttaPlus * _layers[i].BiasCorrectVal.at(0, j), 50.0);
				else if (cur_mult < 0.0)
					cur_correct = max(EttaMinus * _layers[i].BiasCorrectVal.at(0, j), 1e-6);

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
						cur_correct = min(EttaPlus * _layers[i].CorrectVal.at(j, k), 50.0);
					else if (cur_mult < 0.0)
						cur_correct = max(EttaMinus * _layers[i].CorrectVal.at(j, k), 1e-6);

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
					cur_correct = min(EttaPlus * _layers[i].BiasCorrectVal.at(0, j), 50.0);
				else if (cur_mult < 0.0)
					cur_correct = max(EttaMinus * _layers[i].BiasCorrectVal.at(0, j), 1e-6);

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
			os << net._layers[0].Axons.at(0, i) << " ";

		os << std::endl;

		os << std::endl << "Output neurons:" << std::endl;
		for (int i = 0; i < net._layers.back().Axons.GetHorizontalSize(); ++i)
			os << net._layers.back().Axons.at(0, i) << " ";
		os << std::endl;
		os << "===============" << std::endl;
		return os;
	}

	ElmanNetwork::ElmanNetwork(int InputCount, int OutputCount, int NeuronCount, AFType HiddenLayerFunction)
	{
		_layers.clear();
		_layers.emplace_back(Layer(InputCount + NeuronCount, 0, LINE));
		_layers.emplace_back(Layer(NeuronCount, InputCount + NeuronCount, HiddenLayerFunction));
		_layers.back().NguenWidrow(-2, 2, -1, 1);
		_layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
		_layers.back().NguenWidrow(-2, 2, -1, 1);
		_layers.emplace_back(Layer(NeuronCount, NeuronCount, HiddenLayerFunction));
		_layers.back().NguenWidrow(-2, 2, -1, 1);
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
	void ElmanNetwork::AddTest(const vector<double> &ideal) const
	{
		if (TrainingSet.size() + 1 == TEST_COUNT)
			TrainingSet.pop_front();
		auto pr = Problem();
		pr.inputs = _layers[0].States;
		pr.outputs = ideal;
		TrainingSet.emplace_back(pr);
	}
	void ElmanNetwork::AddTest(Matrix2d& ideal) const
	{
		if (TrainingSet.size() + 1 == TEST_COUNT)
			TrainingSet.pop_front();
		auto pr = Problem();
		pr.inputs = _layers[0].States;
		pr.outputs = ideal;
		TrainingSet.emplace_back(pr);
	}
}

NeuroNet::Matrix2d vr(1, 1);
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

SYSTEMTIME st;
bool FirstStep = true;
vector<NeuroNet::ElmanNetwork> nets;
//NeuroNet::ElmanNetwork net;
void MyPlayer::Init()
{
	SetName(L"NeuroPlayer");

	nets.resize(Actions::COUNT);
	for (int i = 0; i < nets.size(); ++i)
	{
		nets[i] = NeuroNet::ElmanNetwork(INPUT_NEURON_COUNT, OUTPUT_NEURON_COUNT, HIDDEN_NEURON_COUNT, NeuroNet::AFType::TANH);
	}
	//net.InitFill(INPUT_NEURON_COUNT, OUTPUT_NEURON_COUNT, HIDDEN_NEURON_COUNT, NeuroNet::AFType::TANH);
}

bool check(Element *player, double x, double y, double angle, double step_angle)
{
	double vx = 10 * cos(angle);
	double vy = 10 * sin(angle);
	double vx2 = x - player->GetX();
	double vy2 = y - player->GetY();
	double cosa = (vx*vx2 + vy*vy2) *1.0 /
		(
			sqrt(vx*vx + vy*vy)
			*sqrt(vx2*vx2 + vy2*vy2)
			);
	cosa = max(-1.0 + EPS, cosa);
	cosa = min(1.0 - EPS, cosa);
	double angleB = acos(cosa);
	if (angleB < -FLT_MAX)
		int y = 0;
	return abs(angleB) <= step_angle / 2.0 + EPS;
}


bool check(Element *player, Element *elem, double angle, double step_angle)
{
	return check(player, elem->GetX(), elem->GetY(), angle, step_angle);
}


enum ElementType { TFOOD, TENEMY, TBLOCK, TCOUNT };
string elty[] = { "FOOD", "ENEMY", "BLOCK", "COUNT" };

pair<double, ElementType> getDistanceOnWall(Player *me, World *w, double angle, double step_angle)
{
	double x, y;
	double x0 = me->GetX();
	double y0 = me->GetY();
	double k = tan(angle);
	double angleB;

	//Y = 0, x = 0..getw
	double dist = DBL_MAX;
	double min_dist = DBL_MAX;
	ElementType cur_type = ElementType::TENEMY;
	if (abs(angle - M_PI / 2) <= DBL_EPSILON || abs(angle - 3 * M_PI / 2) <= DBL_EPSILON)
	{
		//пересечение только с горизонтальными сторонами
		dist = min(dist, me->GetY());
		dist = min(dist, w->GetHeight() - me->GetY());
		if (min_dist > dist)
		{
			min_dist = dist;
			cur_type = ElementType::TBLOCK;
		}
	}
	else if (abs(angle) <= DBL_EPSILON || abs(angle - M_PI) <= DBL_EPSILON)
	{
		//пересечение только с вертилкальными сторонами
		dist = min(dist, me->GetX());
		dist = min(dist, w->GetWidth() - me->GetX());
		if (min_dist > dist)
		{
			min_dist = dist;
			cur_type = ElementType::TBLOCK;
		}
	}
	else
	{
		//случайные пересечения
		/*
		система:
		y = 0, при 0.0 <= x <= getw();
		y = geth(), при 0.0 <= x <= getw();

		x = 0, при 0.0 <= y <= geth();
		x = getw(), при 0.0 <= y <= geth();
		*/

		//y = 0
		double b = y0 - k*x0;
		double y = 0.0;
		double x = (y - b) / k;
		if (x >= -EPS && x <= w->GetWidth() + EPS && check(me, x, y, angle, step_angle))
		{
			dist = min(dist, sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0)));
		}
		//y=geth
		y = w->GetHeight();
		x = (y - b) / k;
		if (x >= -EPS && x <= w->GetWidth() + EPS && check(me, x, y, angle, step_angle))
		{
			dist = min(dist, sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0)));
		}

		//x=0
		x = 0.0;
		y = k*x + b;
		if (y >= -EPS && y <= w->GetHeight() + EPS && check(me, x, y, angle, step_angle))
		{
			dist = min(dist, sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0)));
		}

		//x=getw
		x = w->GetWidth();
		y = k*x + b;
		if (y >= -EPS && y <= w->GetHeight() + EPS && check(me, x, y, angle, step_angle))
		{
			dist = min(dist, sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0)));
		}

		if (min_dist > dist)
		{
			min_dist = dist;
			cur_type = ElementType::TBLOCK;
		}
	}
	return make_pair(min_dist, cur_type);
}

void setInput(NeuroNet::Matrix2d &inputs, Player *me, const int eyes, World *w)
{
	const double step_angle = M_PI*2.0 / eyes;
	double angle = 0.0;

	vector<pair<double, ElementType>> sensors(eyes, make_pair(DBL_MAX, ElementType::TCOUNT));
	int eye = 0;
	while (angle < 2 * M_PI)
	{
		double min_dist = DBL_MAX;
		double dist = DBL_MAX;
		ElementType cur_type = ElementType::TENEMY;
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

		if (cur_type == TENEMY)
			int y = 0;

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
		if (cur_type == TENEMY)
			int y = 0;


		auto res = getDistanceOnWall(me, w, angle, step_angle);
		if (res.second == TENEMY)
			int y = 0;

		if (min_dist > res.first)
		{
			min_dist = res.first;
			cur_type = res.second;
		}
		if (cur_type == TENEMY)
			int y = 0;
		inputs.at(0, eye) = min_dist;
		inputs.at(0, eye + 1) = cur_type;

		eye += 2;
		angle += step_angle;
	}


}
NeuroNet::Matrix2d inputs(1, INPUT_NEURON_COUNT);

int lasttest = -1;
int tick = -1;
void MyPlayer::Move()
{
	tick++;

	/////////////////////////////////////////////////////////////////////
	///////////////////////// InitFill input vector /////////////////////////
	inputs.Fill(-1.0);
	setInput(inputs, this, SENSOR_COUNT, GetWorld());
	/////////////////////////////////////////////////////////////////////
	vr.at(0, 0) = GetFullness();
	int action = -1;
	int rnd = dist(eng);

	if (!FirstStep)
	{
		int test = TrainingSet.size();
		nets[lastAction].AddTest(vr);
		if (tick % TRAIN_PERIOD == 0)
		{
			int Epoch = TRAIN_EPOCH;
			while (Epoch--)
			{
				if (nets[lastAction].RunTrainingSetOffline() < TRAIN_EPS)
					break;
			}
		}
	}

	if ((dist(eng) % 17) < 3)
	{
		Q = vr.at(0, 0);
		action = rnd % Actions::COUNT;
	}
	else
	{

		Q = -DBL_MAX;
		for (int i = 0; i < nets.size(); ++i)
		{
			nets[i].Run(inputs);
			double curQ = nets[i].GetOut().sum();

			if (Q < curQ)
			{
				Q = curQ;
				action = i;
			}
		}
	}

	DoAction(this, (Actions)action);
	FirstStep = false;
	lastAction = action;
	lastQ = Q;
	if (tick % 1000 == 0)
	{
		debug << endl << endl << "======================tick " << tick << "========================================================" << endl;
		for (int i = 0; i < Actions::COUNT; ++i)
		{
			nets[i].debuginfo();
		}
		debug.flush();
	}
}