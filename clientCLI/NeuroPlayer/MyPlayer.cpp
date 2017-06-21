#define _CRT_SECURE_NO_WARNINGS
#define _CRTDBG_MAP_ALLOC #include <stdlib.h> #include <crtdbg.h>
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
#include "Matrix2d.h"
#include "ElmanNetwork.h"

using namespace std;

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

deque<NeuroNet::Problem> TrainingSet;

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
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
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
		nets[lastAction].AddTest(TrainingSet, vr);
		if (tick % TRAIN_PERIOD == 0)
		{
			int Epoch = TRAIN_EPOCH;
			while (Epoch--)
			{
				if (nets[lastAction].RunTrainingSetOffline(TrainingSet) < TRAIN_EPS)
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
			nets[i].debuginfo(debug);
		}
		debug.flush();
	}
}