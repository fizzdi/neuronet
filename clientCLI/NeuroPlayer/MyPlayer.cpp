#include <vector>
#include <iostream>
#include <algorithm>
#include "World.h"
#include "MyPlayer.h"
#include <fstream>
#include <ctime>
#include <climits>
#include <omp.h>
#include <memory>
#include <deque>
#include "Matrix2d.h"
#include "ElmanNetwork.h"
#include <queue>
using namespace std;

ofstream debug("neurodebug.txt");

enum Actions { FORWARD, BACKWARD, LEFTSTEP, RIGHTSTEP, COUNT };

//q learning
const double QL_LEARN = 0.9;
const int COUNT_REWARD = 50;

int lastAction = 0;

deque<NeuroNet::Problem> TrainingSet;

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
NeuroNet::ElmanNetwork net;

bool check(Element *player, double x, double y, double r, double angle, double step_angle)
{
	double xfe = cos(angle);
	double yfe = sin(angle);
	double xse = cos(angle + step_angle);
	double yse = sin(angle + step_angle);
	double xo = x - player->GetX();
	double yo = y - player->GetY();

	//angle first eye - object
	double acos_val = (xfe*xo + yfe*yo) / sqrt(xo*xo + yo*yo);
	double angleA = acos(acos_val);
	if (angleA < 0.0)
		angleA += M_PI;
	//angle second eye - object
	acos_val = (xse*xo + yse*yo) / sqrt(xo*xo + yo*yo);
	double angleB = acos(acos_val);
	if (angleB < 0.0)
		angleB += M_PI;

	return angleA - step_angle <= DBL_EPSILON && angleB - step_angle <= DBL_EPSILON;
}

bool check(Element *player, Element *elem, double angle, double step_angle)
{
	return check(player, elem->GetX(), elem->GetY(), elem->GetR(), angle, step_angle);
}

enum ElementType { TFOOD, TENEMY, TBLOCK, TCOUNT };

pair<double, ElementType> getDistanceOnWall(Player *me, World *w, double angle, double step_angle)
{
	double x0 = me->GetX();
	double y0 = me->GetY();
	double k = -tan(angle);

	//Y = 0, x = 0..getw
	double dist = DBL_MAX;
	double min_dist = DBL_MAX;
	ElementType cur_type = ElementType::TENEMY;
	if (abs(angle - M_PI / 2) <= DBL_EPSILON)
	{
		//пересечение только с горизонтальными сторонами
		dist = min(dist, me->GetY());
		if (min_dist > dist)
		{
			min_dist = dist;
			cur_type = ElementType::TBLOCK;
		}
	}
	else if (abs(angle - 3 * M_PI / 2) <= DBL_EPSILON)
	{
		//пересечение только с горизонтальными сторонами
		dist = min(dist, w->GetHeight() - me->GetY());
		if (min_dist > dist)
		{
			min_dist = dist;
			cur_type = ElementType::TBLOCK;
		}
	}
	else if (abs(angle) <= DBL_EPSILON)
	{
		//пересечение только с вертилкальными сторонами
		dist = min(dist, w->GetWidth() - me->GetX());
		if (min_dist > dist)
		{
			min_dist = dist;
			cur_type = ElementType::TBLOCK;
		}
	}
	else if (abs(angle - M_PI) <= DBL_EPSILON)
	{
		//пересечение только с вертилкальными сторонами
		dist = min(dist, me->GetX());
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

		double origin_vx = cos(angle), origin_vy = -sin(angle);

		double vx = x - x0, vy = y - y0;
		double cosv = (origin_vx*vx + origin_vy * vy) / sqrt(vx*vx + vy*vy);
		if (x >= -EPS && x <= w->GetWidth() + EPS && /*check(me, x, y, angle, step_angle)*/ abs(cosv - 1.0) <= EPS)
		{
			dist = min(dist, sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0)));
		}
		//y=geth
		y = w->GetHeight();
		x = (y - b) / k;
		vx = x - x0, vy = y - y0;
		cosv = (origin_vx*vx + origin_vy * vy) / sqrt(vx*vx + vy*vy);
		if (x >= -EPS && x <= w->GetWidth() + EPS && abs(cosv - 1.0) <= EPS)
		{
			dist = min(dist, sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0)));
		}

		//x=0
		x = 0.0;
		y = k*x + b;
		vx = x - x0, vy = y - y0;
		cosv = (origin_vx*vx + origin_vy * vy) / sqrt(vx*vx + vy*vy);
		if (y >= -EPS && y <= w->GetHeight() + EPS && abs(cosv - 1.0) <= EPS)
		{
			dist = min(dist, sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0)));
		}

		x = w->GetWidth();
		y = k*x + b;
		vx = x - x0, vy = y - y0;
		cosv = (origin_vx*vx + origin_vy * vy) / sqrt(vx*vx + vy*vy);
		if (y >= -EPS && y <= w->GetHeight() + EPS && abs(cosv - 1.0) <= EPS)
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
	double angle = M_PI / 100;

	vector<pair<double, ElementType>> sensors(eyes, make_pair(DBL_MAX, ElementType::TCOUNT));
	int eye = 0;
	while (angle < 2 * M_PI)
	{
		double min_dist = DBL_MAX;
		double dist = DBL_MAX;
		ElementType cur_type = ElementType::TENEMY;
		//food
		int food_counter = 0;
		for each (auto cur in w->GetFood())
		{
			dist = me->GetDistanceTo(cur);
			if (check(me, cur, angle, step_angle))
			{
				if (dist < min_dist)
				{
					min_dist = dist;
					cur_type = ElementType::TFOOD;
				}
				food_counter++;
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

		inputs.at(0, 2 * eye) = min_dist / 800;
		inputs.at(0, 2 * eye + 1) = cur_type;

		//inputs.at(0, 3 * eye) = min_dist;
		/*inputs.at(0, 3 * eye) = min_dist / (640 * 480);
		inputs.at(0, 3 * eye + 1) = cur_type;
		/*int val = 0;
		switch (cur_type)
		{
		case ElementType::TFOOD:
			val = food_counter;
		define:
			val = 1;
		}
		inputs.at(0, 3 * eye + 2) = val;*/

		eye++;
		angle += step_angle;
	}


}
NeuroNet::Matrix2d inputs(1, INPUT_NEURON_COUNT);
NeuroNet::Matrix2d last_inputs;

int tick = -1;
NeuroNet::Matrix2d lastQ;
double last_full = 0.0;
void MyPlayer::Init()
{
	SetName(L"NeuroPlayer");
	SetEyeCount(SENSOR_COUNT);
	net = NeuroNet::ElmanNetwork(INPUT_NEURON_COUNT, OUTPUT_NEURON_COUNT, HIDDEN_NEURON_COUNT, NeuroNet::AFType::TANH);
}
void MyPlayer::Move()
{
	tick++;

	/////////////////////////////////////////////////////////////////////
	///////////////////////// InitFill input vector /////////////////////////
	last_inputs = inputs;
	setInput(inputs, this, SENSOR_COUNT, GetWorld());
	/////////////////////////////////////////////////////////////////////
	int rnd = NeuroNet::myrand();
	int action = rnd%Actions::COUNT;
	double dist_on_wall = min(GetX(), min(GetY(), min(GetWorld()->GetHeight() - GetY(), GetWorld()->GetWidth() - GetX())));
	double wallPenalty = +0.5 * exp(-dist_on_wall / 50.0);
	double reward = this->GetFullness() - last_full + wallPenalty;
	//debug << "rew " << reward << endl;
	last_full = this->GetFullness();
	if (!FirstStep)
	{
		NeuroNet::Matrix2d Context = net.GetContext();
		net.Run(inputs);
		NeuroNet::Matrix2d Q = net.GetOut();
		double tmp = -DBL_MAX;
		for (int i = 0; i < Q.GetHorizontalSize(); ++i)
		{
			if (tmp < Q.at(0, i))
			{
				tmp = Q.at(0, i);
				action = i;
			}
		}
		Q.at(0, lastAction) = reward + QL_LEARN*tmp;
		net.SetContext(Context);
		NeuroNet::AddTest(TrainingSet, last_inputs, Q);
		if (tick % TRAIN_PERIOD == 0)
		{
			int Epoch = TRAIN_EPOCH;
				double error;
			while (Epoch--)
			{
				if ((error = net.RMSTraining(TrainingSet)) < TRAIN_EPS)
					break;
			//debug << error << std::endl;
			}
			//debug << std::endl;
		}
		net.SetContext(Context);
	}

	NeuroNet::Matrix2d Q;
	if (tick < 5000 && (rnd % 28) < 5)
	{
		action = rnd % Actions::COUNT;
		lastAction = action;
	}
	else
	{
		net.Run(inputs);
		Q = net.GetOut();
		double tmp = -DBL_MAX;
		for (int i = 0; i < Q.GetHorizontalSize(); ++i)
		{
			if (tmp < Q.at(0, i))
			{
				tmp = Q.at(0, i);
				action = i;
			}
		}
	}

	DoAction(this, (Actions)action);
	FirstStep = false;
	lastAction = action;
	lastQ = Q;
	if (tick > 0 && tick % 3000 == 0)
	{
		debug << std::endl << std::endl << "========================================================" << tick << "========================================================" << std::endl;
		net.PrintFullInfo(debug);
		debug.flush();
	}
}