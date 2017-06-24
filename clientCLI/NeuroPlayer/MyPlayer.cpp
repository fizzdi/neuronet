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
using namespace std;

ofstream debug("neurodebug.txt");

enum Actions { FORWARD, BACKWARD, LEFTSTEP, RIGHTSTEP, COUNT };

//q learning
const double QL_LEARN = 0.9;

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
void MyPlayer::Init()
{
	SetName(L"NeuroPlayer");
	SetEyeCount(SENSOR_COUNT);
	net = NeuroNet::ElmanNetwork(INPUT_NEURON_COUNT, OUTPUT_NEURON_COUNT, HIDDEN_NEURON_COUNT, NeuroNet::AFType::SIGM);
}

bool check(Element *player, double x, double y, double r, double angle, double step_angle)
{
	double x0 = player->GetX(), y0 = player->GetY();
	double x2 = x0 + 1.0 * cos(angle), y2 = y0 - 1.0*sin(angle);
	double k = (y2 - y0) / (x2 - x0), b = -x0*k + y0;
	if (abs(x2 - x0) <= EPS) {
		double Ds = r*r - (y0 - y)*(y0 - y);
		if (Ds >= 0.0) {
			double cross_x1 = sqrt(Ds) + x, cross_y1 = y0;
			double cosv = ((cross_x1 - x0)*(x2 - x0) + (cross_y1 - y0)*(y2 - y0)) / (sqrt((cross_x1 - x0)*(cross_x1 - x0) + (cross_y1 - y0)*(cross_y1 - y0)));
			if (cosv >= -EPS) {
				return true;
			}
			return false;
		}
	}
	else {
		double Ds = (x - k*(b - y))*(x - k*(b - y)) - (1 + k*k)*(x*x + (b - y)*(b - y) - r*r);
		if (Ds >= 0.0) {
			double cross_x1 = ((x - k*(b - y)) + sqrt(Ds)) / (1 + k*k), cross_y1 = k*cross_x1 + b;
			double cosv = ((cross_x1 - x0)*(x2 - x0) + (cross_y1 - y0)*(y2 - y0)) / (sqrt((cross_x1 - x0)*(cross_x1 - x0) + (cross_y1 - y0)*(cross_y1 - y0)));
			if (cosv >= -EPS) {
				return true;
			}
			return false;
		}
	}
	return false;
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
	/*inputs.at(0, 0) = me->GetX();
	inputs.at(0, 1) = me->GetY();
	int i = 1;
	for each (auto cur in w->GetFood())
	{
		if (i == FOOD_COUNT)
			break;

		inputs.at(0, 2*i) = cur->GetX();
		inputs.at(0, 2*i + 1) = cur->GetY();
	}



	return;*/
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
		inputs.at(0, eye) = min_dist / (640 * 480);
		inputs.at(0, eye + 1) = cur_type;

		eye += 2;
		angle += step_angle;
	}


}
NeuroNet::Matrix2d inputs(1, INPUT_NEURON_COUNT);
NeuroNet::Matrix2d last_inputs;

int tick = -1;
NeuroNet::Matrix2d lastQ;
double last_full = 0.0;
void MyPlayer::Move()
{
	//try use context in training!!! 
	////debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;

	tick++;

	/////////////////////////////////////////////////////////////////////
	///////////////////////// InitFill input vector /////////////////////////
	last_inputs = inputs;
	//inputs.Fill(0.0);
	setInput(inputs, this, SENSOR_COUNT, GetWorld());
	/////////////////////////////////////////////////////////////////////
	int rnd = NeuroNet::myrand();
	int action = rnd%Actions::COUNT;
	double reward = this->GetFullness() - last_full;
	NeuroNet::Matrix2d Context = net.GetContext();
	last_full = this->GetFullness();
	////debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;
	if (!FirstStep)
	{
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
			while (Epoch--)
			{
				double error;
				if ((error = net.RMSTraining(TrainingSet)) < TRAIN_EPS)
					break;
			}
		}
	}

	NeuroNet::Matrix2d Q;
	if ((rnd % 28) < 5)
	{
		action = rnd % Actions::COUNT;
		lastAction = action;
	}
	else
	{
		net.SetContext(Context);
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
	if (tick > 0 && tick % 5000 == 0)
	{
		////debug << ERRORDEF << " " << std::string(__FILE__) << "(" << __LINE__ << "):" << std::string(__FUNCTION__) << std::endl;
		//debug << std::endl << std::endl << "========================================================" << tick << "========================================================" << std::endl;
		net.PrintFullInfo(debug);
		debug.flush();
	}
}