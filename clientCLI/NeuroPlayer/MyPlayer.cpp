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
string elty[] = { "FOOD", "ENEMY", "BLOCK", "COUNT" };

pair<double, ElementType> getDistanceOnWall(Player *me, World *w, double angle, double step_angle)
{
	double x0 = me->GetX();
	double y0 = me->GetY();
	double k = -tan(angle);
	double angleB;

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
	double angle = M_PI/100;

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
NeuroNet::Matrix2d last_inputs;

int tick = -1;
void MyPlayer::Move()
{
	tick++;

	/////////////////////////////////////////////////////////////////////
	///////////////////////// InitFill input vector /////////////////////////
	last_inputs = inputs;
	inputs.Fill(0.0);
	setInput(inputs, this, SENSOR_COUNT, GetWorld());
	/////////////////////////////////////////////////////////////////////
	//debug << inputs << endl;
	vr.at(0, 0) = GetFullness();
	int action = -1;
	int rnd = NeuroNet::myrand();

	if (!FirstStep)
	{
		int test = TrainingSet.size();
		nets[lastAction].AddTest(TrainingSet, last_inputs, vr);
		if (tick % TRAIN_PERIOD == 0)
		{
			int Epoch = TRAIN_EPOCH;
			debug << "LA " << lastAction << endl;
			while (Epoch--)
			{
				double error;
				if ((error = nets[lastAction].RMSTraining(TrainingSet)) < TRAIN_EPS)
					//if (nets[lastAction].RunTrainingSetOffline(TrainingSet) < TRAIN_EPS)
					break;
				debug << "tick " << tick  << " ideal " << vr.at(0,0) << " error: " << error << endl;
			}
		}
	}

	if (tick % RANDOM_ACTION_PERIOD == 0)
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
			debug << curQ << endl;
			if (Q < curQ)
			{
				Q = curQ;
				action = i;
			}
		}
	}
	debug <<"SELECT " << Q << " i " << action << endl;

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