#include "MyPlayer.h"

extern "C" __declspec (dllexport) MyPlayer* GetMyPlayer()
{
	return new MyPlayer;
}