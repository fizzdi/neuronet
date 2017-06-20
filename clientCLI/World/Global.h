#ifndef _GLOBAL
#define _GLOBAL
#pragma once

#include <windows.h>
#include <windowsx.h>
#include <fstream>

extern HWND				hWnd;
extern std::wofstream	debugout;
extern bool				stopFlag;


extern "C" __declspec (dllexport) void RunFromCLI();
#endif