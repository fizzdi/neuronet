#include "Global.h"
#include "DirectXRender.h"
#include "World.h"
#include <ctime>

DWORD WINAPI ThreadProc(LPVOID lpParam);
BOOL RegisterDLLWindowClass(wchar_t szClassName[]);
LRESULT CALLBACK DLLWindowProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam);
DWORD WINAPI DoAll(LPVOID lpParam);

HINSTANCE  inj_hModule;          //Injected Modules Handle
HWND       hWnd;
HANDLE	   hRunThread;

using namespace std;

std::wofstream debugout("debug.txt");

BOOL APIENTRY DllMain(HMODULE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved)
{
	srand(time(nullptr));
	if (ul_reason_for_call == DLL_PROCESS_ATTACH) {
		//debugout << "Start" << endl;
		inj_hModule = hModule;
		CreateThread(0, NULL, ThreadProc, (LPVOID)L"Code hunger games challenge", NULL, NULL);
	}
	return TRUE;
}

BOOL RegisterDLLWindowClass(wchar_t szClassName[])
{
	WNDCLASSEX wc;
	wc.hInstance = inj_hModule;
	wc.lpszClassName = (LPCWSTR)L"CodeGameChallengeClass";
	wc.lpszClassName = (LPCWSTR)szClassName;
	wc.lpfnWndProc = DLLWindowProc;
	wc.style = CS_DBLCLKS;
	wc.cbSize = sizeof(WNDCLASSEX);
	wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	wc.hIconSm = LoadIcon(NULL, IDI_APPLICATION);
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);
	wc.lpszMenuName = NULL;
	wc.cbClsExtra = 0;
	wc.cbWndExtra = 0;
	wc.hbrBackground = (HBRUSH)COLOR_BACKGROUND;
	if (!RegisterClassEx(&wc))
		return 0;
}

DWORD WINAPI ThreadProc(LPVOID lpParam)
{
	//debugout << "ThreadProc" << endl;
	MSG messages;
	wchar_t *pString = reinterpret_cast<wchar_t *> (lpParam);
	RegisterDLLWindowClass(L"CodeGameChallengeClass");
	HWND hwnd = CreateWindowEx(
		0,
		L"CodeGameChallengeClass", pString,
		WS_VISIBLE | WS_SYSMENU | WS_MAXIMIZEBOX | WS_MINIMIZEBOX,
		CW_USEDEFAULT, CW_USEDEFAULT,
		640, 480, NULL, NULL, inj_hModule, NULL
	);
	hWnd = hwnd;
	ShowWindow(hwnd, SW_SHOWNORMAL);

	initD3D(hwnd);

	hRunThread = CreateThread(0, NULL, DoAll, (LPVOID)hwnd, NULL, NULL);

	while (GetMessage(&messages, NULL, 0, 0))
	{
		TranslateMessage(&messages);
		DispatchMessage(&messages);
	}
	//debugout << "End ThreadProc" << endl;
	return 1;
}

//Our new windows proc
LRESULT CALLBACK DLLWindowProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	switch (message)
	{
	case WM_DESTROY:
		//debugout << "Destroy window" << endl;
		stopFlag = true;
		if (hRunThread != INVALID_HANDLE_VALUE) {
			auto res = WaitForSingleObject(hRunThread, INFINITE);
			//debugout << "WaitForSingleObject result: " << res << endl;
		}
		cleanD3D();
		PostQuitMessage(0);
		break;
	default:
		return DefWindowProc(hwnd, message, wParam, lParam);
	}
	return 0;
}

DWORD WINAPI DoAll(LPVOID lpParam)
{
	World *CurWorld = new World;
	try {
		RunWorld(CurWorld, 640, 480);
	}
	catch (...) {
		//debugout << "Clean D3D - exception" << endl;
		delete CurWorld;
		return 1;
	}
	//debugout << "Delete world" << endl;
	delete CurWorld;
	return 1;
}
