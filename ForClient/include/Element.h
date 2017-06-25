#ifndef _ELEMENT
#define _ELEMENT
#pragma once

#include "DirectXRender.h"

class  __declspec (dllexport) Element {
private:
	int X, Y;
	int R;
public:
	Element() {};
	int GetX() { return X; };
	int GetY() { return Y; };
	int GetR() { return R; };
	void SetCoords(int NewX, int NewY) { X = NewX; Y = NewY; };
	void SetR(int NewR) { R = NewR; };
	~Element() {};
private:
	void Render(LPDIRECT3DDEVICE9 d3ddev) {};
	void RenderFull(LPDIRECT3DDEVICE9 d3ddev) {};
};
#endif
