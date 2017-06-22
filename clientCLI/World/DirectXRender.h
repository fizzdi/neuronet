#ifndef _RENDERER
#define _RENDERER
#pragma once

#include <d3d9.h>
#include "d3dx9.h"
#include <cmath>

#define M_PI       3.14159265358979323846   // pi

// include the Direct3D Library file
#pragma comment (lib, "d3d9.lib")
#pragma comment (lib, "d3dx9.lib")

extern LPDIRECT3D9			d3d;    // the pointer to our Direct3D interface
extern LPDIRECT3DDEVICE9	d3ddev;    // the pointer to the device class

void initD3D(HWND hWnd);

void RenderCircle(LPDIRECT3DDEVICE9 d3ddev, int X, int Y, int R, D3DCOLOR Color);
void DrawTextString(LPDIRECT3DDEVICE9 d3ddev, RECT *TextRect, D3DCOLOR Color, const wchar_t *Str);
void RenderLine(LPDIRECT3DDEVICE9 d3ddev, int X1, int Y1, int X2, int Y2);
void RenderTexture(LPDIRECT3DDEVICE9 d3ddev, int X, int Y, const wchar_t TexturePath[]);

void cleanD3D(void);

#endif
