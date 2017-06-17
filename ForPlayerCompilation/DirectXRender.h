#pragma once

#include <d3d9.h>
#include "d3dx9.h"
#include <cmath>

#define M_PI       3.14159265358979323846   // pi

// include the Direct3D Library file
#pragma comment (lib, "d3d9.lib")
#pragma comment (lib, "d3dx9.lib")

void RenderCircle(LPDIRECT3DDEVICE9 d3ddev, int X, int Y, int R, D3DCOLOR Color);
void DrawTextString(LPDIRECT3DDEVICE9 d3ddev, RECT *TextRect, D3DCOLOR Color, const wchar_t *Str);
void RenderTexture(LPDIRECT3DDEVICE9 d3ddev, int X, int Y, const wchar_t TexturePath[]);
