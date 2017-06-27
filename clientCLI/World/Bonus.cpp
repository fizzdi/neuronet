#include "Bonus.h"

void Food::Render(LPDIRECT3DDEVICE9 d3ddev) {
	D3DCOLOR Color = D3DCOLOR_XRGB(204, 204, 255); //Orange
	RenderCircle(d3ddev, GetX(), GetY(), GetR(), Color);
}

void Food::RenderFull(LPDIRECT3DDEVICE9 d3ddev) {
	if (GetDamage() != 0 ) RenderTexture(d3ddev, GetX(), GetY(), L"poison.png");
	else RenderTexture(d3ddev, GetX(), GetY(), L"food.png");
}

void Knife::Render(LPDIRECT3DDEVICE9 d3ddev) {
	D3DCOLOR Color = D3DCOLOR_XRGB(193, 204, 205); //Light grey
	RenderCircle(d3ddev, GetX(), GetY(), GetR(), Color);
}

void Knife::RenderFull(LPDIRECT3DDEVICE9 d3ddev) {
	RenderTexture(d3ddev, GetX(), GetY(), L"knife.png");
}

void Bow::Render(LPDIRECT3DDEVICE9 d3ddev) {
	D3DCOLOR Color = D3DCOLOR_XRGB(119, 60, 15); //Dark brown
	RenderCircle(d3ddev, GetX(), GetY(), GetR(), Color);
}

void Bow::RenderFull(LPDIRECT3DDEVICE9 d3ddev) {
	RenderTexture(d3ddev, GetX(), GetY(), L"bow.png");
}