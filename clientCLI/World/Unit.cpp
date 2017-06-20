#include "Unit.h"

void Block::Render(LPDIRECT3DDEVICE9 d3ddev) {
	D3DCOLOR Color = D3DCOLOR_XRGB(155, 155, 155); //Grey
	if (GetBlockType() == OCEAN) Color = D3DCOLOR_XRGB(18, 215, 208); // Sky
	else if (GetBlockType() == FOREST) Color = D3DCOLOR_XRGB(8, 233, 68); // Green
	RenderCircle(d3ddev, GetX(), GetY(), GetR(), Color);
}

void Block::RenderFull(LPDIRECT3DDEVICE9 d3ddev) {
	if (GetBlockType() == OCEAN) RenderTexture(d3ddev, GetX(), GetY(), L"lake.png"); // Sky
	else if (GetBlockType() == FOREST) RenderTexture(d3ddev, GetX(), GetY(), L"forest.png"); // Green
}

void Trap::Render(LPDIRECT3DDEVICE9 d3ddev) {
	D3DCOLOR Color = D3DCOLOR_XRGB(246, 15, 100); //Red
	if (!IsActive()) Color = D3DCOLOR_XRGB(94, 84, 236); //Blue
	RenderCircle(d3ddev, GetX(), GetY(), GetR(), Color);
}

void Trap::RenderFull(LPDIRECT3DDEVICE9 d3ddev) {
	if (IsActive()) RenderTexture(d3ddev, GetX(), GetY(), L"trap.png");
	else RenderTexture(d3ddev, GetX(), GetY(), L"inactive_trap.png");
}

void Cornucopia::Render(LPDIRECT3DDEVICE9 d3ddev) {
	D3DCOLOR Color = D3DCOLOR_XRGB(246, 192, 15); //Yellow
	RenderCircle(d3ddev, GetX(), GetY(), GetR(), Color);
}

void Cornucopia::RenderFull(LPDIRECT3DDEVICE9 d3ddev) {
	RenderTexture(d3ddev, GetX(), GetY(), L"cornucopia.png");
}