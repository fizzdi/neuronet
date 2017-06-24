#include "DirectXRender.h"
#include "Global.h"

LPDIRECT3D9			d3d;    // the pointer to our Direct3D interface
LPDIRECT3DDEVICE9	d3ddev;    // the pointer to the device class

using namespace std;

// this function initializes and prepares Direct3D for use
void initD3D(HWND hWnd)
{
	d3d = Direct3DCreate9(D3D_SDK_VERSION);    // create the Direct3D interface

	D3DPRESENT_PARAMETERS d3dpp;    // create a struct to hold various device information

	ZeroMemory(&d3dpp, sizeof(d3dpp));    // clear out the struct for use
	d3dpp.Windowed = TRUE;    // program windowed, not fullscreen
	d3dpp.SwapEffect = D3DSWAPEFFECT_DISCARD;    // discard old frames
	d3dpp.hDeviceWindow = hWnd;    // set the window to be used by Direct3D

								   // create a device class using this information and the info from the d3dpp stuct
	d3d->CreateDevice(D3DADAPTER_DEFAULT,
		D3DDEVTYPE_HAL,
		hWnd,
		D3DCREATE_SOFTWARE_VERTEXPROCESSING,
		&d3dpp,
		&d3ddev
	);
}

void RenderCircle(LPDIRECT3DDEVICE9 d3ddev, int X, int Y, int R, D3DCOLOR Color) {

	const int CIRCLE_RESOLUTION = 64;
	const DWORD FVF = D3DFVF_XYZRHW | D3DFVF_DIFFUSE;

	struct VERTEX_2D_DIF { // transformed colorized
		float x, y, z, rhw;
		D3DCOLOR color;
	};

	VERTEX_2D_DIF verts[CIRCLE_RESOLUTION + 1];

	for (int i = 0; i < CIRCLE_RESOLUTION + 1; i++)
	{
		verts[i].x = X + R*cos(M_PI*(i / (CIRCLE_RESOLUTION / 2.0f)));
		verts[i].y = Y + R*sin(M_PI*(i / (CIRCLE_RESOLUTION / 2.0f)));
		verts[i].z = 0;
		verts[i].rhw = 1;
		verts[i].color = Color;
	}
	d3ddev->SetFVF(FVF);
	d3ddev->DrawPrimitiveUP(D3DPT_TRIANGLEFAN, CIRCLE_RESOLUTION - 1, &verts, sizeof(VERTEX_2D_DIF));
}

void RenderLine(LPDIRECT3DDEVICE9 d3ddev, int X1, int Y1, int X2, int Y2, D3DCOLOR Color) {
	LPD3DXLINE line;
	D3DXCreateLine(d3ddev, &line);
	D3DXVECTOR2 lines[] = { D3DXVECTOR2(float(X1), float(Y1)), D3DXVECTOR2(float(X2), float(Y2)) };
	line->Begin();
	line->Draw(lines, 2, Color);
	line->End();
	line->Release();
}

void DrawTextString(LPDIRECT3DDEVICE9 d3ddev, RECT *TextRect, D3DCOLOR Color, const wchar_t *Str) {
	LPD3DXFONT pFont = 0;
	D3DXCreateFont(d3ddev, 20, 0, FW_BOLD, 0, FALSE, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH | FF_DONTCARE, TEXT("Arial"), &pFont);
	pFont->DrawText(NULL, Str, -1, TextRect, 0, Color);
	pFont->Release();
}


void RenderTexture(LPDIRECT3DDEVICE9 d3ddev, int X, int Y, const wchar_t TexturePath[]) {
	IDirect3DTexture9* g_pTexture;
	D3DXIMAGE_INFO		imageInfo;
	D3DXCreateTextureFromFileEx(
		d3ddev, TexturePath, 
		D3DX_DEFAULT, D3DX_DEFAULT, 1, 0, D3DFMT_UNKNOWN,
		D3DPOOL_MANAGED, D3DX_FILTER_NONE, D3DX_FILTER_NONE, 0,
		&imageInfo,
		NULL, 
		&g_pTexture
	);
	LPD3DXSPRITE Sprite;
	D3DXCreateSprite(d3ddev, &Sprite);

	Sprite->Begin(D3DXSPRITE_ALPHABLEND);
	Sprite->Draw(
		g_pTexture, NULL, 
		&D3DXVECTOR3(imageInfo.Width / 2, imageInfo.Height / 2, 0.0), 
		&D3DXVECTOR3(X, Y, 0.0), 
		0xFFFFFFFF
	);
	Sprite->End();
	Sprite->Release();
	g_pTexture->Release();
}

void cleanD3D(void)
{
	d3ddev->Release();    // close and release the 3D device
	d3d->Release();    // close and release Direct3D
	////debugout << "Clean D3D" << endl;
}
