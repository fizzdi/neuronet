#pragma once

#include "Element.h"
#include "World.h"
#include <vector>

class World;

class  __declspec (dllexport) Player : public Element {
	friend class World;
private:
	std::wstring Name;
	double Angle;
	int Speed;
	int RotateSpeed;
	int Health;
	int Fullness; //сытость
	std::vector<Player*> Enemies;
	bool HasGolograf;
	World* CurrentWorld;
	D3DCOLOR Color;
	int Idx;
private:
	void Render(LPDIRECT3DDEVICE9 d3ddev);
	void RenderFull(LPDIRECT3DDEVICE9 d3ddev);

	void SetSpeed(int NewSpeed) { Speed = NewSpeed; };
	void SetRotateSpeed(int NewRotateSpeed) { RotateSpeed = NewRotateSpeed; };
	void SetHealth(int NewHealth) { Health = NewHealth; };
	void SetFullness(int NewFullness) { Fullness = NewFullness; };
	void SetAngle(double NewAngle) { Angle = NewAngle; };
	void AddEnemy(Player *EP) { Enemies.push_back(EP); };
	void ClearEnemies() { Enemies.clear(); };
	void SetHasGolograf(bool NewHasGolograf) { HasGolograf = NewHasGolograf; };
	void SetWorld(World *NewWorld) { CurrentWorld = NewWorld; };
	D3DCOLOR GetColor() { return Color; };
	void SetColor(D3DCOLOR NewColor) { Color = NewColor; };
public:
	Player() { Speed = 0; Health = 0; Fullness = 0; HasGolograf = false; };
	~Player() {};
	int GetSpeed() { return Speed; };
	int GetHealth() { return Health; };
	int GetFullness() { return Fullness; };
	bool GetHasGolograf() { return HasGolograf; };
	World* GetWorld() { return CurrentWorld; };
	std::wstring GetName() { return Name; };
	double GetAngle() { return Angle; };
	void SetName(std::wstring NewName) { Name = NewName; };

	void Rotate(double RotationAngle);
	
	double GetAngleTo(int ToX, int ToY);
	double GetAngleTo(Element *SomeObject);

	double GetDistanceTo(int ToX, int ToY);
	double GetDistanceTo(Element *SomeObject);

	void MoveTo(int ToX, int ToY);
	void MoveTo(Element *SomeObject);
	void Strike(Player *EP);

	virtual void Init() {};
	virtual void Move() {};
};


