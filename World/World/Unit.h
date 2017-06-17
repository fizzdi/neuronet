#pragma once

#include "Element.h"

#define UNKNOWN_UNIT		0
#define BLOCK_UNIT			1
#define TRAP_BLOCK			2
#define CORNUCOPIA_BLOCK	3

#define UNKNOWN_BLOCK		0
#define OCEAN				1
#define FOREST				2

class  Unit : public Element {
	friend class World;
	friend class Player;
private:
	int Type;
public:
	Unit() { Type = UNKNOWN_UNIT; };
	int GetType() { return Type; };
protected:
	void SetType(int NewType) { Type = NewType; };
};

class  Block : public Unit {
	friend class World;
	friend class Player;
private:
	int BlockType;

	void Render(LPDIRECT3DDEVICE9 d3ddev);
	void RenderFull(LPDIRECT3DDEVICE9 d3ddev);
public:
	Block() { SetType(BLOCK_UNIT); BlockType = UNKNOWN_BLOCK; };
	Block(int CurBlockType) { SetType(BLOCK_UNIT); BlockType = CurBlockType; };
	int GetBlockType() { return BlockType; };
};

class  Trap : public Unit {
	friend class World;
	friend class Player;
private:
	int Damage;
	int LifeTime;
	bool Active;

	void Render(LPDIRECT3DDEVICE9 d3ddev);
	void RenderFull(LPDIRECT3DDEVICE9 d3ddev);
public:
	Trap() { SetType(TRAP_BLOCK); LifeTime = 0; Active = true; Damage = 0; };
	Trap(int NewLifeTime, bool NewIsActive) { SetType(TRAP_BLOCK); LifeTime = NewLifeTime; Active = NewIsActive; }
	bool IsActive() { return Active; };
	int GetDamage() { return Damage; };
private:
	void Activate() { Active = true; };
	void DeActivate() { Active = false; };
	void SetDamage(int NewDamage) { Damage = NewDamage; };
};


class  Cornucopia : public Unit {
	friend class World;
	friend class Player;
private:
	bool HasGolograf;
	int KnifeCount;
	int BowCount;

	void SetHasGolograf(bool NewHasGolograf) { HasGolograf = NewHasGolograf; };
	void SetKnifeCount(bool NewKnifeCount) { KnifeCount = NewKnifeCount; };
	void SetBowCount(bool NewBowCount) { BowCount = NewBowCount; };

	void Render(LPDIRECT3DDEVICE9 d3ddev);
	void RenderFull(LPDIRECT3DDEVICE9 d3ddev);
public:
	Cornucopia() {
		SetType(CORNUCOPIA_BLOCK); KnifeCount = BowCount = 0; HasGolograf = true;
	};
	bool GetHasGolograf() { return HasGolograf; };
	bool HasKnife() { return KnifeCount != 0; };
	bool HasBow() { return BowCount != 0; };
};
