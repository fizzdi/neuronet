#pragma once

#include "Element.h"

#define UNKNOWN_BONUS	0
#define FOOD			1
#define WEAPON			2

#define UNKNOWN_WEAPON_TYPE 0
#define KNIFE				1
#define BOW					2

class  Bonus : public Element {
	friend class World;
	friend class Player;
private:
	int Type;
	int Lifetime;
	int Damage;
	bool Taken;

	virtual void Render(LPDIRECT3DDEVICE9 d3ddev) {};
	virtual void RenderFull(LPDIRECT3DDEVICE9 d3ddev) {};
	
	bool IsTaken() { return Taken; };
	void SetTaken() { Taken = true; }
public:
	Bonus() { Type = UNKNOWN_BONUS; Damage = Lifetime = 0; Taken = false; };
	int GetType() { return Type; };
	int GetLifeTime() { return Lifetime; };
	virtual int GetDamage() { return Damage; };
protected:
	void SetType(int NewType) { Type = NewType; };
	void SetLifetime(int NewLifetime) { Lifetime = NewLifetime; };
	void SetDamage(int NewDamage) { Damage = NewDamage; };
};

class  Food : public Bonus {
	friend class World;
	friend class Player;
private:
	int Fullness;

	void SetFullness(int NewFullness) { Fullness = NewFullness; };

	void Render(LPDIRECT3DDEVICE9 d3ddev);
	void RenderFull(LPDIRECT3DDEVICE9 d3ddev);
public:
	Food() {
		SetType(FOOD); Fullness = 0; 
	};
	int GetFullness() { return Fullness; };
private:
	using Bonus::GetDamage;
};

class  Weapon : public Bonus {
	friend class World;
	friend class Player;
private:
	int Speed;
	int WeaponType;

	void SetSpeed(int NewSpeed) { Speed = NewSpeed; };

	virtual void Render(LPDIRECT3DDEVICE9 d3ddev) {};
	virtual void RenderFull(LPDIRECT3DDEVICE9 d3ddev) {};
public:
	Weapon() {
		SetType(WEAPON); WeaponType = UNKNOWN_WEAPON_TYPE;
	};
	int GetSpeed() { return Speed; };
	int GetWeaponType() { return WeaponType; };
protected:
	void SetWeaponType(int NewWeaponType) { WeaponType = NewWeaponType; };
};

class  Knife : public Weapon {
	friend class World;
	friend class Player;
public:
	Knife() { SetWeaponType(KNIFE); };
private:
	void Render(LPDIRECT3DDEVICE9 d3ddev);
	void RenderFull(LPDIRECT3DDEVICE9 d3ddev);
};

class  Bow : public Weapon {
	friend class World;
	friend class Player;
private:
	int ArrowCount;

	void SetArrowCount(int NewArrowCount) { ArrowCount = NewArrowCount; };

	void Render(LPDIRECT3DDEVICE9 d3ddev);
	void RenderFull(LPDIRECT3DDEVICE9 d3ddev);
public:
	Bow() { SetWeaponType(BOW); ArrowCount = 0; };
	int GetArrowCount() { return ArrowCount; };
};


