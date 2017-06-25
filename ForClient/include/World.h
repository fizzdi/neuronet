#ifndef _MY_WORLD
#define _MY_WORLD

#pragma once
#include "windows.h"

#include "Player.h"
#include "Unit.h"
#include "Bonus.h"

#include <vector>
#include <string>
#include <fstream>

class Player;

typedef Player* PPlayer;
typedef PPlayer(*GetMyPlayerFunction) ();

class  World {
private:
	std::vector<Player*> Players;
	std::vector<Block*> Blocks;
	std::vector<Trap *> Traps;
	Cornucopia *CurrentCornucopia;
	std::vector<Food*> WorldFood;
	std::vector<Weapon*> WorldWeapon;
	int Height, Width;
	int CurTime;

	void LoadMap(std::string MapPath);
	void SetSize(int NewHeight, int NewWidth) { Height = NewHeight; Width = NewWidth; };
	void LoadPlayers();
	void UnloadPlayers();
	void Render();
	void RenderFull();
	void RenderResults();

	void RenderBasicWorldField();

	void Run();
	
	void GenerateBonuses();
	void GenerateFood();
	void GenerateWeapon();

	void UpdateBonuses();
public:
	World() {
		CurrentCornucopia = nullptr;		
	};
	int GetHeight() { return Height - 20; };
	int GetWidth() { return Width; };
	std::vector<Player*> GetPlayers() { return Players; };
	std::vector<Block*> GetBlocks() { return Blocks; };
	std::vector<Trap*> GetTraps() { return Traps; };
	Cornucopia* GetCornucopia() { return CurrentCornucopia; };
	std::vector<Food*> GetFood() { return WorldFood; };
	std::vector<Weapon*> GetWeapon() { return WorldWeapon; };
	int GetCurrentTime() { return CurTime; };
	

	friend void RunWorld(World *CurWorld, int Width, int Height) {
		CurWorld->SetSize(Height, Width);
		CurWorld->Run();
	}

	void FreeObjects() {
		UnloadPlayers();
		for (int i = 0; i < Blocks.size(); ++i) if (Blocks[i] != nullptr) {
			delete Blocks[i];
			Blocks[i] = nullptr;
		}
		for (int i = 0; i < Players.size(); ++i) if (Players[i] != nullptr) {
			delete Players[i];
			Players[i] = nullptr;
		}
		for (int i = 0; i < Traps.size(); ++i) if (Traps[i] != nullptr) {
			delete Traps[i];
			Traps[i] = nullptr;
		}
		for (int i = 0; i < WorldFood.size(); ++i) if (WorldFood[i] != nullptr) {
			delete WorldFood[i];
			WorldFood[i] = nullptr;
		}
		for (int i = 0; i < WorldWeapon.size(); ++i) if (WorldWeapon[i] != nullptr) {
			delete WorldWeapon[i];
			WorldWeapon[i] = nullptr;
		}
		if (CurrentCornucopia != nullptr) {
			delete CurrentCornucopia;
			CurrentCornucopia = nullptr;
		}
	};

	~World() {
		std::wofstream //debugout("debug.txt", std::ios_base::app);
		//debugout << "~World" << std::endl;
		FreeObjects();
	}
};

#endif


