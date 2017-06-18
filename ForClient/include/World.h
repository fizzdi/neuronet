#pragma once
#include "windows.h"

#include "Player.h"
#include "Unit.h"
#include "Bonus.h"

#include <vector>
#include <string>

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
	void Render();
	void RenderFull();

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
	int GetHeight() { return Height; };
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

	~World() {
		for (int i = 0; i < Blocks.size(); ++i) delete Blocks[i];
		for (int i = 0; i < Players.size(); ++i) delete Players[i];
		for (int i = 0; i < Traps.size(); ++i) delete Traps[i];
		for (int i = 0; i < WorldFood.size(); ++i) delete WorldFood[i];
		for (int i = 0; i < WorldWeapon.size(); ++i) delete WorldWeapon[i];
		if (CurrentCornucopia != nullptr) {
			delete CurrentCornucopia;
			CurrentCornucopia = nullptr;
		}
	}
};

