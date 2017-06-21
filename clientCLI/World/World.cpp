#include "Global.h"
#include "World.h"
//#include "MyPlayer.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <ctime>

bool	   stopFlag = false;
const int MOVE_EXCEPTION	=	1;
const int OK				=	0;

using namespace std;

DWORD WINAPI DoMove(LPVOID lpParam) {
	Player* CurPlayer = reinterpret_cast<Player*>(lpParam);
	__try {
		CurPlayer->Move();
	}
	__except (EXCEPTION_EXECUTE_HANDLER) {
		return MOVE_EXCEPTION;
	}
	return OK;
}

void World::Run() {
	srand(time(NULL));
	debugout << "Loading map" << endl;
	LoadMap("maps\\map.txt");
	debugout << "Map is loaded" << endl;
	LoadPlayers();
	debugout << "Players are loaded" << endl;
	GenerateBonuses();
	//debugout << "Bonuses are generated" << endl;
	RenderFull();
	//debugout << "First render" << endl;

	int TickCount = 500000;
	
	for (CurTime = 0; CurTime <= TickCount; ++CurTime) {
		if (stopFlag) return;
		for (int j = 0; j < Players.size(); ++j) {
			if (stopFlag) return;
			if (Players[j]->GetHealth()) {
				//debugout << "Before RunTime and TimeLimit Check" << endl;
				if (!Players[j]->RunTime && !Players[j]->TimeLimit) {
					//debugout << "Run thread" << endl;
					//debugout << "Coords: " << Players[j]->GetX() << " " << Players[j]->GetY() << endl;
					auto PlayerThread = CreateThread(0, NULL, DoMove, (LPVOID)Players[j], NULL, NULL);
					auto WaitResult = WaitForSingleObject(PlayerThread, 500000);
					DWORD PlayerExitCode = OK;
					GetExitCodeThread(PlayerThread, &PlayerExitCode);
					
					//debugout << "WaitResult: " << WaitResult << endl;
					//debugout << "PlayerExitCode: " << PlayerExitCode << endl;
					
					Players[j]->RunTime = PlayerExitCode != OK;
					Players[j]->TimeLimit = WaitResult != WAIT_OBJECT_0;
				}
				if (!Players[j]->GetFullness()) Players[j]->SetHealth(max(Players[j]->GetHealth() - 5, 0));
				else {
					if (Players[j]->GetFullness() > 250 && Players[j]->GetHealth() < 1000)
						Players[j]->SetHealth(min(Players[j]->GetHealth() + 2, 1000));
					Players[j]->SetFullness(max(Players[j]->GetFullness() - 0, 0));
					//Players[j]->SetFullness(max(Players[j]->GetFullness() - 1, 0));
				}
			}

			//debugout << "Check food" << endl;
			auto Food = GetFood();
			for (auto i = 0; i < Food.size(); ++i) {
				if (Players[j]->GetDistanceTo(Food[i]) < double(Players[j]->GetR() + Food[i]->GetR())) {
					Food[i]->SetTaken();
					if (Food[i]->GetDamage() != 0) Players[j]->SetHealth(max(Players[j]->GetHealth() - Food[i]->GetDamage(), 0));
					//else Players[j]->SetFullness(min(Players[j]->GetFullness() + Food[i]->GetFullness(), 300));
					//else Players[j]->SetFullness(min(Players[j]->GetFullness() + Food[i]->GetFullness(), INT_MAX));
				else Players[j]->SetFullness(min(Players[j]->GetFullness() + 1, INT_MAX));
				}
			}

			//debugout << "Check traps" << endl;
			auto Traps = GetTraps();
			for (auto i = 0; i < Traps.size(); ++i) {
				if (!Traps[i]->IsActive()) continue;
				if (Players[j]->GetDistanceTo(Traps[i]) < double(Players[j]->GetR() + Traps[i]->GetR())) {
					Players[j]->SetHealth(max(Players[j]->GetHealth() - Traps[i]->GetDamage(), 0));
					Traps[i]->DeActivate();
				}
			}

			//debugout << "Check cornucopia" << endl;
			auto WrldCornucopia = GetCornucopia();
			if (WrldCornucopia == nullptr) continue;
			if (Players[j]->GetDistanceTo(WrldCornucopia) < double(Players[j]->GetR() + WrldCornucopia->GetR())) {
				if (WrldCornucopia->GetHasGolograf()) {
					Players[j]->HasGolograf = true;
					WrldCornucopia->SetHasGolograf(false);
				}
			}
		}
		UpdateBonuses();
		Render();
	}
	RenderResults();
}

void World::LoadMap(std::string MapPath) {
	ifstream fmap(MapPath);
	int line_num = 0;
	while (!fmap.eof()) {
		string cur_line;
		getline(fmap, cur_line);
		if (stopFlag) {
			return;
		}
		for (int i = 0; i < cur_line.length(); ++i) {
			if (stopFlag) {
				return;
			}
			switch (cur_line[i]) {
			case '~':
			{
				Block *NewOceanBlock = new Block(OCEAN);
				NewOceanBlock->SetCoords(i, line_num);
				NewOceanBlock->SetR(30);
				Blocks.push_back(NewOceanBlock);
				break;
			}
			case '*':
			{
				Trap *NewTrap = new Trap();
				NewTrap->SetCoords(i, line_num);
				NewTrap->SetR(rand() % 10 + 10);
				NewTrap->SetDamage(rand() % 30);
				Traps.push_back(NewTrap);
				break;
			}
			case 'o':
			{
				Block *NewForestBlock = new Block(FOREST);
				NewForestBlock->SetCoords(i, line_num);
				NewForestBlock->SetR(30);
				Blocks.push_back(NewForestBlock);
				break;
			}
			case '|':
			{
				Cornucopia *NewCornucopia = new Cornucopia();
				NewCornucopia->SetCoords(i, line_num);
				NewCornucopia->SetR(30);
				CurrentCornucopia = NewCornucopia;
				break;
			}
			default:
				break;
			}
		}
		line_num++;
	}
}

void World::LoadPlayers() {
	////debugout << "Loading players" << endl;

	WIN32_FIND_DATA ffd;
	wchar_t szDir[] = L"players\\*";
	HANDLE hFind = INVALID_HANDLE_VALUE;
	vector<wstring> PlayerFiles;

	hFind = FindFirstFile(szDir, &ffd);

	if (INVALID_HANDLE_VALUE == hFind) return;
	do
	{
		if (ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) continue;
		else PlayerFiles.push_back(wstring(ffd.cFileName));
	} while (FindNextFile(hFind, &ffd) != 0);

	FindClose(hFind);

	////debugout << "Players count: " << PlayerFiles.size() << endl;

	D3DCOLOR Colors[] = { D3DCOLOR_XRGB(255,0,0), D3DCOLOR_XRGB(0,0,255), D3DCOLOR_XRGB(0, 255, 0), D3DCOLOR_XRGB(252, 131, 10) };
	pair<int, int> Coords[] = { {320, 240}, {540, 100}, {100, 380}, {540, 380} };

	for (int i = 0; i < PlayerFiles.size(); ++i) {
		auto hLib = LoadLibrary((L"players\\" + PlayerFiles[i]).c_str());
		debugout << "Loading library: " << PlayerFiles[i] << endl;
		if (hLib != NULL) 
		{
			debugout << "Player library is loaded" << endl;
			auto GetMyPlayer = (GetMyPlayerFunction)GetProcAddress((HMODULE)hLib, "GetMyPlayer");
			if (GetMyPlayer != NULL) 
			{
				Player *NewPlayer = GetMyPlayer();
				NewPlayer->SetCoords(Coords[i].first, Coords[i].second);
				NewPlayer->SetR(20);
				NewPlayer->SetWorld(this);
				//NewPlayer->SetSpeed(10);
				NewPlayer->SetSpeed(20);
				NewPlayer->SetAngle(M_PI*0.5);
				//NewPlayer->SetFullness(300);
				NewPlayer->SetFullness(1);
				//NewPlayer->SetHealth(1000);
				NewPlayer->SetHealth(INT_MAX);
				NewPlayer->SetColor(Colors[i]);
				//debugout << "Player " << PlayerFiles[i] << " library begin init" << endl;
				NewPlayer->Init();
				//debugout << "Player " << PlayerFiles[i] << " library after init" << endl;
				NewPlayer->Idx = i;
				//NewPlayer->PlayerLib = hLib;
				Players.push_back(NewPlayer);
			}
		}
	}
	////debugout << "Players are loaded: " << Players.size() << endl;
}

void World::UnloadPlayers() {
	for (int i = 0; i < Players.size(); ++i) {
		if (Players[i]->PlayerLib != nullptr && !Players[i]->TimeLimit) {
			FreeLibrary(Players[i]->PlayerLib);
			Players[i]->PlayerLib = nullptr;
		}
	}
}


void World::GenerateFood() {
	int FoodCount = rand() % 30;
	for (int i = 0; i < FoodCount; ++i) {
		Food *NewFood = new Food();
		NewFood->SetCoords(rand() % 600 + 20, rand() % 425 + 20);
		NewFood->SetR(15);
		NewFood->SetLifetime(rand() % 200 + 10);
		if (rand() % 20 == 1) NewFood->SetDamage(rand() % 500);
		else NewFood->SetFullness(30);
		WorldFood.push_back(NewFood);
	}
}

void World::GenerateWeapon() {
	int WeaponCount = rand() % 2;

	for (int i = 0; i < WeaponCount; ++i) {
		if (rand() % 10 > 2) {
			Knife *NewKnife = new Knife();
			NewKnife->SetCoords(rand() % 600 + 20, rand() % 425 + 20);
			NewKnife->SetDamage(20);
			NewKnife->SetLifetime(rand() % 200 + 10);
			NewKnife->SetR(5);
			WorldWeapon.push_back(NewKnife);
		}
		else {
			Bow *NewBow = new Bow();
			NewBow->SetCoords(rand() % 640, rand() % 480);
			NewBow->SetDamage(50);
			NewBow->SetLifetime(rand() % 10 + 10);
			NewBow->SetSpeed(40);
			NewBow->SetR(5);
			WorldWeapon.push_back(NewBow);
		}

	}
}

void World::GenerateBonuses() {
	GenerateFood();
	//GenerateWeapon();
}

void World::UpdateBonuses() {
	//debugout << "Update Bonuses" << endl;

	for (auto it = WorldFood.begin(); it != WorldFood.end(); ) {
		if (!(*it)->IsTaken()) {
			(*it)->SetLifetime((*it)->GetLifeTime() - 1);
			if ((*it)->GetLifeTime() == 0) {
				delete (*it);
				it = WorldFood.erase(it);
			}
			else ++it;
		}
		else {
			delete (*it);
			it = WorldFood.erase(it);
		}
	}

	for (auto it = WorldWeapon.begin(); it != WorldWeapon.end(); ) {
		if (!(*it)->IsTaken()) {
			(*it)->SetLifetime((*it)->GetLifeTime() - 1);
			if ((*it)->GetLifeTime() == 0) {
				delete (*it);
				it = WorldWeapon.erase(it);
			}
			else ++it;
		}
		else {
			delete (*it);
			it = WorldWeapon.erase(it);
		}
	}

	if (WorldFood.size() <= 30) GenerateFood();
	//if (WorldWeapon.size() <= 2) GenerateWeapon();
	//debugout << "End Update Bonuses" << endl;
}

void World::Render() {
	////debugout << "Render start" << endl;
	d3ddev->Clear(0, NULL, D3DCLEAR_TARGET, D3DCOLOR_XRGB(10, 10, 10), 1.0f, 0);

	d3ddev->BeginScene();    // begins the 3D scene

	////debugout << "Render objects" << endl;
	for (int i = 0; i < Blocks.size(); ++i) Blocks[i]->Render(d3ddev);
	////debugout << "Render blocks" << endl;
	for (int i = 0; i < Traps.size(); ++i) Traps[i]->Render(d3ddev);
	////debugout << "Render traps" << endl;
	if (CurrentCornucopia != nullptr) CurrentCornucopia->Render(d3ddev);
	////debugout << "Render cornucopia" << endl;
	for (int i = 0; i < Players.size(); ++i) Players[i]->Render(d3ddev);
	////debugout << "Render players" << endl;
	for (int i = 0; i < WorldFood.size(); ++i) WorldFood[i]->Render(d3ddev);
	////debugout << "Render food" << endl;
	for (int i = 0; i < WorldWeapon.size(); ++i) WorldWeapon[i]->Render(d3ddev);
	////debugout << "Render weapon" << endl;

	wchar_t CurMessage[1000];
	RECT CurCoord;
	int i;
	for (i = 0; i < Players.size(); ++i) {
		wstring CurName = Players[i]->GetName();
		////debugout << "Name: " << CurName << endl;
		CurCoord = { 0, i * 20, 500, i * 20 + 20 };
		wsprintf(CurMessage, L"%s: H%d F%d", CurName.c_str(), Players[i]->GetHealth(), Players[i]->GetFullness());
		DrawTextString(d3ddev, &CurCoord, Players[i]->GetColor(), CurMessage);
	}

	//CurCoord = { GetWidth() / 2 - 50, GetHeight() - 100, GetWidth() / 2 + 50, GetHeight() };
	CurCoord = { GetWidth() - 50 - 10, GetHeight() - 70, GetWidth() - 10, GetHeight() - 50 };
	wsprintf(CurMessage, L"%d", CurTime);
	DrawTextString(d3ddev, &CurCoord, D3DCOLOR_XRGB(255, 255, 255), CurMessage);

	//RenderTexture(d3ddev, 0, 0);

	d3ddev->EndScene();    // ends the 3D scene

	d3ddev->Present(NULL, NULL, NULL, NULL);   // displays the created frame on the screen
	////debugout << "Render finish" << endl;
}

void World::RenderFull() {
	//debugout << "Render start" << endl;
	d3ddev->Clear(0, NULL, D3DCLEAR_TARGET, D3DCOLOR_XRGB(255, 255, 255), 1.0f, 0);

	d3ddev->BeginScene();    // begins the 3D scene

	//debugout << "Render objects" << endl;
	for (int i = 0; i < Blocks.size(); ++i) Blocks[i]->RenderFull(d3ddev);
	//debugout << "Render blocks" << endl;
	for (int i = 0; i < Traps.size(); ++i) Traps[i]->RenderFull(d3ddev);
	//debugout << "Render traps" << endl;
	if (CurrentCornucopia != nullptr) CurrentCornucopia->RenderFull(d3ddev);
	//debugout << "Render cornucopia" << endl;
	for (int i = 0; i < Players.size(); ++i) {
		Players[i]->RenderFull(d3ddev);
	}
	////debugout << "Render players" << endl;
	
	for (int i = 0; i < WorldFood.size(); ++i) WorldFood[i]->RenderFull(d3ddev);
	////debugout << "Render food" << endl;
	for (int i = 0; i < WorldWeapon.size(); ++i) WorldWeapon[i]->RenderFull(d3ddev);
	////debugout << "Render weapon" << endl;

	wchar_t CurMessage[1000];
	RECT CurCoord;
	int i;
	for (i = 0; i < Players.size(); ++i) {
		wstring CurName = Players[i]->GetName();
		//debugout << "Name: " << CurName << endl;
		CurCoord = { 0, i * 20, 500, i * 20 + 20 };
		if (Players[i]->TimeLimit) {
			wsprintf(CurMessage, L"%s*: H%d F%d", CurName.c_str(), Players[i]->GetHealth(), Players[i]->GetFullness());
		}
		else if (Players[i]->RunTime) {
			wsprintf(CurMessage, L"%s**: H%d F%d", CurName.c_str(), Players[i]->GetHealth(), Players[i]->GetFullness());
		}
		else {
			wsprintf(CurMessage, L"%s: H%d F%d", CurName.c_str(), Players[i]->GetHealth(), Players[i]->GetFullness());
		}
		DrawTextString(d3ddev, &CurCoord, Players[i]->GetColor(), CurMessage);
	}

	CurCoord = { GetWidth() - 50 - 10, GetHeight() - 70, GetWidth() - 10, GetHeight() - 50 };
	wsprintf(CurMessage, L"%d", CurTime);
	DrawTextString(d3ddev, &CurCoord, D3DCOLOR_XRGB(0, 0, 0), CurMessage);
	
	d3ddev->EndScene();    // ends the 3D scene

	d3ddev->Present(NULL, NULL, NULL, NULL);   // displays the created frame on the screen
	//debugout << "Render finish" << endl;
}

void World::RenderResults() {
	d3ddev->Clear(0, NULL, D3DCLEAR_TARGET, D3DCOLOR_XRGB(255, 255, 255), 1.0f, 0);

	d3ddev->BeginScene();    // begins the 3D scene

	wchar_t CurMessage[1000];
	RECT CurCoord;
	int i;
	sort(Players.begin(), Players.end(), [](Player* CurPlayer1, Player* CurPlayer2) {
		return (CurPlayer1->GetFullness() + CurPlayer1->GetHealth()) >= (CurPlayer2->GetFullness() + CurPlayer2->GetHealth());
	});
	for (i = 0; i < Players.size(); ++i) {
		wstring CurName = Players[i]->GetName();
		CurCoord = { 320 - 100, i * 20 + 100, 320 + 100, i * 20 + 100 + 20 };
		wsprintf(CurMessage, L"%d. %s: %d", i + 1, CurName.c_str(), Players[i]->GetHealth() + Players[i]->GetFullness());
		DrawTextString(d3ddev, &CurCoord, Players[i]->GetColor(), CurMessage);
	}

	d3ddev->EndScene();    // ends the 3D scene

	d3ddev->Present(NULL, NULL, NULL, NULL);   // displays the created frame on the screen
}


void World::RenderBasicWorldField() {

}





