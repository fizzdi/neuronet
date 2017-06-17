#include "Player.h"

#include <fstream>

using namespace std;
#define M_PI       3.14159265358979323846   // pi

void Player::Render(LPDIRECT3DDEVICE9 d3ddev) {
	RenderCircle(d3ddev, GetX(), GetY(), GetR(), Color);
}

void Player::RenderFull(LPDIRECT3DDEVICE9 d3ddev) {
	wstring CurPath[] = {L"player1.png", L"player2.png", L"player3.png", L"player4.png" };
	RenderTexture(d3ddev, GetX(), GetY(), CurPath[Idx].c_str());
}

void Player::Rotate(double RotationAngle) {
	Angle += RotationAngle;
}

double Player::GetAngleTo(int ToX, int ToY) {
	if (GetDistanceTo(ToX, ToY) <= 1e-9) return 0.0;
	auto CalculatedCos = ((ToX - GetX())*cos(Angle) - (ToY - GetY())*sin(Angle)) / GetDistanceTo(ToX, ToY);
	if (CalculatedCos > 1.0) CalculatedCos = 1.0;
	if (CalculatedCos < -1.0) CalculatedCos = -1.0;
	return acos(CalculatedCos);
}

double Player::GetAngleTo(Element *SomeObject) {
	return GetAngleTo(SomeObject->GetX(), SomeObject->GetY());
}

double Player::GetDistanceTo(int ToX, int ToY) {
	if (ToX < 0 || ToY < 0 || ToX > GetWorld()->GetWidth() || ToY > GetWorld()->GetHeight()) return 0.0;
	return sqrt(double(GetX() - ToX)*double(GetX() - ToX) + double(GetY() - ToY)*double(GetY() - ToY));
}

double Player::GetDistanceTo(Element *SomeObject) {
	return GetDistanceTo(SomeObject->GetX(), SomeObject->GetY());
}

void Player::MoveTo(int ToX, int ToY) {
	if (ToX < 0 || ToY < 0 || ToX > GetWorld()->GetWidth() || ToY > GetWorld()->GetHeight()) return;
	if (GetX() == ToX && GetY() == ToY) return;
	double AngleTo = GetAngleTo(ToX, ToY);
	Rotate(AngleTo);
	int NewX = GetX() + (double)GetSpeed()*cos(Angle),
		NewY = GetY() - (double)GetSpeed()*sin(Angle);
	if (NewX - GetR() < 0) NewX = GetR();
	if (NewY - GetR() < 0) NewY = GetR();
	if (NewX + GetR() > GetWorld()->GetWidth()) NewX = GetWorld()->GetWidth() - GetR();
	if (NewY + GetR() > GetWorld()->GetHeight()) NewY = GetWorld()->GetHeight() - GetR();

	int OldX = GetX(), OldY = GetY();
	if (NewX < 0 || NewY < 0 || NewX > GetWorld()->GetWidth() || NewY > GetWorld()->GetHeight()) return;
	SetCoords(NewX, NewY);
	auto Players = GetWorld()->GetPlayers();
	for (auto i = 0; i < Players.size(); ++i) {
		if (Players[i] == this) continue;
		if (GetDistanceTo(Players[i]) < double(GetR() + Players[i]->GetR() + 1e-2)) {
			SetCoords(OldX, OldY);
			return;
		}
	}

	auto Blocks = GetWorld()->GetBlocks();
	for (auto i = 0; i < Blocks.size(); ++i) {
		if (GetDistanceTo(Blocks[i]) < double(GetR() + Blocks[i]->GetR() + +1e-2)) {
			SetCoords(OldX, OldY);
			return;
		}
	}
}

void Player::MoveTo(Element *SomeObject) {
	MoveTo(SomeObject->GetX(), SomeObject->GetY());
}

void Player::Strike(Player *EP) {

}

void Player::StepForward()
{
	int x = 50.0 * cos(Angle);
	int y = 50.0 * sin(Angle);
	MoveTo(x, y);
}

void Player::StepBackward()
{
	int x = 50.0 * cos(Angle + M_PI);
	int y = 50.0 * sin(Angle + M_PI);
	MoveTo(x, y);
}

void Player::StepLeft()
{
	int x = 50.0 * cos(Angle - M_PI/2);
	int y = 50.0 * sin(Angle - M_PI/2);
	MoveTo(x, y);
}

void Player::StepRight()
{
	int x = 50.0 * cos(Angle + M_PI / 2);
	int y = 50.0 * sin(Angle + M_PI / 2);
	MoveTo(x, y);
}

