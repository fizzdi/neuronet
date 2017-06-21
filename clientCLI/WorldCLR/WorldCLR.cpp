// ֳכאגםûי DLL-פאיכ.

#include "stdafx.h"
#include <exception>

#include "WorldCLR.h"


WorldCLR::RunWorld::RunWorld() {
	try {
		RunFromCLI();
	}
	catch (std::exception ex) {
		gcnew Exception(gcnew System::String(ex.what()));
	}
	catch (...)
	{
		gcnew Exception(gcnew System::String("unknown exception"));
	}
}

