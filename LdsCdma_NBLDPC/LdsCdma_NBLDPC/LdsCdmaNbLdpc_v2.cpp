#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>
#include "Simulation.h"
#include "Comm.h"

using namespace std;


int main()
{
	// read profile.txt
	CSimulation Simulation;
	Simulation.Initial("NBLDPC.Profile.txt");
	CComm Comm;
	Comm.Initial(Simulation);

	// show
	Simulation.Show(Screen_Head);
	Simulation.Show(File_Begin_Time);
	Simulation.Show(File_Head);

	// start simulation
	while (Simulation.NextSNR())
	{
		Simulation.ClearSimuCount();
		Comm.SetEbN0(Simulation);
		while (Simulation.SimulateThisSNR())
		{
			// transmmit one frame over the channel
			Comm.Transmission();
			// calculate the err
			Comm.Err(Simulation);
			// output to screen and file
			Simulation.Show(Screen_Sim_Data);
			//			Simulation.Show(File_Sim_Data);
		}
		Simulation.Show(Screen_Sim_End_Data);
		Simulation.Show(File_Sim_End_Data);
	}
	// simulation end
	Simulation.Show(File_End_Time);
	Simulation.fout.close();

	// pause
	system("pause");
}