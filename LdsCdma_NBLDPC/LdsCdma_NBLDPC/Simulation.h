#pragma once

#define BP_DECODE 1
#define EMS_DECODE 2
#define MinMax_DECODE 3
#define T_EMS_DECODE 4
#define T_MinMax_DECODE 5

#define Rand_Seq 1
#define All_0_Seq 0

#define Screen_Logo 0
#define Screen_Conf	1
#define File_Conf	2
#define Screen_Head	3
#define File_Head	4
#define File_Begin_Time	5
#define File_End_Time	6
#define Screen_Sim_Data	7
#define Screen_Sim_End_Data	8
#define File_Sim_Data	9
#define File_Sim_End_Data	10

#include <string>
#include <fstream>
#include <ctime>
using namespace std;

class CSimulation
{
public:
	
	string ProfileFileName;
	int GFq;
	int nQAM;
	int PuntureVarDegree;
	int randomMsg;
	int decodeMethod;
	int ems_nm;
	int ems_nc;
	double ems_factor;
	double ems_offset;
	int maxIter;
	int minSimCycle;
	int minErrFrame;
	int randomseed;
	double snrBegin;
	double snrStep;
	double snrStop;
	string NonBinaryFileName;
	ofstream fout;
	time_t t;
	struct tm* current_time;
	double simCycle;
	double errFrame;
	double errBit;
	double errSym;
	double FER;
	double BER;
	double SER;
	clock_t start;
	clock_t stop;
	clock_t status;
	int showSimFrameStep;
	int tems_nr;
	int tems_nc;
	double tems_factor;
	double tems_offset;
	string ConstellationFileName;
	double EbN0;
	CSimulation(void);
	~CSimulation(void);
	int Initial(string profilename);
	int Show(int mode);
	int ClearSimuCount(void);
	bool NextSNR(void);
	bool SimulateThisSNR(void);
	
};

