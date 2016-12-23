#include "Simulation.h"
#include <iostream>
#include <fstream>
using namespace std;

CSimulation::CSimulation(void)
	: GFq(0)
	, nQAM(0)
	, PuntureVarDegree(0)
	, randomMsg(0)
	, decodeMethod(0)
	, ems_nm(0)
	, ems_nc(0)
	, ems_factor(0)
	, ems_offset(0)
	, maxIter(0)
	, minSimCycle(0)
	, minErrFrame(0)
	, randomseed(0)
	, snrBegin(0)
	, snrStep(0)
	, snrStop(0)
	, current_time(NULL)
	, simCycle(0)
	, errFrame(0)
	, errBit(0)
	, errSym(0)
	, FER(0)
	, BER(0)
	, SER(0)
	, EbN0(0)
	, showSimFrameStep(0)
	, tems_nr(0)
	, tems_nc(0)
	, tems_factor(0)
	, tems_offset(0)
	, U_errFrame(0)
	, U_errBit(0)
	, U_errSym(0)
	, U_FER(0)
	, U_BER(0)
	, U_SER(0)
{
}


CSimulation::~CSimulation(void)
{
}


int CSimulation::Initial(string profilename)
{
	// read profile
	string rub;
	ifstream profile(profilename);
	if (!profile.is_open())
	{
		cerr << "Cannot open " << profilename << endl;
		exit(-1);
	}
	profile >> rub >> GFq;
	profile >> rub >> rub >> NonBinaryFileName;
	profile >> rub >> rub >> PuntureVarDegree;
	
	profile >> rub >> rub >> decodeMethod;
	profile >> rub >> rub >> maxIter;
	
	// ems
	profile >> rub >> rub >> ems_nm;
	profile >> rub >> rub >> ems_nc;
	profile >> rub >> rub >> ems_factor;
	profile >> rub >> rub >> ems_offset;

	// t-ems
	profile >> rub >> rub >> tems_nr;
	profile >> rub >> rub >> tems_nc;
	profile >> rub >> rub >> tems_factor;
	profile >> rub >> rub >> tems_offset;

	// channel
	profile >> rub >> rub >> snrBegin;
	profile >> rub >> rub >> snrStep;
	profile >> rub >> rub >> snrStop;
	
	// modulation
	profile >> rub >> nQAM;
	profile >> rub >> ConstellationFileName;

	profile >> rub >> rub >> randomMsg;
	profile >> rub >> rub >> rub >> minErrFrame;
	profile >> rub >> rub >> rub >> minSimCycle;
	profile >> rub >> rub >> randomseed;
	profile >> rub >> rub >> rub >> showSimFrameStep;
	profile.close();

	// record file
	fout.open("NBLDPC.Decode.Sim.Record.txt", ios::app);

	// record time
	t = time(NULL);
	current_time = localtime(&t);

	// EbN0
	EbN0 = snrBegin - snrStep;
	return 0;
}


int CSimulation::Show(int mode)
{
	switch(mode)
	{
	case Screen_Logo:
		{
			cout << "******************************************************************************\n"
				<< "*******       NBLDPC Decoding Simulation, Powered by Xiaoshi Sha       *******\n"
				<< "******************************************************************************\n" << endl;
			break;
		}
	case Screen_Conf:
		{
			cout << "Configuration lists as follows:\n"
				<< "Code: " << NonBinaryFileName << "\nGF: " << GFq << "\tPuncture Variable Degree: " << PuntureVarDegree
				<< "\tMaximum iterations: " << maxIter << endl;
			switch(decodeMethod)
			{
			case BP_DECODE:
				{
					cout << "Algorithm: " << "BP" << endl;
					break;
				}
			case EMS_DECODE:
				{
					cout << "Algorithm: " << "EMS";
					cout << "\tNm: " << ems_nm << "\tNc: " << ems_nc << "\tFactor: " << ems_factor << "\tOffset: " << ems_offset << endl;
					break;
				}
			case MinMax_DECODE:
				{
					cout << "Algorithm: " << "Min Max" << endl;
					break;
				}
			case T_EMS_DECODE:
				{
					cout << "Algorithm: " << "Trellis EMS";
					cout << "\tNr: " << tems_nr << "\tNc: " << tems_nc << "\tFactor: " << tems_factor << "\tOffset: " << tems_offset << endl;
					break;
				}
			case T_MinMax_DECODE:
				{
					cout << "Algorithm: " << "Trellis Min Max" << endl;
					break;
				}
			}
			cout << "Modulation: " << nQAM << "-QAM" << "\tConf: " << ConstellationFileName << endl;
			cout << "SNR: " << snrBegin << ":" << snrStep << ":" << snrStop;
			switch(randomMsg)
			{
			case Rand_Seq:
				{
					cout << "\tRandom Sequence\n" 
						<< "Min Err Frame: " << minErrFrame << "\tMin Sim Frame: " << minSimCycle 
						<< "\tRand Seed: " << randomseed << endl;
					break;
				}
			case All_0_Seq:
				{
					cout << "\tALL 0 sequence\n" 
						<< "Min Err Frame: " << minErrFrame << "\tMin Sim Frame: " << minSimCycle 
						<< "\tRand Seed: " << randomseed << endl;
					break;
				}
			}
			break;
		}
	case File_Conf:
		{
			
			fout << "\n\n\nConfiguration lists as follows:\n"
				<< "Code: " << NonBinaryFileName << "\nGF: " << GFq << "\tPuncture Variable Degree: " << PuntureVarDegree
				<< "\tMaximum iterations: " << maxIter << endl;
			switch(decodeMethod)
			{
			case BP_DECODE:
				{
					fout << "Algorithm: " << "BP" << endl;
					break;
				}
			case EMS_DECODE:
				{
					fout << "Algorithm: " << "EMS";
					fout << "\tNm: " << ems_nm << "\tNc: " << ems_nc << "\tFactor: " << ems_factor << "\tOffset: " << ems_offset << endl;
					break;
				}
			case MinMax_DECODE:
				{
					fout << "Algorithm: " << "Min Max" << endl;
					break;
				}
			case T_EMS_DECODE:
				{
					fout << "Algorithm: " << "Trellis EMS";
					fout << "\tNr: " << tems_nr << "\tNc: " << tems_nc << "\tFactor: " << tems_factor << "\tOffset: " << tems_offset << endl;
					break;
				}
			case T_MinMax_DECODE:
				{
					fout << "Algorithm: " << "Trellis Min Max" << endl;
					break;
				}
			}
			fout << "Modulation: " << nQAM << "-QAM" << "\tConf: " << ConstellationFileName << endl;
			fout << "SNR: " << snrBegin << ":" << snrStep << ":" << snrStop;
			switch(randomMsg)
			{
			case Rand_Seq:
				{
					fout << "\tRandom Sequence\n" 
						<< "Min Err Frame: " << minErrFrame << "\tMin Sim Frame: " << minSimCycle 
						<< "\tRand Seed: " << randomseed  << endl;
					break;
				}
			case All_0_Seq:
				{
					fout << "\tALL 0 sequence\n" 
						<< "Min Err Frame: " << minErrFrame << "\tMin Sim Frame: " << minSimCycle 
						<< "\tRand Seed: " << randomseed << endl;
					break;
				}
			}
			
			break;
		}
	case Screen_Head:
		{
			cout << "EbN0\t" << "Error\t"  << "BER\t\t" << "SER\t\t" << "FER\t\t" << "TIME" << endl;
//			cout << "EbN0\t" << "BER\t\t"  << "FER\t\t" << "U-BER\t\t"  << "U-FER"  << endl;
			break;
		}
	case File_Head:
		{
			fout << "EbN0\t" << "Error\t"  
				<< "BER\t\t" << "SER\t\t" << "FER\t\t" 
				<< "U-BER\t\t" << "U-SER\t\t" << "U-FER\t\t" 
				<< "SIM\t\t" << "TIME" << endl;
			break;
		}
	case File_Begin_Time:
		{
			fout << "\nBegin in: " << current_time -> tm_year + 1900 << "-" << current_time -> tm_mon + 1 
				<< "-" << current_time -> tm_mday << "," << current_time -> tm_hour << ":" 
				<< current_time -> tm_min << ":" << current_time -> tm_sec << endl;
			break;
		}
	case File_End_Time:
		{
			fout << "End in: " << current_time -> tm_year + 1900 << "-" << current_time -> tm_mon + 1 
				<< "-" << current_time -> tm_mday << "," << current_time -> tm_hour << ":" 
				<< current_time -> tm_min << ":" << current_time -> tm_sec << "\n\n" << endl;
			break;
		}
	case Screen_Sim_Data:
		{
			if(long(simCycle) % showSimFrameStep == 0)
			{
				cout << defaultfloat << EbN0 << '\t' << errFrame<< '\t'; 
				cout << scientific << BER << '\t' << SER << '\t' 
					<< FER  << '\t' << simCycle << '\r';
//				cout << defaultfloat << EbN0 << '\t'; 
//				cout << scientific << BER << '\t' << FER << '\t' << U_BER << '\t' << U_FER << '\r';
			}
			break;
		}
	case Screen_Sim_End_Data:
		{
			stop = clock();
			cout << defaultfloat << EbN0 << '\t' << errFrame << '\t'; 
			cout << scientific << BER << '\t'  << SER << '\t' 
				<< FER << '\t' << 1.0 * (stop - start) / CLOCKS_PER_SEC << endl;
//			cout << defaultfloat << EbN0 << '\t'; 
//			cout << scientific << BER << '\t' << FER << '\t'
//					<< U_BER << '\t' << U_FER << endl;
			break;
		}
	case File_Sim_Data:
		{
			/*
			if(long(simCycle) % (showSimFrameStep * 100000) == 0)
			{
				status = clock();
				if(int(1.0 * (status - start) / CLOCKS_PER_SEC / 450) % 2 == 0)
				{
					fout << defaultfloat << EbN0 << '\t' << errFrame<< '\t'; 
					fout << scientific << BER << '\t' << SER << '\t' 
						<< FER  << '\t' << simCycle 
						<< '\t' << 1.0 * (status - start) / CLOCKS_PER_SEC << endl;
				}
			}
			*/
			break;
		}
	case File_Sim_End_Data:
		{
			stop = clock();
			fout << defaultfloat << EbN0 << '\t' << errFrame << '\t'; 
			fout << scientific << BER << '\t' << SER << '\t' << FER << '\t'
				<< U_BER << '\t' << U_SER << '\t' << U_FER << '\t'
				<< simCycle << '\t' << 1.0 * (stop - start) / CLOCKS_PER_SEC << " **** "<< endl;
				
			break;
		}
	}
	return 0;
}


int CSimulation::ClearSimuCount(void)
{
	simCycle = 0;
	errFrame = errBit = errSym = 0;
	U_errFrame = U_errBit = U_errSym = 0;
	start = clock();
	return 0;
}


bool CSimulation::NextSNR(void)
{
//	EbN0 = Simulation.snrBegin; EbN0 <= Simulation.snrStop; EbN0 += Simulation.snrStep
	EbN0 += snrStep;
	if(EbN0 > snrStop)
		return false;
	else
		return true;
}


bool CSimulation::SimulateThisSNR(void)
{
	return (simCycle ++ <= minSimCycle) && (errFrame <= minErrFrame);
//	return (simCycle ++ <= minSimCycle) || (U_errFrame <= minErrFrame);
}