#pragma once
#include "Complex.h"
#include "Rand.h"
#include <string>
#include <iostream>
#include <fstream>
#include "LDPC.h"
#include <cmath>
//#include "Simulation.h"
using namespace std;
static int LdsStr[12][16] =
	{
		{ 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0 },
		{ 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0 },
		{ 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 },
		{ 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1 },
		{ 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 },
		{ 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0 },
		{ 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1 },
		{ 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1 },
		{ 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0 },
		{ 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0 },
		{ 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0 }
	};
static int LdsStrV[16][3] = { { 5,8,11 },{ 2,7,10 },{ 0,4,9 },{ 1,3,6 },{ 2,9,11 } ,{ 3,5,7 },{ 1,2,8 },{ 0,6,11 } ,{ 4,6,10 } ,{ 0,2,5 } ,{ 3,10,11 },{ 6,7,9 },{ 0,1,10 },{ 3,8,9 } ,{ 1,4,5 } ,{ 4,7,8 } };
static const int LdsStrH[12][4] = { { 2,7,9,12 },{ 3,6,12,14 },{ 1,4,6,9 },{ 3,5,10,13 } ,{ 2,8,14,15 } ,{ 0,5,9,14 },{ 3,7,8,11 } ,{ 1,5,11,15 } ,{ 0,6,13,15 },{ 2,4,11,13 } ,{ 1,8,10,12 },{ 0,4,7,10 } };
static int status[4][16] = { { -1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1 },{ -1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1 },{ -1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1 },{ -1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1 } };
class CComm
{
public:
	CComm(void);
	~CComm(void);
	
	double** ReLdsStr;
	double** ImLdsStr;
	

	// tx-msg
	int ** TX_MSG;
	// tx-code
	int** TX_CODE;
	CComplex ** TX_MOD_SYM;
	CComplex ** RX_MOD_SYM;
	//CComplex ** LDPC_RX;
	double ** RX_LLR;
	int ** LDS_DECODE;
	int** RX_DECODE;


	int MSG_LEN;
	int CODE_LEN;
	int ite;

	int randomMsg; 
	int randomseed;
	double sigma_n;
	double CodeRate;

	int errBit, errFrm;
	int LDSerrBit, LDSerrFrm;
public:
	CLDPC LDPC;
	CRand Rand;
	int regPN[11];

	//	bool Initial(int GFq, int nQAM, string NBLDPCFileName, int PunctureVarDegree, int maxIter, int randMsg, int deocdemethod, int ems_nm, int ems_nc, double ems_factor, double ems_offset, int randseed);
	bool Initial();
	int Transmission(void);
	int GenerateMessage(void);
	int Encode(void);
	int Spread(void);
	int Channel_AWGN(void);
	int DeSpread(void);
	int Decode(void);
	int GenPN(void);
	double SetEbN0(double EbNodb);
	int err(void);
};

