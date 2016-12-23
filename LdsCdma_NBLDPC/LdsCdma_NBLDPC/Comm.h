#pragma once
#include "Simulation.h"
#include "Complex.h"
#include "GF.h"
#include "Rand.h"
#include <string>
#include <iostream>
#include <fstream>
#include "nbldpc.h"
#include <cmath>
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
static int LdsStrH[12][4] = { { 2,7,9,12 },{ 3,6,12,14 },{ 1,4,6,9 },{ 3,5,10,13 } ,{ 2,8,14,15 } ,{ 0,5,9,14 },{ 3,7,8,11 } ,{ 1,5,11,15 } ,{ 0,6,13,15 },{ 2,4,11,13 } ,{ 1,8,10,12 },{ 0,4,7,10 } };
static int status[4][16] = { { -1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1 },{ -1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1 },{ -1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1 },{ -1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1 } };
class CComm
{
public:
	CComm(void);
	~CComm(void);
	
	double** ReLdsStr;
	double** ImLdsStr;
	

	// tx-msg
	int** TX_MSG_SYM;
	int** TX_MSG_BIT;
	// tx-code
	int** TX_CODE_SYM;
	int** TX_CODE_BIT;
	CComplex ** TX_MOD_SYM;
	CComplex ** RX_MOD_SYM;
	double ** RX_LLR_BIT;
	double *** RX_LLR_SYM;
	int** RX_DECODE_SYM;
	int** RX_DECODE_BIT;

	int MSG_BIT_LEN;
	int MSG_SYM_LEN;
	int CODE_BIT_PER_SYM;
	int CODE_BIT_LEN;
	int CODE_SYM_LEN;
	int codeOrder;


public:
	CNBLDPC NBLDPC;
	CRand Rand;
	//	bool Initial(int GFq, int nQAM, string NBLDPCFileName, int PunctureVarDegree, int maxIter, int randMsg, int deocdemethod, int ems_nm, int ems_nc, double ems_factor, double ems_offset, int randseed);
	bool Initial(CSimulation &sim);
	int Transmission(void);
	int GenerateMessage(void);
	int Encode(void);
	int Modulate(void);
	int Channel_AWGN(void);
	int Demodulate(void);
	int Decode(void);
	int regPN[11];
	int GenPN(void);
	// noise
	double sigma_n;
	double CodeRate;
	double SetEbN0(CSimulation &sim);
	// compute the error frame, bit, symbol
	//	int Err(double simcycle, double& errFrame, double& errBit, double& errSym, double& fer, double& ber, double& ser);
	int Err(CSimulation &sim);
	int randomMsg;

	double **RX_DEMOD_LLR_SYM;
	bool DecodeCorrect;


};

