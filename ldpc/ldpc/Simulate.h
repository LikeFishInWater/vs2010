#pragma once

#include "LDPC.h"
#include "Channel.h"
#include "Modulate.h"

class CSimulate
{
public:
	CSimulate();
	~CSimulate();
	void Initial(CLDPC & ldpc, CChannel & channel, CModulate & modulate);
	void ReadProfile(CLDPC & ldpc, CChannel & channel, CModulate & modulate);
	double startSNR;
	double stepSNR;
	double stopSNR;
	string CodeFileName;
	unsigned long errorFrames;
	unsigned long errorBits;
	unsigned long simulationCycles;
	double BER;
	double FER;
	double EbN0;
	int *MsgSeq;
	void GenMsgSeq();
	double sigma;
	int GenPN();
	int regPN[11];
	unsigned long MsgLen;
};

