#include "Complex.h"

#pragma once

struct RandSeed
{
	unsigned long IX;
	unsigned long IY;
	unsigned long IZ;
};

class CChannel
{
public:
	CChannel();
	~CChannel();
	CComplex *SymbolSeq;
	unsigned long SymbolLen;
//	double startSNR;
//	double stepSNR;
//	double stopSNR;
//	double Sigma;
//	int ChannelType;
	RandSeed RS;
	int RandomSeed;


	void Initial(unsigned long len);
//	double *Real;
//	double *Image;
	double Random_Uniform(RandSeed & rs);
	double Random_Norm(double sigma, RandSeed & rs);
	void AWGNChannel(CComplex * inSymbolSeq, double sigma);
};



