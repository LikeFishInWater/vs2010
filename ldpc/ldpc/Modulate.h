#include "Complex.h"

#pragma once

class CModulate
{
public:
	CModulate();
	~CModulate();
	int ModualtionType;
	unsigned long SourceLen;
	unsigned long SymbolLen;
	void Initial(unsigned long codelen);
//	double *Real;
//	double *Image;
	double *DemodSeq;
	CComplex * ModSeq;
	void Modulation(int * SourceSeq);
	void Demodulation(CComplex * ReciveSeq, double sigma);
};

