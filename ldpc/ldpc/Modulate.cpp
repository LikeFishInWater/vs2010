
#include "Modulate.h"
#include <math.h>

CModulate::CModulate()
: ModualtionType(0)
, SourceLen(0)
, SymbolLen(0)
, ModSeq(nullptr)
, DemodSeq(nullptr)
{
}


CModulate::~CModulate()
{
	delete[] ModSeq;
	delete[] DemodSeq;
}


void CModulate::Initial(unsigned long codelen)
{
	SourceLen = codelen;
	DemodSeq = new double[SourceLen];
	switch (ModualtionType)
	{
	case 0://BPSK
	{
			   SymbolLen = SourceLen;
			   ModSeq = new CComplex[SymbolLen];
			   break;
	}
	case 1://QPSK
	{
			   SymbolLen = SourceLen / 2;
			   ModSeq = new CComplex[SymbolLen];
			   break;
	}
	default://not define
	{
			SymbolLen = SourceLen;
			ModSeq = nullptr;
			break;
	}
	}
}


void CModulate::Modulation(int * SourceSeq)
{
	switch (ModualtionType)
	{
	case 0://BPSK
	{
			   for (unsigned long i = 0; i < SymbolLen; i++)
			   {
				   ModSeq[i].Real = 1 - 2 * double(SourceSeq[i]);
				   ModSeq[i].Image = 0;
			   }
			   break;
	}
	case 1://QPSK
	{
			   for (unsigned long i = 0; i < SymbolLen; i++)
			   {
				   if (SourceSeq[2 * i] == 0 && SourceSeq[2 * i + 1] == 0)
				   {
					   ModSeq[i].Real = 1;
					   ModSeq[i].Image = 0;
				   }
				   else if (SourceSeq[2 * i] == 1 && SourceSeq[2 * i + 1] == 0)
				   {
					   ModSeq[i].Real = 0;
					   ModSeq[i].Image = 1;
				   }
				   else if (SourceSeq[2 * i] == 0 && SourceSeq[2 * i + 1] == 1)
				   {
					   ModSeq[i].Real = 0;
					   ModSeq[i].Image = -1;
				   }
				   else if (SourceSeq[2 * i] == 1 && SourceSeq[2 * i + 1] == 1)
				   {
					   ModSeq[i].Real = -1;
					   ModSeq[i].Image = 0;
				   }
			   }
			   break;
	}
	default:
		break;
	}
}


void CModulate::Demodulation(CComplex * ReciveSeq, double sigma)
{	
	switch (ModualtionType)
	{
	case 0://BPSK
	{
			  for (unsigned long i = 0; i < SymbolLen; i++)
				  DemodSeq[i] = 2 * ReciveSeq[i].Real / (sigma * sigma);
			  break;
	}
	case 1://QPSK
	{
			   double p00, p01, p10, p11;
			   for (unsigned long i = 0; i < SymbolLen; i++)
			   {
				   p00 = exp(0 - ((ReciveSeq[i].Real - 1)* (ReciveSeq[i].Real - 1) + (ReciveSeq[i].Image - 0) * (ReciveSeq[i].Image - 0)) / (2 * sigma * sigma));
				   p10 = exp(0 - ((ReciveSeq[i].Real - 0)* (ReciveSeq[i].Real - 0) + (ReciveSeq[i].Image - 1) * (ReciveSeq[i].Image - 1)) / (2 * sigma * sigma));
				   p01 = exp(0 - ((ReciveSeq[i].Real - 0)* (ReciveSeq[i].Real - 0) + (ReciveSeq[i].Image + 1) * (ReciveSeq[i].Image + 1)) / (2 * sigma * sigma));
				   p11 = exp(0 - ((ReciveSeq[i].Real + 1)* (ReciveSeq[i].Real + 1) + (ReciveSeq[i].Image - 0) * (ReciveSeq[i].Image - 0)) / (2 * sigma * sigma));
				   DemodSeq[2 * i] = log((p00 + p01) / (p10 + p11));
				   DemodSeq[2 * i + 1] = log((p00 + p10) / (p01 + p11));
			   }
			   break;
	}
	default:
		break;
	}
}