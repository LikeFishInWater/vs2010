#include "Channel.h"
#include <math.h>
#include <fstream>
using namespace std;

CChannel::CChannel()
:RandomSeed(0)
, SymbolSeq(nullptr)
{
}


CChannel::~CChannel()
{
}


double CChannel::Random_Uniform(RandSeed & rs)
{
	double temp = 0.0;
	rs.IX = (rs.IX * 249) % 61967;
	rs.IY = (rs.IY * 251) % 63443;
	rs.IZ = (rs.IZ * 252) % 63599;
	temp = (((double)rs.IX) / ((double)61967)) + (((double)rs.IY) / ((double)63443))
		+ (((double)rs.IZ) / ((double)63599));
	temp -= (int)temp;
	return temp;
}


double CChannel::Random_Norm(double sigma, RandSeed & rs)
{
	double u1, u2, u;
	u1 = Random_Uniform(rs);
	u2 = Random_Uniform(rs);
	u = sigma* cos(2 * 3.1415926535897932384626433832795 * u2) * sqrt(-2.0 * log(1.0 - u1));
	return u;
}


void CChannel::Initial(unsigned long len)
{
	RS.IX = RS.IY = RS.IZ = RandomSeed;
	SymbolLen = len;
	SymbolSeq = new CComplex[SymbolLen];
}


void CChannel::AWGNChannel(CComplex * inSymbolSeq, double sigma)
{
//	ofstream fout("modsynbol.txt");
	for (unsigned long i = 0; i < SymbolLen; i++)
	{
//		fout << inSymbolSeq[i].Real << endl;
		SymbolSeq[i].Real = Random_Norm(sigma, RS) + inSymbolSeq[i].Real;
		SymbolSeq[i].Image = Random_Norm(sigma, RS) + inSymbolSeq[i].Image;
	}
//	fout.close();
}
