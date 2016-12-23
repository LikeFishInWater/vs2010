#pragma once

#include <string>
using namespace std;

class CNode
{
public:
	CNode();
	~CNode();
	unsigned int Degree;
	//	int maxDegree;
	unsigned long *Index;
	unsigned long *EdgeNo;
};

class CLDPC
{
public:
	CLDPC();
	~CLDPC();
	unsigned long CodeLen;
	unsigned long MsgLen;
	unsigned long ChkLen;
	double Rate;
	int maxIteration;
	CNode *checkNode;
	CNode *varaiableNode;
	unsigned int maxVarDegree;
	unsigned int maxChkDegree;
	double *Q;
	double *R;
	unsigned long EdgeSum;
	int *SourceSeq;
	int *EncodingSeq;
	int *DecodingSeq;
	void Encode(int * MsgSeq);
	unsigned long CalErrorBits();
	//	double *L;
	void Initial(string CodeFileName);
	void ConstructGenMatrix();
	long *ExchangedinEncoding;
	CNode *encodingNode;
	bool EncodeCheck();
	int DecodeAlgorithm;
	double MS_alpha;
	void Decode(double * LLR);
	void Decode_BP(double * LLR);
	void Decode_MS(double * LLR, double alpha);
	void Decode_BF(double * LLR);
	unsigned long *BF_varSeq;
	unsigned long *BF_chkSeq;
};