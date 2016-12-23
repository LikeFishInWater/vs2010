#include <iostream>
#include <fstream>
#include <string>
#include "Simulate.h"
using namespace std;

CSimulate::CSimulate()
: startSNR(0)
, stepSNR(0)
, stopSNR(0)
, errorFrames(0)
, errorBits(0)
, simulationCycles(0)
, BER(0)
, FER(0)
, EbN0(0)
, MsgSeq(nullptr)
, sigma(0)
, MsgLen(0)
{
}


CSimulate::~CSimulate()
{
	delete[] MsgSeq;
}


void CSimulate::Initial(CLDPC & ldpc, CChannel & channel, CModulate & modulate)
{
	//读Profile
	ReadProfile(ldpc, channel, modulate);
	//初始化ldpc
	ldpc.Initial(CodeFileName);
	ldpc.ConstructGenMatrix();//生成G矩阵
	//初始化modulate
	modulate.Initial(ldpc.CodeLen);
	//初始化channel
	channel.Initial(modulate.SymbolLen);
	//分配消息空间
	MsgSeq = new int[ldpc.MsgLen];
	regPN[0] = 1;	regPN[1] = 0;	regPN[2] = 1;	regPN[3] = 0;
	regPN[4] = 0;	regPN[5] = 0;	regPN[6] = 1;	regPN[7] = 1;
	regPN[8] = 0;	regPN[9] = 0;	regPN[10] = 1;
	MsgLen = ldpc.MsgLen;
	//输出初始化结果
	cout << "  Code File: " << CodeFileName << "\n  MaxInteration: " << ldpc.maxIteration
		<< "\tStart: " << startSNR << "  Stop: " << stopSNR << "  Step: " << stepSNR << '\n'
		<< "  Modulation Type:  " << modulate.ModualtionType << "  Decoder Algorithm:  " << ldpc.DecodeAlgorithm << '\n' << endl;
}


void CSimulate::ReadProfile(CLDPC & ldpc, CChannel & channel, CModulate & modulate)
{
	/*打开Profile文件*/
	ifstream fin("Profile.txt");
	if (!fin.is_open())
	{
		cerr << "Cannot open Profile.txt" << endl;
		exit(0);
	}
	string rub;
	/*读取仿真参数*/
	fin >> rub >> rub;
	fin >> rub >> rub >> startSNR;
	fin >> rub >> rub >> stepSNR;
	fin >> rub >> rub >> stopSNR;
	fin >> rub >> rub >> channel.RandomSeed;
	fin >> rub >> rub >> ldpc.maxIteration;
	/*读取ldpc码文件*/
	fin >> rub >> rub;
	fin >> rub >> rub >> CodeFileName;
	/*读取调制参数*/
	fin >> rub >> rub;
	fin >> rub >> modulate.ModualtionType;
	fin >> rub >> rub;
	fin >> rub >> ldpc.DecodeAlgorithm;
	fin >> rub >> rub >> ldpc.MS_alpha;
	fin.close();
}


void CSimulate::GenMsgSeq()
{
	for (unsigned long i = 0; i < MsgLen; i++)
		MsgSeq[i] = GenPN();
}


int CSimulate::GenPN()
{
	/*shift the shift register*/
	for (int i = 10; i >= 1; i--)
	{
		regPN[i] = regPN[i - 1];
	}
	/*calculate the output*/
	regPN[0] = regPN[10] ^ regPN[3];
	/*output the result*/
	return regPN[10];
}
