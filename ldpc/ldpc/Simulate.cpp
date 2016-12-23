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
	//��Profile
	ReadProfile(ldpc, channel, modulate);
	//��ʼ��ldpc
	ldpc.Initial(CodeFileName);
	ldpc.ConstructGenMatrix();//����G����
	//��ʼ��modulate
	modulate.Initial(ldpc.CodeLen);
	//��ʼ��channel
	channel.Initial(modulate.SymbolLen);
	//������Ϣ�ռ�
	MsgSeq = new int[ldpc.MsgLen];
	regPN[0] = 1;	regPN[1] = 0;	regPN[2] = 1;	regPN[3] = 0;
	regPN[4] = 0;	regPN[5] = 0;	regPN[6] = 1;	regPN[7] = 1;
	regPN[8] = 0;	regPN[9] = 0;	regPN[10] = 1;
	MsgLen = ldpc.MsgLen;
	//�����ʼ�����
	cout << "  Code File: " << CodeFileName << "\n  MaxInteration: " << ldpc.maxIteration
		<< "\tStart: " << startSNR << "  Stop: " << stopSNR << "  Step: " << stepSNR << '\n'
		<< "  Modulation Type:  " << modulate.ModualtionType << "  Decoder Algorithm:  " << ldpc.DecodeAlgorithm << '\n' << endl;
}


void CSimulate::ReadProfile(CLDPC & ldpc, CChannel & channel, CModulate & modulate)
{
	/*��Profile�ļ�*/
	ifstream fin("Profile.txt");
	if (!fin.is_open())
	{
		cerr << "Cannot open Profile.txt" << endl;
		exit(0);
	}
	string rub;
	/*��ȡ�������*/
	fin >> rub >> rub;
	fin >> rub >> rub >> startSNR;
	fin >> rub >> rub >> stepSNR;
	fin >> rub >> rub >> stopSNR;
	fin >> rub >> rub >> channel.RandomSeed;
	fin >> rub >> rub >> ldpc.maxIteration;
	/*��ȡldpc���ļ�*/
	fin >> rub >> rub;
	fin >> rub >> rub >> CodeFileName;
	/*��ȡ���Ʋ���*/
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
