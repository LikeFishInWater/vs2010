#include "Comm.h"


CComm::CComm(void)
	: TX_MSG(nullptr)
	, TX_CODE(nullptr)
	, TX_MOD_SYM(nullptr)
	, RX_MOD_SYM(nullptr)
	, RX_LLR(nullptr)
	, RX_DECODE(nullptr)
	, ReLdsStr(nullptr)
	, ImLdsStr(nullptr)
	, sigma_n(0)
	, CodeRate(0)
	, randomMsg(0)
{
}

CComm::~CComm(void)
{
}

bool CComm::Initial()
{
	LDPC.Initial("MN_4000.txt");
	LDPC.ConstructGenMatrix();
	CODE_LEN = LDPC.CodeLen;
	MSG_LEN = CODE_LEN/2;
	CodeRate = 1 / 2.0;
	
	ite = 7;
	randomMsg = 1;
	randomseed = 173;
	TX_MSG = new int*[16];
	for (int i = 0; i < 16; i++)
	{
		TX_MSG[i] = new int[MSG_LEN];
	}
	TX_CODE = new int*[16];
	for (int i = 0; i < 16; i++)
	{
		TX_CODE[i] = new int[CODE_LEN];
	}
	TX_MOD_SYM = new CComplex*[12];
	for (int i = 0; i < 12; i++)
	{
		TX_MOD_SYM[i] = new CComplex[CODE_LEN];
	}
	RX_MOD_SYM = new CComplex*[12];
	for (int i = 0; i < 12; i++)
	{
		RX_MOD_SYM[i] = new CComplex[CODE_LEN];
	}
	/*LDPC_RX = new CComplex*[16];
	for (int i = 0; i < 16; i++)
	{
		LDPC_RX[i] = new CComplex[CODE_LEN];
	}*/
	RX_LLR = new double*[16];
	for (int i = 0; i < 16; i++)
	{
		RX_LLR[i] = new double[CODE_LEN];
	}
	RX_DECODE = new int*[16];
	for (int i = 0; i < 16; i++)
	{
		RX_DECODE[i] = new int[CODE_LEN];
	}
	LDS_DECODE = new int*[16];
	for (int i = 0; i < 16; i++)
	{
		LDS_DECODE[i] = new int[CODE_LEN];
	}

	ReLdsStr = new double*[12];
	ImLdsStr = new double*[12];
	for (int i = 0; i<12; i++)
	{
		ReLdsStr[i] = new double[16];
		ImLdsStr[i] = new double[16];
	}
	for (int i = 0; i < 12; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			if (LdsStr[i][j] == 1)
			{
				ReLdsStr[i][j] = 1 / 2.0 * cos(3.14159 / 16.0 * j);
				ImLdsStr[i][j] = 1 / 2.0 * sin(3.14159 / 16.0 * j);
			}
			else
			{
				ReLdsStr[i][j] = 0;
				ImLdsStr[i][j] = 0;
			}
		}
	}

	regPN[0] = 1;	regPN[1] = 0;	regPN[2] = 1;	regPN[3] = 0;
	regPN[4] = 0;	regPN[5] = 0;	regPN[6] = 1;	regPN[7] = 1;
	regPN[8] = 0;	regPN[9] = 0;	regPN[10] = 1;
	return true;
}


int CComm::Transmission(void)
{
	GenerateMessage();
	Encode();
	Spread();
	Channel_AWGN();
	DeSpread();
	Decode();
	err();
	return 0;
}

int CComm::GenerateMessage(void)
{
	// generate TX_MSG_BIT & TX_MSG_SYM
	for (int i = 0; i<16; i++)        //TX_MSG_BIT size:16*msg_bit_len
	{
		for (int j = 0; j < MSG_LEN; j++)
		{
			if (randomMsg)
			{
				TX_MSG[i][j] = GenPN();
			}
			else
			{
				TX_MSG[i][j] = 0;
			}
		}
	}
	return 0;
}

int CComm::Encode(void)
{
	for (int i = 0; i < 16; i++)
	{
		LDPC.Encode(TX_MSG[i]);
		for (int j = 0; j < CODE_LEN; j++)
		{
			TX_CODE[i][j] = 2*LDPC.EncodingSeq[j]-1;
		}		
	}
	return 0;
}

int CComm::Spread(void)
{
	for (int i = 0; i < 12; i++)
	{
		for (int j = 0; j < CODE_LEN; j++)
		{
			TX_MOD_SYM[i][j].Real = TX_CODE[LdsStrH[i][0]][j] * ReLdsStr[i][LdsStrH[i][0]] + TX_CODE[LdsStrH[i][1]][j] * ReLdsStr[i][LdsStrH[i][1]] + TX_CODE[LdsStrH[i][2]][j] * ReLdsStr[i][LdsStrH[i][2]] + TX_CODE[LdsStrH[i][3]][j] * ReLdsStr[i][LdsStrH[i][3]];
			TX_MOD_SYM[i][j].Image = TX_CODE[LdsStrH[i][0]][j] * ImLdsStr[i][LdsStrH[i][0]] + TX_CODE[LdsStrH[i][1]][j] * ImLdsStr[i][LdsStrH[i][1]] + TX_CODE[LdsStrH[i][2]][j] * ImLdsStr[i][LdsStrH[i][2]] + TX_CODE[LdsStrH[i][3]][j] * ImLdsStr[i][LdsStrH[i][3]];
		}
	}
	return 1;
}

int CComm::Channel_AWGN(void)
{
	// only support bpsk now
	for (int i = 0; i < 12; i++)
	{
		for (int s = 0; s < CODE_LEN; s++)
		{
			RX_MOD_SYM[i][s].Real = TX_MOD_SYM[i][s].Real + Rand.Rand_Norm(0, sigma_n);
			RX_MOD_SYM[i][s].Image = TX_MOD_SYM[i][s].Image + Rand.Rand_Norm(0, sigma_n);
		}
	}
	//for (int i = 0; i < 16; i++)
	//{
	//	for (int s = 0; s < CODE_LEN; s++)
	//	{
	//		LDPC_RX[i][s].Real = double(TX_CODE[i][s]) + Rand.Rand_Norm(0, sigma_n);
	//		LDPC_RX[i][s].Image = Rand.Rand_Norm(0, sigma_n);
	//	}
	//}
	return 0;
}

int CComm::DeSpread(void)
{
	double he, ptem, ntem;
	int rpc, rnc;
	double D[12][16] = {};
	double Dtemp[12][16] = {};
	double rp[8];
	double rn[8];
	for (int ss = 0; ss < CODE_LEN; ss++)
	{
		for (int mm = 0; mm < 12; mm++)
			for (int nn = 0; nn < 16; nn++)
				D[mm][nn] = 0;
		for (int i = 0; i<ite; i++)
		{
			for (int m = 0; m<16; m++)
			{
				he = 0;
				for (int n = 0; n<12; n++)
					he += D[n][m];
				for (int n = 0; n < 12; n++)
				{
					if (LdsStr[n][m])
						D[n][m] = he - D[n][m];
				}	
			}
			for (int mm = 0; mm < 12; mm++)
				for (int nn = 0; nn < 16; nn++)
					Dtemp[mm][nn] = D[mm][nn];
			for (int n = 0; n<12; n++)
			{
				for (int c = 0; c<4; c++)
				{
					rpc = 0; rnc = 0;
					for (int t = 0; t < 16; t++)
					{
						if (status[c][t] == 1)
						{
							rp[rpc] = status[0][t] / 2.0 * D[n][LdsStrH[n][0]] + status[1][t] / 2.0 * D[n][LdsStrH[n][1]] + status[2][t] / 2.0 * D[n][LdsStrH[n][2]] + status[3][t] / 2.0 * D[n][LdsStrH[n][3]] - status[c][t] / 2.0 * D[n][LdsStrH[n][c]] - 1 / (2 * pow(sigma_n, 2.0))*(pow(RX_MOD_SYM[n][ss].Real - status[0][t] * ReLdsStr[n][LdsStrH[n][0]] - status[1][t] * ReLdsStr[n][LdsStrH[n][1]] - status[2][t] * ReLdsStr[n][LdsStrH[n][2]] - status[3][t] * ReLdsStr[n][LdsStrH[n][3]], 2.0) + pow(RX_MOD_SYM[n][ss].Image - status[0][t] * ImLdsStr[n][LdsStrH[n][0]] - status[1][t] * ImLdsStr[n][LdsStrH[n][1]] - status[2][t] * ImLdsStr[n][LdsStrH[n][2]] - status[3][t] * ImLdsStr[n][LdsStrH[n][3]], 2.0));
							rpc += 1;
						}
						else
						{
							rn[rnc] = status[0][t] / 2.0 * D[n][LdsStrH[n][0]] + status[1][t] / 2.0 * D[n][LdsStrH[n][1]] + status[2][t] / 2.0 * D[n][LdsStrH[n][2]] + status[3][t] / 2.0 * D[n][LdsStrH[n][3]] - status[c][t] / 2.0 * D[n][LdsStrH[n][c]] - 1 / (2 * pow(sigma_n, 2.0))*(pow(RX_MOD_SYM[n][ss].Real - status[0][t] * ReLdsStr[n][LdsStrH[n][0]] - status[1][t] * ReLdsStr[n][LdsStrH[n][1]] - status[2][t] * ReLdsStr[n][LdsStrH[n][2]] - status[3][t] * ReLdsStr[n][LdsStrH[n][3]], 2.0) + pow(RX_MOD_SYM[n][ss].Image - status[0][t] * ImLdsStr[n][LdsStrH[n][0]] - status[1][t] * ImLdsStr[n][LdsStrH[n][1]] - status[2][t] * ImLdsStr[n][LdsStrH[n][2]] - status[3][t] * ImLdsStr[n][LdsStrH[n][3]], 2.0));
							rnc += 1;
						}
					}	
					ptem = rp[0];
					ntem = rn[0];
					for (int t = 1; t < 8; t++)
					{
						ptem = (ptem > rp[t] ? ptem : rp[t]);
						ntem = (ntem > rn[t] ? ntem : rn[t]);
					}
					Dtemp[n][LdsStrH[n][c]] = ptem - ntem;
				}
			}
			for (int mm = 0; mm < 12; mm++)
				for (int nn = 0; nn < 16; nn++)
					D[mm][nn] = Dtemp[mm][nn];
		}
		for (int i = 0; i < 16; i++)
		{
			RX_LLR[i][ss] = 0;
			for (int m = 0; m < 12; m++)
			{
				RX_LLR[i][ss] += D[m][i];				
			}
			RX_LLR[i][ss] = -1*RX_LLR[i][ss];
		}
		for (int i = 0; i < 16; i++)
		{
			if (RX_LLR[i][ss] > 0)
				LDS_DECODE[i][ss] = -1;
			else
				LDS_DECODE[i][ss] = 1;
		}
	}
	return 1;
}

int CComm::Decode(void)
{
	/*double *temp;
	temp = new double[CODE_LEN];*/
	
	for (int i = 0; i < 16; i++)
	{
		/*for (int j = 0; j < CODE_LEN; j++)
		{
			temp[j] = -2 * LDPC_RX[i][j].Real / sigma_n / sigma_n;	
		}
		LDPC.Decode(temp);*/
		LDPC.Decode(RX_LLR[i]);
		for (int j = 0; j < CODE_LEN; j++)
		{
			RX_DECODE[i][j] = LDPC.DecodingSeq[j];
		}
	}
	//delete[]temp;
	return 1;
}

int CComm::err(void)
{
	errBit = 0;
	errFrm = 0;
	LDSerrBit = 0;
	bool flag;
	for (int i = 0; i < 16; i++)
	{
		flag = false;
		for (int j = 0; j < MSG_LEN; j++)
		{
			if (RX_DECODE[i][j] != TX_MSG[i][j])
			{
				flag = true;
				errBit += 1;
			}		
		}
		if (flag == true)
		{
			errFrm += 1;
		}
	}
	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < CODE_LEN; j++)
		{
			if (LDS_DECODE[i][j] != TX_CODE[i][j])
				LDSerrBit += 1;
		}
	}
	return 1;
}

int CComm::GenPN(void)
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

double CComm::SetEbN0(double EbNodb)
{
	// reset the seed
	Rand.IX = Rand.IY = Rand.IZ = randomseed;
	// reset the PN register
	regPN[0] = 1;	regPN[1] = 0;	regPN[2] = 1;	regPN[3] = 0;
	regPN[4] = 0;	regPN[5] = 0;	regPN[6] = 1;	regPN[7] = 1;
	regPN[8] = 0;	regPN[9] = 0;	regPN[10] = 1;
	// set sigma
	double EbN0 = pow(10.0, EbNodb / 10.0);
	//sigma_n = 0.00001;
	sigma_n = 1.0 / sqrt(2 * 1  * EbN0);
	return sigma_n;
}