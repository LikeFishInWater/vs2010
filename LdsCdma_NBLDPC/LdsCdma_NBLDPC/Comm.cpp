#include "Comm.h"


CComm::CComm(void)
	: TX_MSG_SYM(nullptr)
	, TX_MSG_BIT(nullptr)
	, TX_CODE_SYM(nullptr)
	, TX_CODE_BIT(nullptr)
	, TX_MOD_SYM(nullptr)
	, RX_MOD_SYM(nullptr)
	, RX_LLR_SYM(nullptr)
	, RX_LLR_BIT(nullptr)
	, RX_DECODE_SYM(nullptr)
	, RX_DECODE_BIT(nullptr)
	, ReLdsStr(nullptr)
	, ImLdsStr(nullptr)
	, sigma_n(0)
	, CodeRate(0)
	, randomMsg(0)
	, RX_DEMOD_LLR_SYM(NULL)
	, DecodeCorrect(false)
{
}


CComm::~CComm(void)
{
}


//bool CComm::Initial(int GFq, int nQAM, string NBLDPCFileName, int PunctureVarDegree, int maxIter, int randMsg, int deocdemethod, int ems_nm, int ems_nc, double ems_factor, double ems_offset, int randseed)
bool CComm::Initial(CSimulation &sim)
{
	NBLDPC.Initial(sim);
	MSG_SYM_LEN= NBLDPC.CodeLen/2;
	CODE_BIT_PER_SYM = 4;
	MSG_BIT_LEN = MSG_SYM_LEN * CODE_BIT_PER_SYM;
	CODE_BIT_LEN= MSG_BIT_LEN * 2;
	CODE_SYM_LEN= MSG_SYM_LEN * 2;
	codeOrder=16;
	// transfer parameters
	randomMsg = sim.randomMsg;

	// intial NBLDPC
	
	// allocate memory
	TX_MSG_BIT = new int*[16];
	for (int i = 0; i < 16; i++)
	{
		TX_MSG_BIT[i] = new int[MSG_BIT_LEN];
	}
	TX_MSG_SYM = new int*[16];
	for (int i = 0; i < 16; i++)
	{
		TX_MSG_SYM[i] = new int[MSG_SYM_LEN];
	}
	TX_CODE_SYM = new int*[16];
	for (int i = 0; i < 16; i++)
	{
		TX_CODE_SYM[i] = new int[CODE_SYM_LEN];
	}
	TX_CODE_BIT = new int*[16];
	for (int i = 0; i < 16; i++)
	{
		TX_CODE_BIT[i] = new int[CODE_BIT_LEN];
	}
	TX_MOD_SYM = new CComplex*[12];
	for (int i = 0; i < 12; i++)
	{
		TX_MOD_SYM[i] = new CComplex[CODE_BIT_LEN];
	}
	RX_MOD_SYM = new CComplex*[12];
	for (int i = 0; i < 12; i++)
	{
		RX_MOD_SYM[i] = new CComplex[CODE_BIT_LEN];
	}
	RX_LLR_BIT = new double*[16];
	for (int i = 0; i < 16; i++)
	{
		RX_LLR_BIT[i] = new double[CODE_BIT_LEN];
	}
	RX_LLR_SYM = new double**[16];
	for (int i = 0; i < 16; i++)
	{
		RX_LLR_SYM[i] = new double*[CODE_SYM_LEN];
	}
	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < CODE_SYM_LEN; j++)
		{
			RX_LLR_SYM[i][j] = new double[codeOrder-1];
		}
	}
	RX_DECODE_SYM = new int*[16];
	for (int i = 0; i < 16; i++)
	{
		RX_DECODE_SYM[i] = new int[CODE_SYM_LEN];
	}
	RX_DECODE_BIT = new int*[16];
	for (int i = 0; i < 16; i++)
	{
		RX_DECODE_BIT[i] = new int[CODE_BIT_LEN];
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
				ReLdsStr[i][j] = 1 / 2.0 * cos(3.14159 / 16 * j) ;
				ImLdsStr[i][j] = 1 / 2.0 * sin(3.14159 / 16 * j) ;
			}
			else
			{
				ReLdsStr[i][j] = 0;
				ImLdsStr[i][j] = 0;
			}
		}
	}
	// intial PN register
	regPN[0] = 1;	regPN[1] = 0;	regPN[2] = 1;	regPN[3] = 0;
	regPN[4] = 0;	regPN[5] = 0;	regPN[6] = 1;	regPN[7] = 1;
	regPN[8] = 0;	regPN[9] = 0;	regPN[10] = 1;
	// intial random generator
	Rand.IX = Rand.IY = Rand.IZ = sim.randomseed;
	return false;
}


int CComm::Transmission(void)
{
	GenerateMessage();
	Encode();
	Modulate();
	Channel_AWGN();
	Demodulate();
	Decode();
	return 0;
}


int CComm::GenerateMessage(void)
{
	// generate TX_MSG_BIT & TX_MSG_SYM
	for (int i = 0; i<16; i++)        //TX_MSG_BIT size:16*msg_bit_len
	{
		for (int j = 0; j < MSG_BIT_LEN; j++)
		{
			if (randomMsg)
			{
				TX_MSG_BIT[i][j] = GenPN();
			}
			else
			{
				TX_MSG_BIT[i][j] = 0;
			}
		}
	}
	for (int i = 0; i < 16; i++)
	{
		for (int s = 0; s < MSG_SYM_LEN; s++)
		{
			TX_MSG_SYM[i][s] = 0;
			for (int b_p_s = 0; b_p_s < CODE_BIT_PER_SYM; b_p_s++)
			{
				TX_MSG_SYM[i][s] += (TX_MSG_BIT[i][s * CODE_BIT_PER_SYM + b_p_s] << (CODE_BIT_PER_SYM - 1 - b_p_s));
			}
		}
	}
	
	
	return 0;
}


int CComm::Encode(void)
{
	// TX_MSG_SYM --> TX_CODE_SYM
	for (int i = 0; i < 16; i++)
	{
		NBLDPC.Encode(TX_MSG_SYM[i], TX_CODE_SYM[i]);
	}
	// TX_CODE_SYM --> TX_CODE_BIT
	for (int i = 0; i < 16; i++)
	{
		for (int s = 0; s < CODE_SYM_LEN; s++)
		{
			for (int b_p_s = 0; b_p_s < CODE_BIT_PER_SYM; b_p_s++)
			{
				TX_CODE_BIT[i][s * CODE_BIT_PER_SYM + b_p_s]
					= ((TX_CODE_SYM[i][s] & (1 << (CODE_BIT_PER_SYM - 1 - b_p_s))) == 0) ? 0 : 1;
				TX_CODE_BIT[i][s * CODE_BIT_PER_SYM + b_p_s] = TX_CODE_BIT[i][s * CODE_BIT_PER_SYM + b_p_s] * 2 - 1;
			}
		}
	}
	return 0;
}

int CComm::Modulate(void)
{
	//
	for (int i = 0; i < 12; i++)
	{
		for (int j = 0; j < CODE_BIT_LEN; j++)
		{
			TX_MOD_SYM[i][j].Real = TX_CODE_BIT[LdsStrH[i][0]][j] * ReLdsStr[i][LdsStrH[i][0]] + TX_CODE_BIT[LdsStrH[i][1]][j] * ReLdsStr[i][LdsStrH[i][1]] + TX_CODE_BIT[LdsStrH[i][2]][j] * ReLdsStr[i][LdsStrH[i][2]] + TX_CODE_BIT[LdsStrH[i][3]][j] * ReLdsStr[i][LdsStrH[i][3]];
			TX_MOD_SYM[i][j].Image = TX_CODE_BIT[LdsStrH[i][0]][j] * ImLdsStr[i][LdsStrH[i][0]] + TX_CODE_BIT[LdsStrH[i][1]][j] * ImLdsStr[i][LdsStrH[i][1]] + TX_CODE_BIT[LdsStrH[i][2]][j] * ImLdsStr[i][LdsStrH[i][2]] + TX_CODE_BIT[LdsStrH[i][3]][j] * ImLdsStr[i][LdsStrH[i][3]];
		}
	}
	return 0;
}


int CComm::Channel_AWGN(void)
{
	// only support bpsk now
	for (int i = 0; i < 12; i++)
	{
		for (int s = 0; s < CODE_BIT_LEN; s++)
		{
			RX_MOD_SYM[i][s].Real = TX_MOD_SYM[i][s].Real + Rand.Rand_Norm(0, sigma_n);
			RX_MOD_SYM[i][s].Image = TX_MOD_SYM[i][s].Image + Rand.Rand_Norm(0, sigma_n);
		}
	}
	return 0;
}


int CComm::Demodulate(void)
{
	double he,ptem,ntem;
	int rpc, rnc;
	double D[12][16] = {};
	double Dtemp[12][16] = {};
	double rp[8];
	double rn[8];
	for (int ss = 0; ss < CODE_BIT_LEN; ss++)
	{
		for (int mm = 0; mm < 12; mm++)
			for (int nn = 0; nn < 16; nn++)
				D[mm][nn] = 0;
		for (int i = 0; i<7; i++)
		{
			for (int m = 0; m<16; m++)
			{
				he = 0;
				for (int n = 0; n<12; n++)
					he += D[n][m];
				for (int n = 0; n<12; n++)
					if (LdsStr[n][m])
						D[n][m] = he - D[n][m];
			}
			for (int mm = 0; mm < 12; mm++)
				for (int nn = 0; nn < 16; nn++)
					Dtemp[mm][nn] = D[mm][nn];
			for (int n = 0; n<12; n++)
			{
				for (int c = 0; c<4; c++)
				{
					rpc = 0; rnc = 0;
					for (int t = 0; t<16; t++)
						if (status[c][t] == 1)
						{
							rp[rpc] = status[0][t] / 2.0 * D[n][LdsStrH[n][0]] + status[1][t] / 2.0 * D[n][LdsStrH[n][1]] + status[2][t] / 2.0 * D[n][LdsStrH[n][2]] + status[3][t] / 2.0 * D[n][LdsStrH[n][3]] - status[c][t] / 2.0 * D[n][LdsStrH[n][c]] - 1 / (2 * pow(sigma_n, 2))*(pow(RX_MOD_SYM[n][ss].Real - status[0][t] * ReLdsStr[n][LdsStrH[n][0]] - status[1][t] * ReLdsStr[n][LdsStrH[n][1]] - status[2][t] * ReLdsStr[n][LdsStrH[n][2]] - status[3][t] * ReLdsStr[n][LdsStrH[n][3]], 2) + pow(RX_MOD_SYM[n][ss].Image - status[0][t] * ImLdsStr[n][LdsStrH[n][0]] - status[1][t] * ImLdsStr[n][LdsStrH[n][1]] - status[2][t] * ImLdsStr[n][LdsStrH[n][2]] - status[3][t] * ImLdsStr[n][LdsStrH[n][3]], 2));
							rpc += 1;
						}
						else
						{
							rn[rnc] = status[0][t] / 2.0 * D[n][LdsStrH[n][0]] + status[1][t] / 2.0 * D[n][LdsStrH[n][1]] + status[2][t] / 2.0 * D[n][LdsStrH[n][2]] + status[3][t] / 2.0 * D[n][LdsStrH[n][3]] - status[c][t] / 2.0 * D[n][LdsStrH[n][c]] - 1 / (2 * pow(sigma_n, 2))*(pow(RX_MOD_SYM[n][ss].Real - status[0][t] * ReLdsStr[n][LdsStrH[n][0]] - status[1][t] * ReLdsStr[n][LdsStrH[n][1]] - status[2][t] * ReLdsStr[n][LdsStrH[n][2]] - status[3][t] * ReLdsStr[n][LdsStrH[n][3]], 2) + pow(RX_MOD_SYM[n][ss].Image - status[0][t] * ImLdsStr[n][LdsStrH[n][0]] - status[1][t] * ImLdsStr[n][LdsStrH[n][1]] - status[2][t] * ImLdsStr[n][LdsStrH[n][2]] - status[3][t] * ImLdsStr[n][LdsStrH[n][3]], 2));
							rnc += 1;
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
			RX_LLR_BIT[i][ss] = 0;
			for (int m = 0; m < 12; m++)
			{
				RX_LLR_BIT[i][ss] += D[m][i];
			}
		}	
	}	
	//RX_LLR_BIT to RX_LLR_SYM
	for (int ss = 0; ss < CODE_SYM_LEN; ss++)
	{
		for (int i = 0; i < 16; i++)
		{
			for (int q = 1; q < codeOrder; q++)
			{
				RX_LLR_SYM[i][ss][q - 1] = 0;
				for (int b_p_s = 0; b_p_s < CODE_BIT_PER_SYM; b_p_s++)
				{
					if ((q & (1 << (CODE_BIT_PER_SYM - 1 - b_p_s))) != 0)
					{
						RX_LLR_SYM[i][ss][q - 1] += RX_LLR_BIT[i][ss * CODE_BIT_PER_SYM + b_p_s];
					}
				}
			}
		}
	}
	return 0;
}


int CComm::Decode(void)
{

	// NBLDPC deocding
	for (int i = 0; i < 16; i++)
	{
		DecodeCorrect = NBLDPC.Decoding(RX_LLR_SYM[i], RX_DECODE_SYM[i]);
	}
	for (int i = 0; i < 16; i++)
	{
		for (int s = 0; s < CODE_SYM_LEN; s++)
		{
			for (int b_p_s = 0; b_p_s < CODE_BIT_PER_SYM; b_p_s++)
			{
				RX_DECODE_BIT[i][s * CODE_BIT_PER_SYM + b_p_s]
					= ((RX_DECODE_SYM[i][s] & (1 << (CODE_BIT_PER_SYM - 1 - b_p_s))) == 0) ? 0 : 1;
				//RX_DECODE_BIT[i][s * 4 + b_p_s] = RX_DECODE_BIT[i][s * 4 + b_p_s] * 2 - 1;
			}
		}
	}	
	return 0;
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


double CComm::SetEbN0(CSimulation &sim)
{
	// reset the seed
	Rand.IX = Rand.IY = Rand.IZ = sim.randomseed;
	// reset the PN register
	regPN[0] = 1;	regPN[1] = 0;	regPN[2] = 1;	regPN[3] = 0;
	regPN[4] = 0;	regPN[5] = 0;	regPN[6] = 1;	regPN[7] = 1;
	regPN[8] = 0;	regPN[9] = 0;	regPN[10] = 1;
	// set sigma
	double EbN0 = pow(10.0, sim.EbN0 / 10.0);
	//sigma_n = 0.6;
	sigma_n = 1.0 / sqrt(2 * 1  * EbN0);
	return sigma_n;
}


//int CComm::Err(double simcycle, double& errFrame, double& errBit, double& errSym, double& fer, double& ber, double& ser)
int CComm::Err(CSimulation &sim)
{
	// BER, SER, FER
	double errSym, errBit, errFrame;
	errSym = errBit = errFrame = 0;
	bool flag;
	for (int i = 0; i < 16; i++)
	{
		flag = false;
		for (int j = 0; j < MSG_BIT_LEN; j++)
		{
			if (TX_MSG_BIT[i][j] != RX_DECODE_BIT[i][j])
			{
				errBit++;
				flag = true;
			}
		}
		if (flag)
			errFrame = errFrame + 1;
	}
	sim.errBit += errBit;
	sim.errFrame += errFrame;

	/*if (DecodeCorrect && (errSym != 0))
	{
		sim.U_errSym += errSym;
		sim.U_errBit += errBit;
		sim.U_errFrame += 1;
	}*/

	
	sim.BER = sim.errBit / (sim.simCycle * 16 * MSG_BIT_LEN);
	sim.FER = sim.errFrame / (sim.simCycle * 16);

	/*sim.U_SER = sim.U_errSym / (sim.simCycle * 16 * 2026);
	sim.U_BER = sim.U_errBit / (sim.simCycle * 16 * 2026);
	sim.U_FER = sim.U_errFrame / sim.simCycle;*/
	cout <<sim.simCycle <<"  "<<sim.errBit << "  " << sim.errFrame << "  " << sim.BER << "  " << sim.FER << endl;
	return 0;
}
