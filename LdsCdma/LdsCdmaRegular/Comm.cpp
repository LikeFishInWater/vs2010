#include "Comm.h"


CComm::CComm(void)
	: TX_MSG(nullptr)
	, TX_MOD_SYM(nullptr)
	, RX_MOD_SYM(nullptr)
	, RX_LLR(nullptr)
	, RX_DECODE(nullptr)
	, ReLdsStr(nullptr)
	, ImLdsStr(nullptr)
	, sigma_n(0)
	, randomMsg(0)
{
}

CComm::~CComm(void)
{
}

bool CComm::Initial()
{
	ifstream profile("MatszhN240_120_1_2.txt");
	profile>>K>>N;
	profile>>dv>>dc;
	LdsStr=new int*[N];
	for(int i=0;i<N;i++)
	{
		LdsStr[i]=new int[K];
		for(int j=0;j<K;j++)
		{
			LdsStr[i][j]=0;
		}
	}
	LdsStrV=new int*[K];
	for(int i=0;i<K;i++)
	{
		LdsStrV[i]=new int[dv];
	}
	LdsStrH=new int*[N];
	for(int i=0;i<N;i++)
	{
		LdsStrH[i]=new int[dc];
	}
	status=new int*[dc];
	for(int i=0;i<dc;i++)
	{
		status[i]=new int[int(pow(double(2),double(dc)))];
	}
	int temp;
	for(int i=0;i<int(pow(double(2),double(dc)));i++)
	{
		temp=i;
		for(int j=0;j<dc;j++)
		{
			status[j][i]=2*(temp%2)-1;
			temp=temp/2;
		}
	}

	int rub;
	for(int i=0;i<K+N;i++)
	{
		profile>>rub;
	}
	for(int i=0;i<K;i++)
	{
		for(int j=0;j<dv;j++)
		{
			profile>>rub;
			LdsStrV[i][j]=rub-1;
			LdsStr[rub-1][i]=1;
		}
	}
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<dc;j++)
		{
			profile>>rub;
			LdsStrH[i][j]=rub-1;
		}
	}
	profile.close();
	//K=16;N=12;dc=4;dv=3;

	ite = 1;
	randomMsg = 1;
	randomseed = 173;
	TX_MSG = new int[K];
	TX_MOD_SYM = new CComplex[N];
	RX_MOD_SYM = new CComplex[N];
	RX_LLR = new double[K];
	RX_DECODE = new int[K];

	ReLdsStr = new double*[N];
	ImLdsStr = new double*[N];
	for (int i = 0; i<N; i++)
	{
		ReLdsStr[i] = new double[K];
		ImLdsStr[i] = new double[K];
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < K; j++)
		{
			if (LdsStr[i][j] == 1)
			{
				ReLdsStr[i][j] =  sqrt(1/double(dc))*cos(3.14159 / double(K)* (j+1));
				ImLdsStr[i][j] =  sqrt(1/double(dc))*sin(3.14159 / double(K)* (j+1));
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
	for (int i = 0; i<K; i++)        //TX_MSG_BIT size:K*msg_bit_len
	{
		if (randomMsg)
		{
			TX_MSG[i] = GenPN();
		}
		else
		{
			TX_MSG[i] = 0;
		}
		TX_MSG[i] = 2 * TX_MSG[i] - 1;
	}
	return 0;
}

int CComm::Spread(void)
{
	for (int i = 0; i < N; i++)
	{
		TX_MOD_SYM[i].Real=0;
		TX_MOD_SYM[i].Image=0;
		for(int j=0;j<dc;j++)
		{
			TX_MOD_SYM[i].Real+=TX_MSG[LdsStrH[i][j]] * ReLdsStr[i][LdsStrH[i][j]];
			TX_MOD_SYM[i].Image+=TX_MSG[LdsStrH[i][j]] * ImLdsStr[i][LdsStrH[i][j]];
		}
	}
	return 1;
}

int CComm::Channel_AWGN(void)
{
	// only support bpsk now
	for (int i = 0; i < N; i++)
	{
		RX_MOD_SYM[i].Real = TX_MOD_SYM[i].Real + Rand.Rand_Norm(0, sigma_n);
		RX_MOD_SYM[i].Image = TX_MOD_SYM[i].Image + Rand.Rand_Norm(0, sigma_n);
	}
	return 0;
}

int CComm::DeSpread(void)
{
	double he, ptem, ntem;
	double sum1,sum2,sum3;
	int rpc, rnc;
	double** D;
	double** Dtemp;
	D=new double*[N];
	Dtemp=new double*[N];
	for(int i=0;i<N;i++)
	{
		D[i]=new double[K];
		Dtemp[i]=new double[K];
	}
	double *rp=new double[int(pow(double(2),double(dc-1)))];
	double *rn=new double[int(pow(double(2),double(dc-1)))];

	for (int mm = 0; mm < N; mm++)
		for (int nn = 0; nn < K; nn++)
			D[mm][nn] = 0;
	for (int i = 0; i<ite; i++)
	{
		for (int m = 0; m<K; m++)
		{
			he = 0;
			for (int n = 0; n<N; n++)
				he += D[n][m];
			for (int n = 0; n < N; n++)
			{
				if (LdsStr[n][m])
					D[n][m] = he - D[n][m];
			}
		}
		for (int mm = 0; mm < N; mm++)
			for (int nn = 0; nn < K; nn++)
				Dtemp[mm][nn] = D[mm][nn];
		for (int n = 0; n<N; n++)
		{
			for (int c = 0; c<dc; c++)
			{
				rpc = 0; rnc = 0;
				for (int t = 0; t < int(pow(double(2),double(dc))); t++)
				{
					if (status[c][t] == 1)
					{
						sum1=- status[c][t] / 2.0 * D[n][LdsStrH[n][c]];
						sum2=0;
						sum3=0;
						for(int tt=0;tt<dc;tt++)
						{
							sum1=sum1+status[tt][t] / 2.0 * D[n][LdsStrH[n][tt]];
							sum2=sum2+status[tt][t] * ReLdsStr[n][LdsStrH[n][tt]];
							sum3=sum3+status[tt][t] * ImLdsStr[n][LdsStrH[n][tt]];
						}
						rp[rpc]=sum1- 1 / (2 * pow(sigma_n, 2))*(pow(RX_MOD_SYM[n].Real-sum2,2)+pow(RX_MOD_SYM[n].Image-sum3,2));
						rpc += 1;
					}
					else
					{
						sum1=- status[c][t] / 2.0 * D[n][LdsStrH[n][c]];
						sum2=0;
						sum3=0;
						for(int tt=0;tt<dc;tt++)
						{
							sum1=sum1+status[tt][t] / 2.0 * D[n][LdsStrH[n][tt]];
							sum2=sum2+status[tt][t] * ReLdsStr[n][LdsStrH[n][tt]];
							sum3=sum3+status[tt][t] * ImLdsStr[n][LdsStrH[n][tt]];
						}
						rn[rnc]=sum1- 1 / (2 * pow(sigma_n, 2))*(pow(RX_MOD_SYM[n].Real-sum2,2)+pow(RX_MOD_SYM[n].Image-sum3,2));
						rnc += 1;
					}
				}
				ptem = rp[0];
				ntem = rn[0];
				for (int t = 1; t < int(pow(double(2),double(dc-1))); t++)
				{
					ptem = (ptem > rp[t] ? ptem : rp[t]);
					ntem = (ntem > rn[t] ? ntem : rn[t]);
				}
				Dtemp[n][LdsStrH[n][c]] = ptem - ntem;
			}
		}
		for (int mm = 0; mm < N; mm++)
			for (int nn = 0; nn < K; nn++)
				D[mm][nn] = Dtemp[mm][nn];
	}
	for (int i = 0; i < K; i++)
	{
		RX_LLR[i] = 0;
		for (int m = 0; m < N; m++)
		{
			RX_LLR[i] += D[m][i];
		}
	}
	for(int i=0;i<N;i++)
	{
		delete []D[i];
		delete []Dtemp[i];
	}
	delete []D;
	delete []Dtemp;
	delete []rp;
	delete []rn;
	return 1;
}

int CComm::Decode(void)
{
	for (int i = 0; i < K; i++)
	{
		if(RX_LLR[i]>0)
			RX_DECODE[i] = 1;
		else
			RX_DECODE[i]= -1;			
	}
	return 1;
}

int CComm::err(void)
{
	errBit = 0;
	for (int i = 0; i < K; i++)
	{
		if (RX_DECODE[i] != TX_MSG[i])
			errBit += 1;
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
	sigma_n = 1.0 / sqrt(2 * 1 * EbN0);
	return sigma_n;
}