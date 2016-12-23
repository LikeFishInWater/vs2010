#include "Comm.h"


CComm::CComm(void)
	: TX_MSG(nullptr)
	, TX_MSG_MOD(nullptr)
	, TX_MOD_SYM(nullptr)
	, RX_MOD_SYM(nullptr)
	, decode(nullptr)
	, decodeBit(nullptr)
	, LdsStrC(nullptr)
	, sigma_n(0)
	, randomMsg(0)
{
}

CComm::~CComm(void)
{
}

bool CComm::Initial()
{
	ifstream profile("MatszhN16_12_3_4_sample.txt");
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
	status0=new CComplex*[dc];
	for(int i=0;i<dc;i++)
	{
		status[i]=new int[int(pow(double(4),double(dc)))];
		status0[i]=new CComplex[int(pow(double(4),double(dc)))];
	}
	int temp;
	for(int i=0;i<int(pow(double(4),double(dc)));i++)
	{
		temp=i;
		for(int j=0;j<dc;j++)
		{
			status[j][i]=temp%4;
			temp=temp/4;
		}
	}
	for (int i=0;i<pow(double(4),double(dc));i++)
		for (int j=0;j<dc;j++)
		{
			status0[j][i].Real=cos(status[j][i]*3.14159/2.0);
			status0[j][i].Image=sin(status[j][i]*3.14159/2.0);
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
	LdsStrC = new CComplex*[N];
	for (int i = 0; i<N; i++)
	{
		LdsStrC[i] = new CComplex[K];
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < K; j++)
		{
			if (LdsStr[i][j] == 1)
			{
				LdsStrC[i][j].Real =  sqrt(1/double(dc))*cos(3.14159 / double(K)* double(j+1));
				LdsStrC[i][j].Image =  sqrt(1/double(dc))*sin(3.14159 / double(K)* double(j+1));
			}
			else
			{
				LdsStrC[i][j].Real = 0;
				LdsStrC[i][j].Image  = 0;
			}
		}
	}

	ite = 7;
	randomMsg = 1;
	randomseed = 173;
	TX_MSG = new int[2*K];
	TX_MSG_MOD=new CComplex[K];
	TX_MOD_SYM = new CComplex[N];
	RX_MOD_SYM = new CComplex[N];
	decode = new int[K];
	decodeBit = new int[2*K];

	

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
	for (int i = 0; i<2*K; i++)        //TX_MSG_BIT size:K*msg_bit_len
	{
		if (randomMsg)
		{
			TX_MSG[i] = GenPN();
		}
		else
		{
			TX_MSG[i] = 0;
		}
		// TX_MSG[i] = 2 * TX_MSG[i] - 1;
	}
	for (int i=0;i<K;i++)
	{
		if (TX_MSG[2*i]==0)
		{
			if (TX_MSG[2*i+1]==0)
			{TX_MSG_MOD[i].Real=1;TX_MSG_MOD[i].Image=0;}
			else
			{TX_MSG_MOD[i].Real=0;TX_MSG_MOD[i].Image=1;}
		}
		else
		{
			if (TX_MSG[2*i+1]==0)
			{TX_MSG_MOD[i].Real=0;TX_MSG_MOD[i].Image=-1;}
			else
			{TX_MSG_MOD[i].Real=-1;TX_MSG_MOD[i].Image=0;}
		}
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
			TX_MOD_SYM[i]=TX_MOD_SYM[i]+TX_MSG_MOD[LdsStrH[i][j]] * LdsStrC[i][LdsStrH[i][j]];
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
	double **mult;
	mult=new double*[K];
	for (int k=0;k<K;k++)
	{
		mult[k]=new double[4];
	}
	double sum0,temp1;
	CComplex temp2;
	double*** D;
	double*** Dtemp;
	D=new double**[N];
	Dtemp=new double**[N];
	for(int i=0;i<N;i++)
	{
		D[i]=new double*[K];
		Dtemp[i]=new double*[K];
		for (int j=0;j<K;j++)
		{
			D[i][j]=new double[4];
			Dtemp[i][j]=new double[4];
		}
	}

	for (int mm = 0; mm < N; mm++)
		for (int nn = 0; nn < K; nn++)
			for(int tt=0;tt<4;tt++)
				D[mm][nn][tt] = 0;
	for(int mm=0;mm<K;mm++)
		for(int nn=0;nn<dv;nn++)
			for(int tt=0;tt<4;tt++)
				D[LdsStrV[mm][nn]][mm][tt]=0.25;
	for (int i = 0; i<ite; i++)
	{
		for (int mm = 0; mm < N; mm++)
			for (int nn = 0; nn < K; nn++)
				for(int tt=0;tt<4;tt++)
					Dtemp[mm][nn][tt] = D[mm][nn][tt];
		for (int m = 0; m<K; m++)
		{
			for (int nn=0;nn<dv;nn++)
			{
				for(int tt=0;tt<4;tt++)
					Dtemp[LdsStrV[m][nn]][m][tt]=0.25;
				for(int tt=0;tt<dv;tt++)
				{
					if(tt!=nn)
					{
						for(int ttt=0;ttt<4;ttt++)
							Dtemp[LdsStrV[m][nn]][m][ttt]=Dtemp[LdsStrV[m][nn]][m][ttt]*D[LdsStrV[m][tt]][m][ttt];
					}
				}
				sum0=0;
				for(int ttt=0;ttt<4;ttt++)
					sum0+=Dtemp[LdsStrV[m][nn]][m][ttt];
				for(int ttt=0;ttt<4;ttt++)
					Dtemp[LdsStrV[m][nn]][m][ttt]=Dtemp[LdsStrV[m][nn]][m][ttt]/sum0;
			}
		}
		for (int mm = 0; mm < N; mm++)
			for (int nn = 0; nn < K; nn++)
				for(int tt=0;tt<4;tt++)
					D[mm][nn][tt] = Dtemp[mm][nn][tt];
		for (int mm = 0; mm < N; mm++)
			for (int nn = 0; nn < K; nn++)
				for(int tt=0;tt<4;tt++)
					Dtemp[mm][nn][tt] = 0;
		for (int m = 0; m<N; m++)
		{
			for (int n = 0; n<dc; n++)
			{
				for (int t = 0; t < int(pow(double(4),double(dc))); t++)
				{
					temp1=1/D[m][LdsStrH[m][n]][status[n][t]];
                    temp2=RX_MOD_SYM[m];
                    for (int tt=0;tt<dc;tt++)
					{
                        temp1*=D[m][LdsStrH[m][tt]][status[tt][t]];
                        temp2=temp2-status0[tt][t]*LdsStrC[m][LdsStrH[m][tt]];
					}
					Dtemp[m][LdsStrH[m][n]][status[n][t]]+=temp1*(1/(sqrt(2*3.14159)*sigma_n)*exp(-1/2.0/pow(sigma_n,2)*pow(temp2.Cabs(),2)));
				}
			}
		}
		for (int mm = 0; mm < N; mm++)
			for (int nn = 0; nn < K; nn++)
				for(int tt=0;tt<4;tt++)
					D[mm][nn][tt] = Dtemp[mm][nn][tt];
	}
	for (int i = 0; i < K; i++)
	{
		for(int tt=0;tt<4;tt++)
		{
			mult[i][tt]=1;
			for (int nn=0;nn<dv;nn++)
			{
				mult[i][tt]*=D[LdsStrV[i][nn]][i][tt];
			}
		}
		decode[i]=0;
		for(int tt=1;tt<4;tt++)
		{
			if(mult[i][tt]>mult[i][decode[i]])
				decode[i]=tt;
		}		
	}
	for(int i=0;i<N;i++)
	{
		delete []D[i];
		delete []Dtemp[i];
	}
	for(int i=0;i<K;i++)
	{
		delete []mult[i];
	}
	delete []mult;
	delete []D;
	delete []Dtemp;
	return 1;
}

int CComm::Decode(void)
{
	for (int i = 0; i < K; i++)
	{
		switch(decode[i])
		{
		case 0:decodeBit[2*i]=0;decodeBit[2*i+1]=0;break;
		case 1:decodeBit[2*i]=0;decodeBit[2*i+1]=1;break;
		case 2:decodeBit[2*i]=1;decodeBit[2*i+1]=1;break;
		case 3:decodeBit[2*i]=1;decodeBit[2*i+1]=0;break;
		default:break;
		}		
	}
	int temp[32];
	for (int i=0;i<2*K;i++)
	{
		temp[i]=decodeBit[i]-TX_MSG[i];
	}
	return 1;
}

int CComm::err(void)
{
	errBit = 0;
	for (int i = 0; i < 2*K; i++)
	{
		if (decodeBit[i] != TX_MSG[i])
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