#include "NBLDPC.h"


CNBLDPC::CNBLDPC(void)
	: GFq(0)
	, CodeLen(0)
	, ChkLen(0)
	, PuncturePosition(nullptr)
	, PunctureLen(0)
	, maxVarDegree(0)
	, maxChkDegree(0)
	, VarDegree(nullptr)
	, ChkDegree(nullptr)
	, VarLink(nullptr)
	, ChkLink(nullptr)
	, VarLinkGFe(NULL)
	, ChkLinkGFe(NULL)
	, maxIter(0)
	, L_v2c(NULL)
	, L_c2v(NULL)
	, LLRV(NULL)
	, L_post(NULL)
	, VarLinkDc(NULL)
	, ChkLinkDv(NULL)
	, L_sigma(NULL)
	, L_rho(NULL)
	, EncodeExchangeColSource(NULL)
	, EncodeExchangeColDest(NULL)
	, EncodeExchangeColNum(0)
	, EncodeVarLink(NULL)
	, EncodeVarLinkGFe(NULL)
	, EncodeVarDv(NULL)
	, EMS_Nm(0)
	, EMS_Nc(0)
	, DecodeMethod(0)
	, EMS_sort_L_v2c(NULL)
	, EMS_sort_Entr_v2c(NULL)
	, EMS_L_c2v(NULL)
	, EMS_Correction_Factor(0)
	, EMS_Correction_Offset(0)
	, TEMS_Nr(0)
	, TEMS_Nc(0)
	, TEMS_deltaU(NULL)
	, TEMS_deltaW(NULL)
	, TEMS_Beta(NULL)
	, TEMS_Eta(NULL)
	, TEMS_Min(NULL)
	, TEMS_inConf_Nr(NULL)
	, TEMS_Eta_Candidate(NULL)
	, TEMS_Lc2v(NULL)
	, TEMS_isUpdate(NULL)
	, TEMS_Syndrome(0)
	, TEMS_Factor(0)
	, TEMS_Offset(0)
	, TEMS_isSelected(NULL)
{
}


CNBLDPC::~CNBLDPC(void)
{
}


//bool CNBLDPC::Initial(int q, string NBFileName, int PunctureVarDegree, int maxiter, int randMsg, int deocdemethod, int ems_nm, int ems_nc, double ems_factor, double ems_offset)
bool CNBLDPC::Initial(CSimulation &sim)
{
	GFq = sim.GFq;
	maxIter = sim.maxIter;
	cout << "1 " << sim.NonBinaryFileName << "1"<<endl;
	ifstream fnb(sim.NonBinaryFileName);
	if (!fnb.is_open())
	{
		cerr << "Cannot open " << sim.NonBinaryFileName << endl;
		exit(-1);
	}
	fnb >> CodeLen >> ChkLen >> GFq;
	fnb >> maxVarDegree >> maxChkDegree;
	VarDegree = new int[CodeLen];
	ChkDegree = new int[ChkLen];

	PunctureLen = 0;
	VarLink = new int*[CodeLen];
	VarLinkGFe = new int*[CodeLen];
	VarLinkDc = new int*[CodeLen];
	for (int col = 0; col < CodeLen; col++)
	{
		fnb >> VarDegree[col];
		VarLink[col] = new int[VarDegree[col]];
		VarLinkGFe[col] = new int[VarDegree[col]];
		VarLinkDc[col] = new int[VarDegree[col]];
		if (VarDegree[col] == sim.PuntureVarDegree)
		{
			PunctureLen++;
		}
	}
	PuncturePosition = new int[PunctureLen];
	int PunctureIndex = 0;
	for (int col = 0; col < CodeLen; col++)
	{
		if (VarDegree[col] == sim.PuntureVarDegree)
		{
			PuncturePosition[PunctureIndex] = col;
			PunctureIndex++;
		}
	}
	ChkLink = new int*[ChkLen];
	ChkLinkGFe = new int*[ChkLen];
	ChkLinkDv = new int*[ChkLen];
	for (int row = 0; row < ChkLen; row++)
	{
		fnb >> ChkDegree[row];
		ChkLink[row] = new int[ChkDegree[row]];
		ChkLinkGFe[row] = new int[ChkDegree[row]];
		ChkLinkDv[row] = new int[ChkDegree[row]];
	}

	for (int col = 0; col < CodeLen; col++)
	{
		for (int d = 0; d < VarDegree[col]; d++)
		{
			fnb >> VarLink[col][d] >> VarLinkGFe[col][d];
			VarLink[col][d] --;
		}
	}
	for (int row = 0; row < ChkLen; row++)
	{
		for (int d = 0; d < ChkDegree[row]; d++)
		{
			fnb >> ChkLink[row][d] >> ChkLinkGFe[row][d];
			ChkLink[row][d] --;
		}
	}
	fnb.close();

	// initial GF
	GF.Initial(GFq);


	// allocate memory
	LLRV = new double[GFq - 1];
	L_post = new double *[CodeLen];
	for (int col = 0; col < CodeLen; col++)
	{
		L_post[col] = new double[GFq - 1];
	}
	L_v2c = new double **[CodeLen];
	for (int col = 0; col < CodeLen; col++)
	{
		L_v2c[col] = new double *[VarDegree[col]];
		for (int d = 0; d < VarDegree[col]; d++)
		{
			L_v2c[col][d] = new double[GFq - 1];
		}
	}
	L_c2v = new double **[ChkLen];
	for (int row = 0; row < ChkLen; row++)
	{
		L_c2v[row] = new double*[ChkDegree[row]];
		for (int d = 0; d < ChkDegree[row]; d++)
		{
			L_c2v[row][d] = new double[GFq - 1];
		}
	}

	// create the link of VarLinkDc & ChkLinkDv
	for (int col = 0; col < CodeLen; col++)
	{
		for (int dv = 0; dv < VarDegree[col]; dv++)
		{
			int row = VarLink[col][dv];
			for (int dc = 0; dc < ChkDegree[row]; dc++)
			{
				if (ChkLink[row][dc] == col)
				{
					VarLinkDc[col][dv] = dc;
				}
			}
		}
	}
	for (int row = 0; row < ChkLen; row++)
	{
		for (int dc = 0; dc < ChkDegree[row]; dc++)
		{
			int col = ChkLink[row][dc];
			for (int dv = 0; dv < VarDegree[col]; dv++)
			{
				if (VarLink[col][dv] == row)
				{
					ChkLinkDv[row][dc] = dv;
				}
			}
		}
	}

	// allocate memory for iterative decoding
	DecodeMethod = sim.decodeMethod;
	if (DecodeMethod == BP_DECODE)
	{
		L_sigma = new double[GFq - 1];
		L_rho = new double[GFq - 1];
	}
	else if (DecodeMethod == EMS_DECODE)
	{
		EMS_Nm = sim.ems_nm;
		EMS_Nc = sim.ems_nc;
		EMS_Correction_Factor = sim.ems_factor;
		EMS_Correction_Offset = sim.ems_offset;

		if (EMS_Nm > GFq)
		{
			cerr << "EMS configuration error! EMS_Nm is too large!" << endl;
			exit(-1);
		}

		EMS_sort_L_v2c = new double *[maxChkDegree];
		EMS_sort_Entr_v2c = new int*[maxChkDegree];
		for (int d = 0; d < maxChkDegree; d++)
		{
			EMS_sort_L_v2c[d] = new double[GFq];
			EMS_sort_Entr_v2c[d] = new int[GFq];
		}
		EMS_L_c2v = new double[GFq];
	}
	else if (DecodeMethod == T_EMS_DECODE)
	{
		TEMS_Nr = sim.tems_nr;
		TEMS_Nc = sim.tems_nc;
		TEMS_Factor = sim.tems_factor;
		TEMS_Offset = sim.tems_offset;

		TEMS_deltaU = new double*[maxChkDegree];
		for (int d = 0; d < maxChkDegree; d++)
		{
			TEMS_deltaU[d] = new double[GFq];
		}
		TEMS_Beta = new int[maxChkDegree];
		TEMS_Min = new int*[GFq];
		TEMS_inConf_Nr = new bool*[GFq];
		for (int q = 0; q < GFq; q++)
		{
			TEMS_Min[q] = new int[maxChkDegree];
			TEMS_inConf_Nr[q] = new bool[maxChkDegree];
		}
		TEMS_deltaW = new double[GFq];
		TEMS_Eta = new int*[GFq];
		for (int q = 0; q < GFq; q++)
		{
			TEMS_Eta[q] = new int[maxChkDegree];
		}
		TEMS_Eta_Candidate = new int[maxChkDegree];
		TEMS_Lc2v = new double[GFq];
		TEMS_isUpdate = new bool[GFq];
		TEMS_isSelected = new bool[GFq];
	}
	// initial encode
	if (sim.randomMsg)
	{
		InitialEncode();
	}

	return true;
}


int CNBLDPC::Decoding_BP(double** L_ch, int* DecodeOutput)
{
	//	ofstream decoderecord("decode.txt");
	// initial
	for (int col = 0; col < CodeLen; col++)
	{
		for (int d = 0; d < VarDegree[col]; d++)
		{
			CopyLLRVector(L_v2c[col][d], L_ch[col]);
		}
	}
	for (int row = 0; row < ChkLen; row++)
	{
		for (int d = 0; d < ChkDegree[row]; d++)
		{
			ClearLLRVector(L_c2v[row][d]);
		}
	}
	// iterative decoding
	int iter_number = 0;
	while (iter_number++ < maxIter)
	{
		// tentative decoding
		for (int col = 0; col < CodeLen; col++)
		{
			CopyLLRVector(L_post[col], L_ch[col]);
			int row, dc;
			for (int d = 0; d < VarDegree[col]; d++)
			{
				row = VarLink[col][d];
				dc = VarLinkDc[col][d];
				AddLLRVector(L_post[col], L_post[col], L_c2v[row][dc]);
			}
			DecodeOutput[col] = DecideLLRVector(L_post[col]);
			//		decoderecord << DecodeOutput[col] << " " << L_post[col][DecodeOutput[col] - 1] << " ";
		}
		//	decoderecord << endl;

		/*
		bool isZero = true;
		for(int col = 0; col < CodeLen; col ++)
		{
		isZero = true;
		for(int q = 0; q < GFq - 1; q ++)
		{
		decoderecord << L_post[col][q] << " \t";
		if(L_post[col][q] != 0)
		{
		isZero = false;
		}
		}
		if(isZero)
		{
		int a = 0;
		}
		decoderecord << endl;
		}
		decoderecord << endl;
		*/

		// check whether satisfy the parity matrix
		bool decode_correct = true;
		for (int row = 0; row < ChkLen; row++)
		{
			int col, multip, syndrome;
			syndrome = 0;
			for (int d = 0; d < ChkDegree[row]; d++)
			{
				col = ChkLink[row][d];
				multip = GF.GFMultiply(ChkLinkGFe[row][d], DecodeOutput[col]);
				syndrome = GF.GFAdd(syndrome, multip);
			}
			if (syndrome != 0)
			{
				decode_correct = false;
				break;
			}
		}
		if (decode_correct)
		{
			return decode_correct;
		}

		// message from var to check
		for (int col = 0; col < CodeLen; col++)
		{
			int row, dc;
			for (int dv = 0; dv < VarDegree[col]; dv++)
			{
				row = VarLink[col][dv];
				dc = VarLinkDc[col][dv];
				MinusLLRVector(L_v2c[col][dv], L_post[col], L_c2v[row][dc]);
			}
		}

		// message from check to var
		for (int row = 0; row < ChkLen; row++)
		{
			for (int d = 0; d < ChkDegree[row]; d++)
			{
				ClearLLRVector(L_sigma);
				ClearLLRVector(L_rho);
				L_Back(row, d - 1);
				L_Forward(row, d + 1);
				int A1, A2;
				A1 = A2 = GF.GFInverse(ChkLinkGFe[row][d]);
				if (d == 0)
				{
					A1 = 0;
				}
				else if (d == ChkDegree[row] - 1)
				{
					A2 = 0;
				}
				LLR_BoxPlus(L_c2v[row][d], L_sigma, L_rho, A1, A2);
			}
		}
	}
	return 0;
}


int CNBLDPC::CopyLLRVector(double* L_d, double* L_s)
{
	for (int q = 0; q < GFq - 1; q++)
	{
		L_d[q] = L_s[q];
	}
	return 0;
}


int CNBLDPC::ClearLLRVector(double* L)
{
	for (int q = 0; q < GFq - 1; q++)
	{
		L[q] = 0;
	}
	return 0;
}


int CNBLDPC::AddLLRVector(double* L, double* L1, double* L2)
{
	for (int q = 0; q < GFq - 1; q++)
	{
		L[q] = L1[q] + L2[q];
	}
	return 0;
}


int CNBLDPC::DecideLLRVector(double* LLR)
{
	double max = 0;
	int alpha_i;
	for (int q = 0; q < GFq - 1; q++)
	{
		if (LLR[q] > max)
		{
			max = LLR[q];
			alpha_i = q + 1;
		}
	}
	if (max <= 0)
	{
		return 0;
	}
	else
	{
		return alpha_i;
	}
}


int CNBLDPC::MinusLLRVector(double* L, double* L1, double* L2)
{
	for (int q = 0; q < GFq - 1; q++)
	{
		L[q] = L1[q] - L2[q];
	}
	return 0;
}


int CNBLDPC::L_Back(int row, int l)
{
	int A1, A2;
	int col, dv;
	if (l < 0)
	{
		return 0;
	}
	else if (l == 0)
	{
		A1 = 0;
		A2 = ChkLinkGFe[row][l];
		col = ChkLink[row][l];
		dv = ChkLinkDv[row][l];
		LLR_BoxPlus(L_sigma, L_sigma, L_v2c[col][dv], A1, A2);
	}
	else
	{
		L_Back(row, l - 1);
		A1 = 1;
		A2 = ChkLinkGFe[row][l];
		col = ChkLink[row][l];
		dv = ChkLinkDv[row][l];
		LLR_BoxPlus(L_sigma, L_sigma, L_v2c[col][dv], A1, A2);
	}
	return 1;
}


int CNBLDPC::L_Forward(int row, int l)
{
	int A1, A2;
	int col, dv;
	if (l >= ChkDegree[row])
	{
		return 0;
	}
	else if (l == ChkDegree[row] - 1)
	{
		A1 = 0;
		A2 = ChkLinkGFe[row][l];
		col = ChkLink[row][l];
		dv = ChkLinkDv[row][l];
		LLR_BoxPlus(L_rho, L_rho, L_v2c[col][dv], A1, A2);
	}
	else
	{
		L_Forward(row, l + 1);
		A1 = 1;
		A2 = ChkLinkGFe[row][l];
		col = ChkLink[row][l];
		dv = ChkLinkDv[row][l];
		LLR_BoxPlus(L_rho, L_rho, L_v2c[col][dv], A1, A2);
	}
	return 1;
}

/*
int CNBLDPC::LLR_BoxPlus(double* L, double* L1, double* L2, int A1, int A2)
{
if(A1 == 0)
{
int v2;
for(int alpha_i = 1; alpha_i < GFq; alpha_i ++)
{
v2 = GF.GFMultiply(GF.GFInverse(A2), alpha_i);
//			LLRV[alpha_i - 1] = (L2[v2 - 1] > LLR_MAX)? LLR_MAX : L2[v2 - 1];
if(L2[v2 - 1] > LLR_MAX)
{
cout << '\n' << "\t1\t" << endl;
system("pause");
}
LLRV[alpha_i - 1] =  L2[v2 - 1];
}
CopyLLRVector(L, LLRV);
}
else if(A2 == 0)
{
int v1;
for(int alpha_i = 1; alpha_i < GFq; alpha_i ++)
{
v1 = GF.GFMultiply(GF.GFInverse(A1), alpha_i);
//			LLRV[alpha_i - 1] = (L1[v1 - 1] > LLR_MAX)? LLR_MAX : L1[v1 - 1];
if(L1[v1 - 1] > LLR_MAX)
{
cout << '\n' << "\t2\t" << endl;
system("pause");
}
LLRV[alpha_i - 1] = L1[v1 - 1];
}
CopyLLRVector(L, LLRV);
}
else
{
long double sum1, sum2;
int v1, v2;
// down
sum2 = 0;
for(int x = 0; x < GFq; x ++)
{
v1 = x;
v2 = GF.GFMultiply(GF.GFInverse(A2), GF.GFMultiply(A1, x));
if(v1 == 0 || v2 == 0)
{
sum2 += 1;
}
else
{

if(L1[v1 - 1] + L2[v2 - 1] > LLR_MAX)
{
cout << '\n' << "\t3\t" << endl;
system("pause");
}

sum2 += (L1[v1 - 1] + L2[v2 - 1] > LLR_MAX)?
exp(LLR_MAX * 1.0) : exp(L1[v1 - 1] + L2[v2 - 1]);
}
}
for(int alpha_i = 1; alpha_i < GFq; alpha_i ++)
{
// up
sum1 = 0;
for(int x = 0; x < GFq; x ++)
{
// up
v1 = x;
v2 = GF.GFMultiply(GF.GFInverse(A2), GF.GFAdd(alpha_i, GF.GFMultiply(x, A1)));
if(v1 == 0)
{
//					sum1 += (L2[v2 - 1] > LLR_MAX)? exp(LLR_MAX * 1.0) : exp(L2[v2 - 1]);
if(L2[v2 - 1] > LLR_MAX)
{
cout << '\n' << "\t4\t" << endl;
system("pause");
}

sum1 += exp(L2[v2 - 1]);
}
else if(v2 == 0)
{
//					sum1 += (L1[v1 - 1] > LLR_MAX)? exp(LLR_MAX * 1.0) : exp(L1[v1 - 1]);

if(L1[v1 - 1] > LLR_MAX)
{
cout << '\n' << "\t5\t" << endl;
system("pause");
}

sum1 += exp(L1[v1 - 1]);
}
else
{

if(L1[v1 - 1] + L2[v2 - 1] > LLR_MAX)
{
cout << '\n' << "\t6\t" << endl;
system("pause");
}

sum1 += (L1[v1 - 1] + L2[v2 - 1] > LLR_MAX)?
exp(LLR_MAX * 1.0) : exp(L1[v1 - 1] + L2[v2 - 1]);
}
}
LLRV[alpha_i - 1] = log(sum1 / sum2);
if(LLRV[alpha_i - 1] > LLR_MAX)
{
cout << '\n' << "\t7\t" << endl;
system("pause");
LLRV[alpha_i - 1] = LLR_MAX;
}
else if(LLRV[alpha_i - 1] < -LLR_MAX)
{
cout << '\n' << "\t8\t" << endl;
system("pause");
LLRV[alpha_i - 1] = -LLR_MAX;
}
}
CopyLLRVector(L, LLRV);
}
return 0;
}
*/

int CNBLDPC::LLR_BoxPlus(double* L, double* L1, double* L2, int A1, int A2)
{
	if (A1 == 0)
	{
		int v2;
		for (int alpha_i = 1; alpha_i < GFq; alpha_i++)
		{
			v2 = GF.GFMultiply(GF.GFInverse(A2), alpha_i);
			LLRV[alpha_i - 1] = L2[v2 - 1];
		}
		CopyLLRVector(L, LLRV);
	}
	else if (A2 == 0)
	{
		int v1;
		for (int alpha_i = 1; alpha_i < GFq; alpha_i++)
		{
			v1 = GF.GFMultiply(GF.GFInverse(A1), alpha_i);
			LLRV[alpha_i - 1] = L1[v1 - 1];
		}
		CopyLLRVector(L, LLRV);
	}
	else
	{
		long double sum1, sum2;
		//		long double sum1_prime, sum2_prime;
		int v1, v2;
		// down
		sum2 = 0;
		//		sum2_prime = 0;
		for (int x = 0; x < GFq; x++)
		{
			v1 = x;
			v2 = GF.GFMultiply(GF.GFInverse(A2), GF.GFMultiply(A1, x));
			//			if(v1 == 0 || v2 == 0)
			//			{
			//				sum2_prime = 0;
			//				sum2 += 1;
			//			}
			//			else
			if (v1 != 0 && v2 != 0)
			{
				//				sum2_prime = log(sum2);
				if (sum2 > L1[v1 - 1] + L2[v2 - 1])
				{
					sum2 = sum2 + log(1 + exp(-1 * (sum2 - (L1[v1 - 1] + L2[v2 - 1]))));
				}
				else
				{
					sum2 = L1[v1 - 1] + L2[v2 - 1] + log(1 + exp(-1 * (L1[v1 - 1] + L2[v2 - 1] - sum2)));
				}
			}
		}
		for (int alpha_i = 1; alpha_i < GFq; alpha_i++)
		{
			// up
			sum1 = 0;
			//			sum1_prime = 0;
			v1 = GF.GFMultiply(GF.GFInverse(A1), alpha_i);
			v2 = GF.GFMultiply(GF.GFInverse(A2), alpha_i);
			if (L1[v1 - 1] > L2[v2 - 1])
			{
				sum1 = L1[v1 - 1] + log(1 + exp(-1 * (L1[v1 - 1] - L2[v2 - 1])));
			}
			else
			{
				sum1 = L2[v2 - 1] + log(1 + exp(-1 * (L2[v2 - 1] - L1[v1 - 1])));
			}
			for (int x = 0; x < GFq; x++)
			{
				// up
				v1 = x;
				v2 = GF.GFMultiply(GF.GFInverse(A2), GF.GFAdd(alpha_i, GF.GFMultiply(x, A1)));
				if (v1 != 0 && v2 != 0)
				{
					//					sum1_prime = log(sum1);
					if (sum1 > L1[v1 - 1] + L2[v2 - 1])
					{
						sum1 = sum1 + log(1 + exp(-1 * (sum1 - (L1[v1 - 1] + L2[v2 - 1]))));
					}
					else
					{
						sum1 = L1[v1 - 1] + L2[v2 - 1] + log(1 + exp(-1 * (L1[v1 - 1] + L2[v2 - 1] - sum1)));
					}
				}
			}
			LLRV[alpha_i - 1] = sum1 - sum2;
		}
		CopyLLRVector(L, LLRV);
	}
	return 0;
}


int CNBLDPC::InitialEncode(void)
{
	// allocate matrix H
	int ** H = new int*[ChkLen];
	for (int row = 0; row < ChkLen; row++)
	{
		H[row] = new int[CodeLen];
	}
	// reset value of H
	for (int row = 0; row < ChkLen; row++)
	{
		for (int col = 0; col < CodeLen; col++)
		{
			H[row][col] = 0;
		}
	}
	for (int row = 0; row < ChkLen; row++)
	{
		for (int d = 0; d < ChkDegree[row]; d++)
		{
			H[row][ChkLink[row][d]] = ChkLinkGFe[row][d];
		}
	}

	// allocate memory for encode exchange
	EncodeExchangeColSource = new int[CodeLen];
	EncodeExchangeColDest = new int[CodeLen];
	EncodeExchangeColNum = 0;

	// GaussEliminate
	/*
	ofstream OriginH("OriginH.txt");
	for(int row = 0; row < ChkLen; row ++)
	{
	for(int col = 0; col < CodeLen; col ++)
	{
	OriginH << H[row][col] << " ";
	}
	OriginH << endl;
	}
	OriginH.close();
	*/
	GaussEliminate(H, EncodeExchangeColSource, EncodeExchangeColDest, EncodeExchangeColNum);
	/*
	ofstream tempH("EliminateH.txt");
	for(int row = 0; row < ChkLen; row ++)
	{
	for(int col = 0; col < CodeLen; col ++)
	{
	tempH << H[row][col] << " ";
	}
	tempH << endl;
	}
	tempH.close();
	*/

	// allocate memory for encode var node
	EncodeVarDv = new int[ChkLen];
	EncodeVarLink = new int*[ChkLen];
	EncodeVarLinkGFe = new int*[ChkLen];

	int all = 0;
	for (int p = 0; p < ChkLen; p++)
	{
		EncodeVarDv[p] = 0;
		for (int col = 0; col < CodeLen - ChkLen + p; col++)
		{
			if (H[p][col] != 0)
			{
				EncodeVarDv[p] ++;
			}
		}
		all += EncodeVarDv[p];
		EncodeVarLink[p] = new int[EncodeVarDv[p]];
		EncodeVarLinkGFe[p] = new int[EncodeVarDv[p]];
	}
	for (int p = 0; p < ChkLen; p++)
	{
		int d = 0;
		for (int col = 0; col < CodeLen - ChkLen + p; col++)
		{
			if (H[p][col] != 0)
			{
				EncodeVarLink[p][d] = col;
				EncodeVarLinkGFe[p][d] = H[p][col];
				d++;
			}
		}
	}

	// release memory of H
	for (int row = 0; row < ChkLen; row++)
	{
		delete[] H[row];
	}
	delete[] H;

	return 0;
}


int CNBLDPC::GaussEliminate(int** H, int* ExchangeColSource, int* ExchangeColDest, int& ExchangeColNum)
{
	for (int row = ChkLen - 1; row >= 0; row--)
	{
		// eliminate row 'row'
		int col = row + CodeLen - ChkLen;
		// make the H[row][col] != 0
		if (H[row][col] == 0)
		{
			bool exchanged = false;
			// search up
			if (!exchanged)
			{
				for (int row_up = row - 1; row_up >= 0; row_up--)
				{
					if (H[row_up][col] != 0)
					{
						ExchRow(H, row, row_up);
						exchanged = true;
						break;
					}
				}
			}
			// if cannot find nonzero ele, search left
			if (!exchanged)
			{
				for (int col_left = col - 1; col_left >= 0; col_left--)
				{
					if (H[row][col_left] != 0)
					{
						ExchCol(H, col, col_left);
						EncodeExchangeColSource[EncodeExchangeColNum] = col;
						EncodeExchangeColDest[EncodeExchangeColNum] = col_left;
						EncodeExchangeColNum++;
						exchanged = true;
						break;
					}
				}
			}
			if (!exchanged)
			{
				cerr << "NB matrix is not full rank" << endl;
				exit(-1);
			}
		}
		// make the H[0:row-1][col] = 0
		int h_inverse = GF.GFInverse(H[row][col]);
		for (int row_up = row - 1; row_up >= 0; row_up--)
		{
			int x;
			if (H[row_up][col] != 0)
			{
				x = GF.GFMultiply(h_inverse, H[row_up][col]);
				AddRow(H, row, row_up, x);
			}
		}
		// make the H[row][col] = 1
		for (int col_left = 0; col_left <= col; col_left++)
		{
			H[row][col_left] = GF.GFMultiply(H[row][col_left], h_inverse);
		}
	}
	return 0;
}


int CNBLDPC::ExchRow(int** H, int source, int destination)
{
	int temp;
	for (int col = 0; col < CodeLen; col++)
	{
		temp = H[source][col];
		H[source][col] = H[destination][col];
		H[destination][col] = temp;
		//		H[source][col] =  H[source][col] + H[destination][col];
		//		H[destination][col] = H[source][col] - H[destination][col];
		//		H[source][col] = H[source][col] - H[destination][col];
	}
	return 0;
}


int CNBLDPC::AddRow(int** H, int source, int destination, int multifactor)
{
	for (int col = 0; col < CodeLen; col++)
	{
		H[destination][col] = GF.GFAdd(H[destination][col], GF.GFMultiply(multifactor, H[source][col]));
	}
	return 0;
}


int CNBLDPC::ExchCol(int** H, int source, int destination)
{
	int temp;
	for (int row = 0; row < ChkLen; row++)
	{
		temp = H[row][source];
		H[row][source] = H[row][destination];
		H[row][destination] = temp;
	}
	return 0;
}


int CNBLDPC::Encode(int* msg_sym, int* code_sym)
{
	// put information symbol
	for (int col = 0; col < CodeLen - ChkLen; col++)
	{
		code_sym[col] = msg_sym[col];
	}
	for (int col = CodeLen - ChkLen; col < CodeLen; col++)
	{
		code_sym[col] = 0;
	}

	// encode
	for (int p = 0; p < ChkLen; p++)
	{
		int col = CodeLen - ChkLen + p;
		int v, h, x;
		for (int d = 0; d < EncodeVarDv[p]; d++)
		{
			v = EncodeVarLink[p][d];
			h = EncodeVarLinkGFe[p][d];
			x = code_sym[v];
			code_sym[col] = GF.GFAdd(code_sym[col], GF.GFMultiply(h, x));
		}
	}

	// permute
	int temp, col_s, col_d;
	for (int num = EncodeExchangeColNum - 1; num >= 0; num--)
	{
		col_s = EncodeExchangeColSource[num];
		col_d = EncodeExchangeColDest[num];
		temp = code_sym[col_s];
		code_sym[col_s] = code_sym[col_d];
		code_sym[col_d] = temp;
	}
	for (int col = 0; col < CodeLen - ChkLen; col++)
	{
		msg_sym[col] = code_sym[col];
	}
	/*
	ofstream encoderecord("encode.txt");
	for(int col = 0; col < CodeLen; col ++)
	{
	encoderecord << code_sym[col] << " ";
	}
	encoderecord << endl;
	encoderecord.close();
	*/	return 0;
}


int CNBLDPC::Decoding_EMS(double** L_ch, int* DecodeOutput)
{
	// initial
	for (int col = 0; col < CodeLen; col++)
	{
		for (int d = 0; d < VarDegree[col]; d++)
		{
			CopyLLRVector(L_v2c[col][d], L_ch[col]);
		}
	}
	for (int row = 0; row < ChkLen; row++)
	{
		for (int d = 0; d < ChkDegree[row]; d++)
		{
			ClearLLRVector(L_c2v[row][d]);
		}
	}
	// iterative decoding
	int iter_number = 0;
	while (iter_number++ < maxIter)
	{
		// tentative decoding
		for (int col = 0; col < CodeLen; col++)
		{
			CopyLLRVector(L_post[col], L_ch[col]);
			int row, dc;
			for (int d = 0; d < VarDegree[col]; d++)
			{
				row = VarLink[col][d];
				dc = VarLinkDc[col][d];
				AddLLRVector(L_post[col], L_post[col], L_c2v[row][dc]);
			}
			DecodeOutput[col] = DecideLLRVector(L_post[col]);
		}

		// check whether satisfy the parity matrix
		bool decode_correct = true;
		for (int row = 0; row < ChkLen; row++)
		{
			int col, multip, syndrome;
			syndrome = 0;
			for (int d = 0; d < ChkDegree[row]; d++)
			{
				col = ChkLink[row][d];
				multip = GF.GFMultiply(ChkLinkGFe[row][d], DecodeOutput[col]);
				syndrome = GF.GFAdd(syndrome, multip);
			}
			if (syndrome != 0)
			{
				decode_correct = false;
				break;
			}
		}
		if (decode_correct)
		{
			return 1;
		}
		// message from var to check
		for (int col = 0; col < CodeLen; col++)
		{
			int row, dc;
			for (int dv = 0; dv < VarDegree[col]; dv++)
			{
				row = VarLink[col][dv];
				dc = VarLinkDc[col][dv];
				MinusLLRVector(L_v2c[col][dv], L_post[col], L_c2v[row][dc]);
			}
		}
		// message from check to var
		for (int row = 0; row < ChkLen; row++)
		{
			// sort to get the Nm maximum LLR
			int col, dv;
			for (int dc = 0; dc < ChkDegree[row]; dc++)
			{
				col = ChkLink[row][dc];
				dv = ChkLinkDv[row][dc];
				SortLLRVector(EMS_sort_L_v2c[dc], EMS_sort_Entr_v2c[dc], GFq, L_v2c[col][dv], GFq - 1);
			}
			// print
			/*
			ofstream fsort("sort.txt");
			for(int dc = 0; dc < ChkDegree[row]; dc ++)
			{
			for(int i = 0; i < GFq; i ++)
			{
			fsort << EMS_sort_L_v2c[dc][i] << " " << EMS_sort_Entr_v2c[dc][i] << " ";
			}
			fsort << endl;
			}
			fsort.close();
			*/
			for (int dc = 0; dc < ChkDegree[row]; dc++)
			{
				// reset the sum store vector to the munimum
				for (int q = 0; q < GFq; q++)
				{
					EMS_L_c2v[q] = -DBL_MAX;
				}
				// recursly exhaustly
				int sumNonele, diff;
				double sumNonLLR;
				// conf(q, 1)
				sumNonele = 0; sumNonLLR = 0; diff = 0;
				ConstructConf(GFq, 1, sumNonele, sumNonLLR, diff, 0, dc, ChkDegree[row] - 1, row);
				// conf(nm, nc)
				sumNonele = 0; sumNonLLR = 0; diff = 0;
				ConstructConf(EMS_Nm, EMS_Nc, sumNonele, sumNonLLR, diff, 0, dc, ChkDegree[row] - 1, row);
				// calculate each c2v LLR
				int v = 0;
				for (int k = 1; k < GFq; k++)
				{
					v = GF.GFMultiply(k, ChkLinkGFe[row][dc]);
					L_c2v[row][dc][k - 1] = (EMS_L_c2v[v] - EMS_L_c2v[0]) / EMS_Correction_Factor;
					if (L_c2v[row][dc][k - 1] < -1 * EMS_Correction_Offset)
					{
						L_c2v[row][dc][k - 1] = L_c2v[row][dc][k - 1] + EMS_Correction_Offset;
					}
					else if (L_c2v[row][dc][k - 1] >  EMS_Correction_Offset)
					{
						L_c2v[row][dc][k - 1] = L_c2v[row][dc][k - 1] - EMS_Correction_Offset;
					}
					else
					{
						L_c2v[row][dc][k - 1] = 0;
					}
				}
			}
		}
	}
	return 0;
}


int CNBLDPC::Decoding(double** L_ch, int* DecodeOutput)
{
	int decoderesult;
	switch (DecodeMethod)
	{
	case BP_DECODE:
		decoderesult = Decoding_BP(L_ch, DecodeOutput);
		break;
	case EMS_DECODE:
		decoderesult = Decoding_EMS(L_ch, DecodeOutput);
		break;
	case MinMax_DECODE:
		cout << "Min Max algorithm has not been developed" << endl;
		exit(-1);
		break;
	case T_EMS_DECODE:
		decoderesult = Decoding_TEMS(L_ch, DecodeOutput);
		break;
	case T_MinMax_DECODE:
		cout << "Trellis Min Max algorithm has not been developed" << endl;
		exit(-1);
		break;
	default:
		cout << "Decode Method is not specified!" << endl;
		exit(-1);
	}
	return decoderesult;
}


int CNBLDPC::SortLLRVector(double* L_sorted, int* Enre_sorted, int Len_sorted, double* L_sorting, int Len_sorting)
{
	int sortedCount = 0;
	L_sorted[0] = 0;
	Enre_sorted[0] = 0;
	double cand;
	while (sortedCount < Len_sorting)
	{
		// the candidate value
		cand = L_sorting[sortedCount];
		// put the cand on the tail
		L_sorted[sortedCount + 1] = cand;
		Enre_sorted[sortedCount + 1] = sortedCount + 1;
		// search forward
		for (int i = sortedCount; i >= 0; i--)
		{
			if (cand >= L_sorted[i])
			{
				L_sorted[i + 1] = L_sorted[i];
				L_sorted[i] = cand;
				Enre_sorted[i + 1] = Enre_sorted[i];
				Enre_sorted[i] = sortedCount + 1;
			}
			else
			{
				break;
			}
		}
		sortedCount++;
	}
	return 0;
}


int CNBLDPC::ConstructConf(int Nm, int Nc, int& sumNonele, double& sumNonLLR, int& diff, int begin, int except, int end, int row)
{
	if (begin > end)
	{
		if (sumNonLLR > EMS_L_c2v[sumNonele])
		{
			EMS_L_c2v[sumNonele] = sumNonLLR;
		}
	}
	else if (begin == except)
	{
		ConstructConf(Nm, Nc, sumNonele, sumNonLLR, diff, begin + 1, except, end, row);
		return 0;
	}
	else
	{
		for (int k = 0; k < Nm; k++)
		{
			sumNonele = GF.GFAdd(GF.GFMultiply(EMS_sort_Entr_v2c[begin][k], ChkLinkGFe[row][begin]), sumNonele);
			sumNonLLR = sumNonLLR + EMS_sort_L_v2c[begin][k];
			diff += (k != 0) ? 1 : 0;
			if (diff <= Nc)
			{
				ConstructConf(Nm, Nc, sumNonele, sumNonLLR, diff, begin + 1, except, end, row);
				sumNonele = GF.GFAdd(GF.GFMultiply(EMS_sort_Entr_v2c[begin][k], ChkLinkGFe[row][begin]), sumNonele);
				sumNonLLR = sumNonLLR - EMS_sort_L_v2c[begin][k];
				diff -= (k != 0) ? 1 : 0;
			}
			else
			{
				sumNonele = GF.GFAdd(GF.GFMultiply(EMS_sort_Entr_v2c[begin][k], ChkLinkGFe[row][begin]), sumNonele);
				sumNonLLR = sumNonLLR - EMS_sort_L_v2c[begin][k];
				diff -= (k != 0) ? 1 : 0;
				break;
			}
		}
	}
	return 0;
}


int CNBLDPC::Decoding_TEMS(double** L_ch, int* DecodeOutput)
{
	// initial
	for (int col = 0; col < CodeLen; col++)
	{
		for (int d = 0; d < VarDegree[col]; d++)
		{
			CopyLLRVector(L_v2c[col][d], L_ch[col]);
		}
	}
	for (int row = 0; row < ChkLen; row++)
	{
		for (int d = 0; d < ChkDegree[row]; d++)
		{
			ClearLLRVector(L_c2v[row][d]);
		}
	}
	// iterative decoding
	int iter_number = 0;
	while (iter_number++ < maxIter)
	{
		//		cout << iter_number << " ";
		// tentative decoding
		for (int col = 0; col < CodeLen; col++)
		{
			CopyLLRVector(L_post[col], L_ch[col]);
			int row, dc;
			for (int d = 0; d < VarDegree[col]; d++)
			{
				row = VarLink[col][d];
				dc = VarLinkDc[col][d];
				AddLLRVector(L_post[col], L_post[col], L_c2v[row][dc]);
			}
			DecodeOutput[col] = DecideLLRVector(L_post[col]);
		}

		// check whether satisfy the parity matrix
		bool decode_correct = true;
		for (int row = 0; row < ChkLen; row++)
		{
			int col, multip, syndrome;
			syndrome = 0;
			for (int d = 0; d < ChkDegree[row]; d++)
			{
				col = ChkLink[row][d];
				multip = GF.GFMultiply(ChkLinkGFe[row][d], DecodeOutput[col]);
				syndrome = GF.GFAdd(syndrome, multip);
			}
			if (syndrome != 0)
			{
				decode_correct = false;
				break;
			}
		}
		if (decode_correct)
		{
			return 1;
		}
		// message from var to check
		for (int col = 0; col < CodeLen; col++)
		{
			int row, dc;
			for (int dv = 0; dv < VarDegree[col]; dv++)
			{
				row = VarLink[col][dv];
				dc = VarLinkDc[col][dv];
				MinusLLRVector(L_v2c[col][dv], L_post[col], L_c2v[row][dc]);
			}
		}
		// message from check to var
		for (int row = 0; row < ChkLen; row++)
		{
			// get the most reliable symbol beta
			TEMS_Get_Beta(row);
			// transform to the delta domain
			TEMS_Get_deltaU(row);
			// get the two minimum value of each row and mark the Nr minimum pos
			TEMS_Get_Min(row);
			// mark the conf
			for (int q = 0; q < GFq; q++)
			{
				TEMS_deltaW[q] = DBL_MAX;
				TEMS_isSelected[q] = false;
			}
			int sumNonele, diff;
			double sumNonLLR;
			sumNonele = diff = 0; sumNonLLR = 0.0;
			TEMS_ConstructConf(TEMS_Nr, TEMS_Nc, sumNonele, sumNonLLR, diff, 0, ChkDegree[row] - 1);
			//calculate the output
			for (int dc = 0; dc < ChkDegree[row]; dc++)
			{
				// clear the temp vector TEMS_Lc2v
				for (int q = 0; q < GFq; q++)
				{
					TEMS_Lc2v[q] = DBL_MAX;
					TEMS_isUpdate[q] = false;
				}
				// calculate the dc output
				for (int eta = 0; eta < GFq; eta++)
				{
					int deviation = TEMS_Eta[eta][dc];
					int eta_plus_deviation = GF.GFAdd(eta, deviation);
					double minus = TEMS_deltaW[eta] - TEMS_deltaU[dc][deviation];
					if (TEMS_Lc2v[eta_plus_deviation] > minus)
					{
						TEMS_Lc2v[eta_plus_deviation] = minus;
						TEMS_isUpdate[eta_plus_deviation] = true;
					}
				}
				for (int q = 0; q < GFq; q++)
				{
					if (!TEMS_isUpdate[q])
					{
						TEMS_Lc2v[q] = (dc == TEMS_Min[q][0]) ?
							TEMS_deltaU[TEMS_Min[q][1]][q] : TEMS_deltaU[TEMS_Min[q][0]][q];
					}
				}
				// factor & offset
				/*
				for(int q = 0; q < GFq; q ++)
				{
				if(TEMS_Offset - TEMS_Lc2v[q] > 0)
				{
				TEMS_Lc2v[q] = 0;
				}
				else
				{
				TEMS_Lc2v[q] = (TEMS_Lc2v[q] - TEMS_Offset) / TEMS_Factor;
				}
				}
				*/
				// delta domain --> normal domain & repermutation
				int h_inverse = GF.GFInverse(ChkLinkGFe[row][dc]);
				int beta_syn = GF.GFAdd(TEMS_Syndrome, TEMS_Beta[dc]);
				double L0 = -1.0 * TEMS_Lc2v[beta_syn];
				for (int eta = 0; eta < GFq; eta++)
				{
					if (eta != beta_syn)
					{
						int beta = GF.GFMultiply(h_inverse, GF.GFAdd(eta, beta_syn));
						L_c2v[row][dc][beta - 1] = -1.0 * TEMS_Lc2v[eta] - L0;
						L_c2v[row][dc][beta - 1] = L_c2v[row][dc][beta - 1] / TEMS_Factor;
						if (L_c2v[row][dc][beta - 1] < -1 * TEMS_Offset)
						{
							L_c2v[row][dc][beta - 1] = L_c2v[row][dc][beta - 1] + TEMS_Offset;
						}
						else if (L_c2v[row][dc][beta - 1] >  TEMS_Offset)
						{
							L_c2v[row][dc][beta - 1] = L_c2v[row][dc][beta - 1] - TEMS_Offset;
						}
						else
						{
							L_c2v[row][dc][beta - 1] = 0;
						}
					}
				}
			}
		}
	}

	return 0;
}


int CNBLDPC::TEMS_Get_Beta(int row)
{
	TEMS_Syndrome = 0;
	for (int dc = 0; dc < ChkDegree[row]; dc++)
	{
		int col = ChkLink[row][dc];
		int dv = ChkLinkDv[row][dc];
		int h = ChkLinkGFe[row][dc];
		double max = 0;
		int max_ele = 0;
		for (int q = 1; q < GFq; q++)
		{
			if (L_v2c[col][dv][q - 1] > max)
			{
				max = L_v2c[col][dv][q - 1];
				max_ele = GF.GFMultiply(q, h);
			}
		}
		TEMS_Beta[dc] = max_ele;
		TEMS_Syndrome = GF.GFAdd(TEMS_Syndrome, max_ele);
	}
	return 0;
}


int CNBLDPC::TEMS_Get_deltaU(int row)
{
	for (int dc = 0; dc < ChkDegree[row]; dc++)
	{
		int col = ChkLink[row][dc];
		int dv = ChkLinkDv[row][dc];
		int h_inverse = GF.GFInverse(ChkLinkGFe[row][dc]);

		int beta_p = GF.GFMultiply(h_inverse, TEMS_Beta[dc]);
		double max;
		max = (beta_p == 0) ? 0 : L_v2c[col][dv][beta_p - 1];

		TEMS_deltaU[dc][TEMS_Beta[dc]] = max - 0;
		for (int x = 1; x < GFq; x++)
		{
			int eta = GF.GFAdd(x, TEMS_Beta[dc]);
			TEMS_deltaU[dc][eta] = max - L_v2c[col][dv][GF.GFMultiply(h_inverse, x) - 1];
		}
	}
	return 0;
}


int CNBLDPC::TEMS_Get_Min(int row)
{
	// sort
	for (int q = 0; q < GFq; q++)
	{
		// clear
		for (int dc = 0; dc < ChkDegree[row]; dc++)
		{
			TEMS_Min[q][dc] = dc;
		}
		// search min
		double min = TEMS_deltaU[0][q];
		int min_pos = 0;
		for (int dc = 1; dc < ChkDegree[row]; dc++)
		{
			for (int d_s = dc; d_s >= 1; d_s--)
			{
				double du = TEMS_deltaU[TEMS_Min[q][d_s]][q];
				double du_forward = TEMS_deltaU[TEMS_Min[q][d_s - 1]][q];
				if (du < du_forward)
				{
					int temp;
					temp = TEMS_Min[q][d_s];
					TEMS_Min[q][d_s] = TEMS_Min[q][d_s - 1];
					TEMS_Min[q][d_s - 1] = temp;
				}
				else
				{
					break;
				}
			}
		}

	}
	// set TEMS_inConf_Nr
	for (int dc = 0; dc < ChkDegree[row]; dc++)
	{
		TEMS_inConf_Nr[0][dc] = true;
	}
	for (int q = 1; q < GFq; q++)
	{
		for (int dc = 0; dc < ChkDegree[row]; dc++)
		{
			TEMS_inConf_Nr[q][dc] = false;
		}
		for (int index = 0; index < TEMS_Nr; index++)
		{
			if (index >= ChkDegree[row])
				break;
			int dc = TEMS_Min[q][index];
			TEMS_inConf_Nr[q][dc] = true;
		}
	}

	return 0;
}


int CNBLDPC::TEMS_ConstructConf(int Nr, int Nc, int& sumNonele, double& sumNonLLR, int& diff, int begin, int end)
{
	if (begin > end)
	{
		if (sumNonLLR < TEMS_deltaW[sumNonele])
		{
			TEMS_deltaW[sumNonele] = sumNonLLR;
			for (int d = 0; d <= end; d++)
			{
				TEMS_Eta[sumNonele][d] = TEMS_Eta_Candidate[d];
			}
		}
	}
	else
	{
		for (int q = 0; q < GFq; q++)
		{
			if (TEMS_inConf_Nr[q][begin] && !TEMS_isSelected[q])
			{
				diff += (q != 0) ? 1 : 0;
				if (diff <= Nc)
				{
					TEMS_isSelected[q] = (q != 0) ? true : false;
					sumNonele = GF.GFAdd(sumNonele, q);
					sumNonLLR += TEMS_deltaU[begin][q];
					TEMS_Eta_Candidate[begin] = q;
					TEMS_ConstructConf(Nr, Nc, sumNonele, sumNonLLR, diff, begin + 1, end);
					sumNonele = GF.GFAdd(sumNonele, q);
					sumNonLLR -= TEMS_deltaU[begin][q];
					TEMS_isSelected[q] = false;
					diff -= (q != 0) ? 1 : 0;
				}
				else
				{
					diff -= (q != 0) ? 1 : 0;
					break;
				}

				/*
				else
				{
				sumNonele = GF.GFAdd(sumNonele, q);
				sumNonLLR -= TEMS_deltaU[begin][q];
				diff -= (q != 0)? 1:0;
				break;
				}
				*/
			}
		}
	}
	return 0;
}
