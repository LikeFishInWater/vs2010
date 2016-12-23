#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include "gf.h"
#include "Simulation.h"
using namespace std;

class CNBLDPC
{
public:
	CNBLDPC(void);
	~CNBLDPC(void);
	
	// Galois Field
	int GFq;
	CGF GF;
	
	// Code Parameter
	int CodeLen;
	int ChkLen;
	int* PuncturePosition;
	int PunctureLen;
	int maxVarDegree;
	int maxChkDegree;
	int* VarDegree;
	int* ChkDegree;
	int** VarLink;
	int** ChkLink;
	int** VarLinkGFe;
	int** ChkLinkGFe;
	int** VarLinkDc;
	int** ChkLinkDv;
//	bool Initial(int q, string NBFileName, int PunctureVarDegree, int maxIter, int randMsg, int deocdemethod, int ems_nm, int ems_nc, double ems_factor, double ems_offset);
	bool Initial(CSimulation &sim);

	int maxIter;
	
	// Encode
	int InitialEncode(void);
	int GaussEliminate(int** H, int* ExchangeColSource, int* ExchangeColDest, int& ExchangeColNum);
	int* EncodeExchangeColSource;
	int* EncodeExchangeColDest;
	int EncodeExchangeColNum;
	int ExchRow(int** H, int source, int destination);
	int AddRow(int** H, int source, int destination, int multifactor);
	int ExchCol(int** H, int source, int destination);
	int** EncodeVarLink;
	int** EncodeVarLinkGFe;
	int* EncodeVarDv;
	int Encode(int* msg_sym, int* code_sym);

	// BP
	double* L_sigma;
	double* L_rho;
	int L_Back(int row, int l);
	int LLR_BoxPlus(double* L, double* L1, double* L2, int A1, int A2);
	int L_Forward(int row, int l);
	int Decoding_BP(double** L_ch, int* DecodeOutput);
	double*** L_v2c;
	double*** L_c2v;
	double* LLRV;
	double** L_post;
	int CopyLLRVector(double* L_d, double* L_s);
	int ClearLLRVector(double* L);
	int AddLLRVector(double* L, double* L1, double* L2);
	int DecideLLRVector(double* LLR);
	int MinusLLRVector(double* L, double* L1, double* L2);

	// EMS
	int Decoding_EMS(double** L_ch, int* DecodeOutput);
	int EMS_Nm;
	int EMS_Nc;
	double TEMS_Factor;
	double TEMS_Offset;
	int DecodeMethod;
	int Decoding(double** L_ch, int* DecodeOutput);
	double** EMS_sort_L_v2c;
	int** EMS_sort_Entr_v2c;
	double* EMS_L_c2v;
//	double* EMS_L_c2v_temp;
	int SortLLRVector(double* L_sorted, int* Enre_sorted, int Len_sorted, double* L_sorting, int Len_sorting);
	int ConstructConf(int Nm, int Nc, int& sumNonele, double& sumNonLLR, int& diff, int begin, int except, int end, int row);
	double EMS_Correction_Factor;
	double EMS_Correction_Offset;

	// Trellis EMS
	int Decoding_TEMS(double** L_ch, int* DecodeOutput);
	int TEMS_Nr;
	int TEMS_Nc;
	double** TEMS_deltaU;
	int TEMS_Get_deltaU(int row);
	
	int* TEMS_Beta;
	int TEMS_Syndrome;
	int TEMS_Get_Beta(int row);
	
	int** TEMS_Min;
	bool** TEMS_inConf_Nr;
	int TEMS_Get_Min(int row);
	
	double* TEMS_deltaW;
	int** TEMS_Eta;
	int* TEMS_Eta_Candidate;
	bool* TEMS_isSelected;
	int TEMS_ConstructConf(int Nr, int Nc, int& sumNonele, double& sumNonLLR, int& diff, int begin, int end);
	
	double* TEMS_Lc2v;
//	bool* isUpdate;
	bool* TEMS_isUpdate;
	
};

