#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include "Complex.h"
#include "LDPC.h"
using namespace std;

CNode::CNode()
	: Degree()
	//, maxDegree(0)
	//, Index(NULL)
	//, EdgeNo(NULL)
{
}


CNode::~CNode()
{
}

CLDPC::CLDPC()
	: CodeLen(0)
	, MsgLen(0)
	, ChkLen(0)
	, Rate(0)
	//, checkNode(NULL)
	//, varaiableNode(NULL)
	, maxVarDegree(0)
	, maxChkDegree(0)
	//, Q(NULL)
	//, R(NULL)
	//, EdgeSum(0)
	//, SourceSeq(NULL)
	, EncodingSeq(nullptr)
	, DecodingSeq(nullptr)
	, maxIteration(50)
	//, L(nullptr)
	, ExchangedinEncoding(nullptr)
	, encodingNode(NULL)
	, DecodeAlgorithm(0)
	, MS_alpha(0)
	, BF_varSeq(NULL)
	, BF_chkSeq(NULL)
{
}


CLDPC::~CLDPC()
{
	delete[] checkNode;
	delete[] varaiableNode;
	delete[] Q;
	delete[] R;
	//	delete[] L;
	delete[] SourceSeq;
	delete[] EncodingSeq;
	delete[] DecodingSeq;
	delete[] ExchangedinEncoding;
	delete[] encodingNode;
	if (DecodeAlgorithm == 2)
	{
		delete[] BF_chkSeq;
		delete[] BF_varSeq;
	}
}


void CLDPC::Encode(int * MsgSeq)
{
	//置消息位
	for (unsigned long i = 0; i < MsgLen; i++)
		EncodingSeq[i] = MsgSeq[i];
	//清空校验位
	for (unsigned long i = MsgLen; i < CodeLen; i++)
		EncodingSeq[i] = 0;
	//计算校验位
	for (unsigned long i = MsgLen; i < CodeLen; i++)
	{
		for (unsigned long d = 0; d < encodingNode[i - MsgLen].Degree; d++)
		{
			EncodingSeq[i] = (EncodingSeq[i] + EncodingSeq[encodingNode[i - MsgLen].Index[d]]) % 2;
		}
	}
	/*交换顺序，由于生成G矩阵中有一些列交换*/
	int tempbit;
	for (unsigned long row = 0; row < ChkLen; row++)
	{
		if (ExchangedinEncoding[row] != -1)
		{
			tempbit = EncodingSeq[MsgLen + row];
			EncodingSeq[MsgLen + row] = EncodingSeq[ExchangedinEncoding[row]];
			EncodingSeq[ExchangedinEncoding[row]] = tempbit;
		}
	}
}


unsigned long CLDPC::CalErrorBits()
{
	unsigned long errorbits = 0;
	for (unsigned long i = 0; i < CodeLen; i++)
		errorbits += (EncodingSeq[i] == DecodingSeq[i]) ? 0 : 1;
	return errorbits;
}


void CLDPC::Initial(string CodeFileName)
{
	ifstream fin(CodeFileName);
	if (!fin.is_open())
	{
		cerr << "Cannot open " << CodeFileName << endl;
		exit(0);
	}
	//读入LDPC码长，校验节点个数，消息节点个数，码率
	fin >> CodeLen >> ChkLen;
	MsgLen = CodeLen - ChkLen;
	Rate = double(MsgLen) / double(CodeLen);
	fin >> maxVarDegree >> maxChkDegree;
	//分配未编码序列，编码后序列，解码序列
	SourceSeq = new int[MsgLen]();
	EncodingSeq = new int[CodeLen]();
	DecodingSeq = new int[CodeLen]();
	if (DecodeAlgorithm == 2)
	{
		BF_chkSeq = new unsigned long[ChkLen]();
		BF_varSeq = new unsigned long[CodeLen]();
	}
	//	L = new double[CodeLen]();
	//分配变量节点和校验节点空间
	varaiableNode = new CNode[CodeLen];
	checkNode = new CNode[ChkLen];
	//读入各个节点度
	for (unsigned long i = 0; i < CodeLen; i++)
		fin >> varaiableNode[i].Degree;
	for (unsigned long i = 0; i < ChkLen; i++)
		fin >> checkNode[i].Degree;
	//统计边的总数
	EdgeSum = 0;
	for (unsigned long i = 0; i < CodeLen; i++)
		EdgeSum = EdgeSum + varaiableNode[i].Degree;
	//分配边变量的内存空间
	Q = new double[EdgeSum]();
	R = new double[EdgeSum]();
	//读入与各个节点相连的节点号
	unsigned long tempIndex;
	for (unsigned long i = 0; i < CodeLen; i++)
	{
		varaiableNode[i].Index = new unsigned long[varaiableNode[i].Degree];
		varaiableNode[i].EdgeNo = new unsigned long[varaiableNode[i].Degree];
		for (unsigned int j = 0; j < varaiableNode[i].Degree; j++)
		{
			fin >> tempIndex;
			if (tempIndex > 0)
				varaiableNode[i].Index[j] = tempIndex - 1;
		}
	}
	for (unsigned long i = 0; i < ChkLen; i++)
	{
		checkNode[i].Index = new unsigned long[checkNode[i].Degree];
		checkNode[i].EdgeNo = new unsigned long[checkNode[i].Degree];
		for (unsigned int j = 0; j < checkNode[i].Degree; j++)
		{
			fin >> tempIndex;
			if (tempIndex > 0)
				checkNode[i].Index[j] = tempIndex - 1;
		}
	}
	//计算两两节点间的边的编号
	unsigned long edgeNum = 0;
	for (unsigned long i = 0; i < ChkLen; i++)
	{
		for (unsigned int j = 0; j < checkNode[i].Degree; j++)
			checkNode[i].EdgeNo[j] = edgeNum++;
	}
	for (unsigned long i = 0; i < ChkLen; i++)
	{
		for (unsigned int j = 0; j < checkNode[i].Degree; j++)
			for (unsigned int k = 0; k < varaiableNode[checkNode[i].Index[j]].Degree; k++)
			{
				if (varaiableNode[checkNode[i].Index[j]].Index[k] == i)
					varaiableNode[checkNode[i].Index[j]].EdgeNo[k] = checkNode[i].EdgeNo[j];
			}
	}
	fin.close();
}


void CLDPC::ConstructGenMatrix()
{
	/*构建一个临时矩阵tempH*/
	int * tempH = new int[ChkLen * CodeLen]();
	ExchangedinEncoding = new long[ChkLen]();
	long col, row;
	for (row = 0; row < ChkLen; row++)
	{
		ExchangedinEncoding[row] = -1;
		for (unsigned int d = 0; d < checkNode[row].Degree; d++)
		{
			col = checkNode[row].Index[d];
			tempH[col + row * CodeLen] = 1;
		}
	}
	/*对tempH进行高斯消去，形成tempH=[P|I]*/
	long seekrow, seekcol;
	bool oneisfound = false;
	bool isFullRank = true;
	for (row = ChkLen - 1; row >= 0; row--)
	{
		col = MsgLen + row;
		//如果这一行最后一个元素不是1，那么向上找，直到为1的seek行，然后交换两行seek和row行
		if (tempH[col + row * CodeLen] != 1)
		{
			oneisfound = false;
			for (seekrow = row; seekrow >= 0; seekrow--)
			{
				if (tempH[col + seekrow * CodeLen] == 1)
				{
					//交换找到的第一个为1的行seekrow，交换seekrow和当前行row
					int tempexchange;
					for (unsigned long i = 0; i < CodeLen; i++)
					{
						tempexchange = tempH[i + row * CodeLen];
						tempH[i + row * CodeLen] = tempH[i + seekrow * CodeLen];
						tempH[i + seekrow * CodeLen] = tempexchange;
					}
					oneisfound = true;
					break;
				}
			}
			//最后几行会出现不满秩的情况，待完成
			//例如第2行出现不满秩，从这一行这一列往左数，如果发现有1，就交换两列
			if (!oneisfound)
			{
				for (seekcol = col - 1; seekcol >= 0; seekcol--)
				{
					if (tempH[seekcol + row * CodeLen] == 1)
					{
						int tempexchange;
						for (unsigned long int i = 0; i < ChkLen; i++)
						{
							tempexchange = tempH[seekcol + i * CodeLen];
							tempH[seekcol + i * CodeLen] = tempH[col + i*CodeLen];
							tempH[col + i*CodeLen] = tempexchange;
						}
						/*需要记录下交换了哪些列，在编码结束后，需要将这些列给变回来*/
						ExchangedinEncoding[row] = seekcol;
						break;//找到就退出寻找
						oneisfound = true;
					}
				}
			}
			/*			if (!oneisfound)
			{
			isFullRank = false;
			cout << "Error Not Full Rank Gauss elimate at row = " << row << endl;
			exit(0);
			}
			*/
			//			fout << row << '\n';
		}

		//从此行向上找，直到为1的第seek行，seek行加上row行，使得当前col列row行以上没有1
		for (seekrow = row - 1; seekrow >= 0; seekrow--)
		{
			if (tempH[col + seekrow * CodeLen] == 1)
			{
				//交换找到为1的行seekrow，seekrow和当前行row相加
				for (unsigned long i = 0; i < CodeLen; i++)
					tempH[i + seekrow * CodeLen] = (tempH[i + seekrow * CodeLen] + tempH[i + row * CodeLen]) % 2;
			}
		}
	}
	/*
	ofstream fout("gauss_elimate.txt");
	for (row = 0; row < ChkLen; row++)
	{
	for (col = 0; col < CodeLen; col++)
	{
	fout << int(tempH[col + row * CodeLen]) << " ";
	}
	fout << endl;
	}
	fout.close();
	*/
	encodingNode = new CNode[ChkLen]();
	unsigned long degreesum = 0;
	for (row = 0; row < ChkLen; row++)
	{
		degreesum = 0;
		for (col = 0; col < MsgLen + row; col++)
		{
			if (tempH[col + row*CodeLen] == 1)
				degreesum++;;
		}
		encodingNode[row].Index = new unsigned long[degreesum]();
	}
	for (row = 0; row < ChkLen; row++)
	{
		for (col = 0; col < MsgLen + row; col++)
		{
			if (tempH[col + row*CodeLen] == 1)
			{
				encodingNode[row].Index[encodingNode[row].Degree] = col;
				encodingNode[row].Degree++;
			}
		}
	}
	delete[] tempH;
}


bool CLDPC::EncodeCheck()
{
	int chk;
	bool EncodingCorrect = true;
	for (unsigned long i = 0; i < ChkLen; i++)
	{
		chk = 0;
		for (int d = 0; d < checkNode[i].Degree; d++)
			chk = chk + EncodingSeq[checkNode[i].Index[d]];
		if (chk % 2 != 0)
		{
			EncodingCorrect = false;
			break;
		}
	}
	return EncodingCorrect;
}


void CLDPC::Decode(double * LLR)
{
	switch (DecodeAlgorithm)
	{
	case 0:
		Decode_BP(LLR);
		break;
	case 1:
		Decode_MS(LLR, MS_alpha);
		break;
	case 2:
		Decode_BF(LLR);
		break;
	default:
		break;
	}
}


void CLDPC::Decode_BP(double * LLR)
{
	/*BP LLR*/
	/*计算外信息，对数似然比*/
	/*初始化Q，Q为变量节点的输出消息*/
	unsigned long i;
	unsigned int d;
	for (i = 0; i < CodeLen; i++)
	{
		//d为第i个变量节点的第d条边,该条边在整个图中的编号为varaiableNode[i].EdgeNo[d]
		for (d = 0; d < varaiableNode[i].Degree; d++)
		{
			Q[varaiableNode[i].EdgeNo[d]] = LLR[i];
		}
	}
	int iterationTimes = 0;
	unsigned long DecodingOver = 0;
	while (1)
	{
		/*更新R，R为校验节点的输出消息*/
		for (i = 0; i < ChkLen; i++)
		{
			for (d = 0; d < checkNode[i].Degree; d++)
			{
				bool firstedge = true;
				double L1, L2;
				double sign1, sign2, minmum;
				double minussign, plussign;
				for (unsigned int dseek = 0; dseek < checkNode[i].Degree; dseek++)
				{
					if (dseek != d)
					{
						if (firstedge)
						{
							L2 = Q[checkNode[i].EdgeNo[dseek]];
							firstedge = false;
						}
						else
						{
							L1 = Q[checkNode[i].EdgeNo[dseek]];
							sign1 = (L1 < 0) ? -1 : 1;
							sign2 = (L2 < 0) ? -1 : 1;
							minmum = ((L1 * sign1) >(L2 * sign2)) ? (L2 * sign2) : (L1 * sign1);
							plussign = (L1 + L2 < 0) ? -1 : 1;
							minussign = (L1 - L2 < 0) ? -1 : 1;
							L2 = sign1 * sign2 * minmum + log(1 + exp(-1 * plussign * (L1 + L2))) - log(1 + exp(-1 * minussign * (L1 - L2)));
						}
					}
				}
				R[checkNode[i].EdgeNo[d]] = L2;
			}
		}
		/*
		double sign, sum, temp;
		for (i = 0; i < ChkLen; i++)
		{
		sum = 0;
		sign = 1;
		for (d = 0; d < checkNode[i].Degree; d++)
		{
		temp = Q[checkNode[i].EdgeNo[d]];
		if (temp > 0)
		{
		sum += log((exp(temp) + 1) / (exp(temp) - 1));
		sign *= 1;
		}
		else
		{
		sum += log((exp(0 - temp) + 1) / (exp(0 - temp) - 1));
		sign *= -1;
		}
		}
		for (d = 0; d < checkNode[i].Degree; d++)
		{
		temp = Q[checkNode[i].EdgeNo[d]];
		if (temp > 0)
		{
		temp = sum - log((exp(temp) + 1) / (exp(temp) - 1));
		R[checkNode[i].EdgeNo[d]] = sign * log((exp(temp) + 1) / (exp(temp) - 1));
		}
		else
		{
		temp = sum - log((exp(0 - temp) + 1) / (exp(0 - temp) - 1));
		R[checkNode[i].EdgeNo[d]] = sign * log((exp(temp) - 1) / (exp(temp) + 1));
		}
		}
		}
		*/
		/*更新Q，Q为变量节点的输出消息*/
		/*判决*/
		DecodingOver = 0;
		double log_sum;
		for (i = 0; i < CodeLen; i++)
		{
			//varaiableNode[i]所有输入的R值求和
			//顺便计算把各变量节点判决值
			log_sum = LLR[i];
			for (d = 0; d < varaiableNode[i].Degree; d++)
				log_sum = log_sum + R[varaiableNode[i].EdgeNo[d]];
			DecodingSeq[i] = (log_sum < 0) ? 1 : 0;
			//求varaiableNode[i]各输出边的Q
			for (d = 0; d < varaiableNode[i].Degree; d++)
				Q[varaiableNode[i].EdgeNo[d]] = log_sum - R[varaiableNode[i].EdgeNo[d]];
		}
		/*并校验解码结果是否满足H*/
		int chk;
		for (i = 0; i < ChkLen; i++)
		{
			chk = 0;
			for (d = 0; d < checkNode[i].Degree; d++)
				chk = chk + DecodingSeq[checkNode[i].Index[d]];
			DecodingOver += chk % 2;
		}
		/*如果达到最大迭代次数，或者，译码成功，结束迭代*/
		iterationTimes = iterationTimes + 1;
		if (DecodingOver == 0 || iterationTimes >= maxIteration)
			break;
	}
}


void CLDPC::Decode_MS(double * LLR, double alpha)
{
	/*MS LLR*/
	/*初始化Q，Q为变量节点的输出消息*/
	unsigned long i;
	unsigned int d;
	for (i = 0; i < CodeLen; i++)
	{
		//d为第i个变量节点的第d条边,该条边在整个图中的编号为varaiableNode[i].EdgeNo[d]
		for (d = 0; d < varaiableNode[i].Degree; d++)
		{
			Q[varaiableNode[i].EdgeNo[d]] = LLR[i];
		}
	}
	int iterationTimes = 0;
	unsigned long DecodingOver = 0;
	while (1)
	{
		/*更新R，R为校验节点的输出消息*/
		double sign, min1, min0, tempq, tempsign, tempr;//min1倒数第二小，min0最小
		for (i = 0; i < ChkLen; i++)
		{
			sign = ((Q[checkNode[i].EdgeNo[0]] > 0) ? 1 : -1) * ((Q[checkNode[i].EdgeNo[1]] > 0) ? 1 : -1);
			min1 = (Q[checkNode[i].EdgeNo[1]] > 0) ? Q[checkNode[i].EdgeNo[1]] : 0 - Q[checkNode[i].EdgeNo[1]];
			min0 = (Q[checkNode[i].EdgeNo[0]] > 0) ? Q[checkNode[i].EdgeNo[0]] : 0 - Q[checkNode[i].EdgeNo[0]];
			if (min1 < min0)
			{
				double temp;
				temp = min0;
				min0 = min1;
				min1 = temp;
			}
			//cout << Q[checkNode[i].EdgeNo[0]] << '\t' << min1 << '\t' << min0 << '\t' << sign << endl;
			//找出最小值min0和次小值min1
			for (d = 2; d < checkNode[i].Degree; d++)
			{
				sign *= (Q[checkNode[i].EdgeNo[d]] > 0) ? 1 : -1;
				tempq = (Q[checkNode[i].EdgeNo[d]] > 0) ? Q[checkNode[i].EdgeNo[d]] : 0 - Q[checkNode[i].EdgeNo[d]];
				if (tempq < min1)
				{
					if (tempq < min0)
					{
						min1 = min0;
						min0 = tempq;
					}
					else
					{
						min1 = tempq;
					}
				}
				//cout << Q[checkNode[i].EdgeNo[d]]  << '\t' << min1 << '\t' << min0 << '\t' << sign <<  endl;
			}
			//如果本身就是最小值Min0，那么就是alpha * tempsign * min1，否则，就是alpha * tempsign * min0
			for (d = 0; d < checkNode[i].Degree; d++)
			{
				tempsign = sign * ((Q[checkNode[i].EdgeNo[d]] > 0) ? 1 : -1);
				tempr = (Q[checkNode[i].EdgeNo[d]] > 0) ? Q[checkNode[i].EdgeNo[d]] : 0 - Q[checkNode[i].EdgeNo[d]];
				if (tempr != min0)
				{
					R[checkNode[i].EdgeNo[d]] = alpha * tempsign * min0;
				}
				else
				{
					R[checkNode[i].EdgeNo[d]] = alpha * tempsign * min1;
				}
			}
		}
		/*更新Q，Q为变量节点的输出消息*/
		/*判决*/
		DecodingOver = 0;
		double log_sum;
		for (i = 0; i < CodeLen; i++)
		{
			//varaiableNode[i]所有输入的R值求和
			//顺便计算把各变量节点判决值
			log_sum = LLR[i];
			for (d = 0; d < varaiableNode[i].Degree; d++)
				log_sum = log_sum + R[varaiableNode[i].EdgeNo[d]];
			DecodingSeq[i] = (log_sum < 0) ? 1 : 0;
			//求varaiableNode[i]各输出边的Q
			for (d = 0; d < varaiableNode[i].Degree; d++)
				Q[varaiableNode[i].EdgeNo[d]] = log_sum - R[varaiableNode[i].EdgeNo[d]];
		}
		/*并校验解码结果是否满足H*/
		int chk;
		for (i = 0; i < ChkLen; i++)
		{
			chk = 0;
			for (d = 0; d < checkNode[i].Degree; d++)
				chk = chk + DecodingSeq[checkNode[i].Index[d]];
			DecodingOver += chk % 2;
		}
		/*如果达到最大迭代次数，或者，译码成功，结束迭代*/
		iterationTimes = iterationTimes + 1;
		if (DecodingOver == 0 || iterationTimes >= maxIteration)
			break;
	}
}


void CLDPC::Decode_BF(double * LLR)
{
	/*BF LLR*/
	unsigned long i;
	unsigned int d;
	for (i = 0; i < CodeLen; i++)
	{
		DecodingSeq[i] = (LLR[i] > 0) ? 0 : 1;
	}
	int iterationTimes = 0;
	unsigned long DecodingOver = 0;
	unsigned long maxIndex;
	unsigned int max;
	while (1)
	{
		//校验
		DecodingOver = 0;
		for (i = 0; i < ChkLen; i++)
		{
			BF_chkSeq[i] = 0;
			for (d = 0; d < checkNode[i].Degree; d++)
				BF_chkSeq[i] = (BF_chkSeq[i] + DecodingSeq[checkNode[i].Index[d]]) % 2;
			DecodingOver += BF_chkSeq[i];
		}
		//计算每个变量节点使得几个校验节点不为0
		for (i = 0; i < CodeLen; i++)
		{
			BF_varSeq[i] = 0;
			for (d = 0; d < varaiableNode[i].Degree; d++)
			{
				if (BF_chkSeq[varaiableNode[i].Index[d]] != 0)
				{
					BF_varSeq[i] += 1;
				}
			}
		}
		if (DecodingOver != 0)
		{
			//寻找和最大的变量节点
			max = BF_varSeq[0];
			for (i = 0; i < CodeLen; i++)
			{
				if (BF_varSeq[i] > max)
				{
					max = BF_varSeq[i];
				}
			}
			//翻转和最大的变量节点
			for (i = 0; i < CodeLen; i++)
			{
				if (BF_varSeq[i] == max)
				{
					DecodingSeq[i] = 1 - DecodingSeq[i];
				}
			}
			iterationTimes = iterationTimes + 1;
		}
		if (DecodingOver == 0 || iterationTimes >= maxIteration)
		{
			break;
		}
	}
}
