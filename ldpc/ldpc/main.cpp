#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <time.h>
#include "Channel.h"
#include "LDPC.h"
#include "Modulate.h"
#include "Simulate.h"
using namespace std;

int main()
{
	CSimulate simulate;
	CLDPC ldpc;
	CChannel channel;
	CModulate modulate;
	/*初始化LDPC H矩阵,信道channel,仿真参数*/
	simulate.Initial(ldpc, channel, modulate);
	/*开始仿真*/
	ofstream fout("Result.txt");
	fout << "  Code File: " << simulate.CodeFileName << "\n  MaxInteration: " << ldpc.maxIteration
		<< "\tStart: " << simulate.startSNR << "  Stop: " << simulate.stopSNR << "  Step: " << simulate.stepSNR << '\n'
		<< "  Modulation Type:  " << modulate.ModualtionType << "  Decoder Algorithm:  " << ldpc.DecodeAlgorithm << '\n' << endl;
	clock_t start, stop;
	int steps = int((simulate.stopSNR - simulate.startSNR) / simulate.stepSNR) + 1;
	cout << setw(5) << "SNR" << setw(15) << "BER" << setw(15) << "FER"
		<< setw(10) << "ErrB" << setw(10) << "ErrF" << setw(10) << "TestF" << endl;
	fout << setw(5) << "SNR" << setw(15) << "BER" << setw(15) << "FER"
		<< setw(10) << "ErrB" << setw(10) << "ErrF" << setw(10) << "TestF" << endl;
	for (int i = 0; i < steps; i++)
	{
		start = clock();
		/*清空上一个EbN0仿真数据*/
		simulate.errorFrames = simulate.errorBits = simulate.simulationCycles = 0;
		/*更新此次仿真EbN0*/
		simulate.EbN0 = (simulate.startSNR + i * simulate.stepSNR);
		simulate.sigma = sqrt(1.0 / (2 * ldpc.Rate * pow(2.0, modulate.ModualtionType) * pow(10.0, simulate.EbN0 / 10)));
		/*仿真1000次以上，或者错帧达到100帧*/
		while ((simulate.errorFrames < 30) && (simulate.simulationCycles < 1000))
		{
			/*产生随机消息序列*/
			simulate.GenMsgSeq();//产生PN序列
			/*LDPC编码*/
			ldpc.Encode(simulate.MsgSeq);//进行编码
			/*调制*/
			modulate.Modulation(ldpc.EncodingSeq);
			/*经过信道*/
			channel.AWGNChannel(modulate.ModSeq, simulate.sigma);
			/*解调*/
			modulate.Demodulation(channel.SymbolSeq, simulate.sigma);
			/*LDPC译码*/
			ldpc.Decode(modulate.DemodSeq);
			/*对比encoding,decoding序列，统计错误*/
			simulate.errorFrames += (ldpc.CalErrorBits() > 0) ? 1 : 0;
			simulate.errorBits += ldpc.CalErrorBits();
			simulate.simulationCycles++;
			simulate.BER = double(simulate.errorBits) / double(simulate.simulationCycles * ldpc.CodeLen);
			simulate.FER = double(simulate.errorFrames) / double(simulate.simulationCycles);
			if (simulate.simulationCycles % 10 == 0)
				cout << defaultfloat << setw(5) << simulate.EbN0 << setw(15) << scientific << simulate.BER << setw(15) << scientific << simulate.FER
				<< setw(10) << simulate.errorBits << setw(10) << simulate.errorFrames << setw(10) << simulate.simulationCycles << '\r';
		}
		stop = clock();
		/*最后一帧错误，只出现了一次，减回去*/
//		simulate.errorFrames -= (ldpc.CalErrorBits() > 0) ? 1 : 0;
//		simulate.errorBits -= ldpc.CalErrorBits();
//		simulate.simulationCycles --;
		/*计算BER，FER*/
		simulate.BER = double(simulate.errorBits) / double(simulate.simulationCycles * ldpc.CodeLen);
		simulate.FER = double(simulate.errorFrames) / double(simulate.simulationCycles);
		cout << defaultfloat << setw(5) << simulate.EbN0 << setw(15) << scientific << simulate.BER << setw(15) << scientific << simulate.FER
			<< setw(10) << simulate.errorBits << setw(10) << simulate.errorFrames << setw(10) << simulate.simulationCycles
			<< setw(10) << defaultfloat << double(stop - start) / CLOCKS_PER_SEC << " s\r" << endl;
		fout << defaultfloat << setw(5) << simulate.EbN0 << setw(15) << scientific << simulate.BER << setw(15) << scientific << simulate.FER
			<< setw(10) << simulate.errorBits << setw(10) << simulate.errorFrames << setw(10) << simulate.simulationCycles
			<< setw(10) << defaultfloat << double(stop - start) / CLOCKS_PER_SEC << " s\r" << endl;
	}
	cout << endl;
	fout << endl;
	fout.close();
	system("pause");
	return 0;
}