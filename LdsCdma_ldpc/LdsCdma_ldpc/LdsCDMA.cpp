// LdsCDMA.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "Comm.h"

int main()
{
	CComm Comm;
	Comm.Initial();
	double EbNodb;
	int total, errBit, errFrm;
	int LDSerrBit;
	bool flag;
	double ber, fer;
	double LDSber;
	for (double snr = 0.5; snr < 5; snr += 0.2)
	{
		total = 0;
		errBit = 0;
		errFrm = 0;
		LDSerrBit = 0;
		flag = true;
		EbNodb = snr;
		Comm.SetEbN0(EbNodb);
		while (flag)
		{
			total += 1;
			Comm.Transmission();
			errBit += Comm.errBit;
			errFrm += Comm.errFrm;
			LDSerrBit += Comm.LDSerrBit;
			if (errFrm > 20&&total>20)
			{
				flag = false;
			}
		}
		LDSber = LDSerrBit / double(total) / 16.0 / 4000;
		ber = errBit / double(total) / 16.0 / 2000;
		fer = errFrm / double(total) / 16.0;
		cout << "EbNodb:" << EbNodb << "   BER:" << ber << "   FER:" << fer <<"   LDSber"<<LDSber<< "   errBit:" << errBit << "   errFrm:" << errFrm << "   total:" << total << endl;
	}
	while (1) {};
    return 0;
}

