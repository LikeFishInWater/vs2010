// LdsCDMA_v2.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "Comm.h"

int main()
{
	CComm Comm;
	Comm.Initial();
	double EbNodb;
	int total, errBit, errFrm;
	bool flag;
	double ber, fer;
	for (double snr = 1; snr < 8; snr += 0.5)
	{
		total = 0;
		errBit = 0;
		errFrm = 0;
		flag = true;
		EbNodb = snr;
		Comm.SetEbN0(EbNodb);
		while (flag)
		{
			total += 1;
			Comm.Transmission();
			Comm.err();
			errBit += Comm.errBit;
			errFrm += Comm.errFrm;
			if (errFrm > 30 || total>16 * 1000)
			{
				flag = false;
			}
		}
		ber = errBit / 16.0 / 2000 /total;
		fer = errFrm / 16.0 / total;
		cout << "EbNodb:" << EbNodb << "   BER:" << ber << "   FER:" << fer << "   errBit:" << errBit << "   errFrm:" << errFrm << "   total:" << total << endl;
	}
	while (1) {};
    return 0;
}

