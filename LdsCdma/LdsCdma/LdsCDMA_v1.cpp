// LdsCDMA_v1.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "Comm.h"

int main()
{
	CComm Comm;
	Comm.Initial();
	bool flag;
	int total, errFrm, errBit, errTemp;
	double ber, fer;
	double EbNodb;
	for (double snr = 1; snr < 10; snr=snr+1)
	{
		flag = true;
		total = 0;
		errFrm = 0;
		errBit = 0;
		EbNodb = snr;
		Comm.SetEbN0(EbNodb);
		while (flag)
		{
			total = total + 1;		
			Comm.Transmission();
			if (Comm.errBit > 0)
			{
				errBit += Comm.errBit;
				errFrm += 1;
			}
			if (errFrm == 100)
				flag = false;
		}
		ber = errBit / 16.0 / total;
		fer = errFrm / double(total);
		cout <<"EbNodb:"<<EbNodb<< "   BER:" << ber << "   FER:" << fer <<"   errBit:"<<errBit<<"   errFrm:"<<errFrm<<"   total:"<<total<< endl;
	}
	while (1) {};
    return 0;
}

