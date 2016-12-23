#pragma once
class CNode
{
public:
	CNode();
	~CNode();
	unsigned int Degree;
//	int maxDegree;
	unsigned long *Index;
	unsigned long *EdgeNo;
};

