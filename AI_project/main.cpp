#define _CRT_SECURE_NO_WARNINGS
#define _CRTDBG_MAP_ALLOC

#include <stdio.h>
#include "ProcBuffers.h"
#include "sigproc.h"
#include <iostream>
#include "SNMF.h"
#include <stdlib.h>
#include <crtdbg.h>
#include <time.h>
#include "header.h"

#pragma comment(lib,"Winmm.lib")

#ifdef _DEBUG 
#ifndef DBG_NEW 
#define DBG_NEW new (_NORMAL_BLOCK, __FILE__, __LINE__) 
#define new DBG_NEW 
#endif 
#endif  // _DEBUG

using namespace std;

void usage(void)
{
	// Error function in case of incorrect command-line
	// argument specifications
	std::cout << "\nuseage: AI_project S(N)_ch1.wav S(N)_ch2.wav Source_type \n";
	std::cout << "N is # of mixing signal.\n";
	std::cout << "Source_type : 0 -> Female & Male.\n";
	std::cout << "              1 -> Male & Male.\n";
	exit(0);
}


int main(int argc, char *argv[])
{
	if (argc<4)
	{
		usage();
	}
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

	//_CrtSetBreakAlloc(2129);
	clock_t start, end;
	start = clock();
	int ch, Nframe, i, j;
	ProcBuffers *proc;
	proc = new ProcBuffers();

	int type = (int)atoi(argv[3]);

	//souce & noise training
	proc->Training(type);

	char file_name[Nch][50];
	double *signal[Nch];

	//frame ¥‹¿ß input
	for (ch = 0; ch < Nch; ch++)
	{
		sprintf(file_name[ch], argv[ch + 1]);
		std::cout << "reading " << file_name[ch] << std::endl;
		signal[ch] = wavread(file_name[ch]);
	}

	int input_length = length(signal[1]);
	double **Input;
	Input = new double *[Nch];
	for (int i = 0; i < Nch; i++)
	{
		Input[i] = new double[BufferSize];
		for (int j = 0; j < BufferSize; j++)
		{
			Input[i][j] = 0;
		}
	}

	int iteration = input_length / BufferSize;
	for (Nframe = 0; Nframe <= iteration; Nframe++)
	{
		for (i = 0; i < Nch; i++)
		{
			for (j = 0; j < BufferSize; j++)
			{
				if (Nframe*BufferSize + j > input_length)
				{
					Input[i][j] = 0;
				}
				else
				{
					Input[i][j] = signal[i][Nframe*BufferSize + j + 1];
				}
			}

		}
		proc->Process(Input, type);

	}
	for (ch = 0; ch < Nch; ch++)
	{
		delete[] Input[ch];
	}
	delete[] Input;

	for (int i = 0; i < Nch; i++)
	{
		free(signal[i]);
	}
	end = clock();
	std::cout << "1frame : " << end - start << "ms" << std::endl;
	delete proc;
}
