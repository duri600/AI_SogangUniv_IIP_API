#define _CRT_SECURE_NO_WARNINGS

#include "VAD.h"
#include "header.h"
#include "sigproc.h"
#include <float.h>


oldVAD::oldVAD(int BuffSize, int micnum, int SamFreq)
{
	Buffer_Size = BuffSize;
	Nol = Buffer_Size * 3 / 4;
	hanning_win = hanning(Buffer_Size);
	Nvar_mean = 0;;
	count = 0;
	CNThangup = 0;
	CNThangdown = 0;
	Mic_num = micnum;
	isSpeech = 0;
	Freq = SamFreq;
	input_short = new short[Buffer_Size];
	EPDstate = 0;
}

oldVAD::~oldVAD()
{
	free(hanning_win);
	delete[] input_short;
}

void oldVAD::Estimate_Noise(double *input)
{
	int i;
	double N_r_m, N_i_m, Nvar;


	N = stft((input), Buffer_Size, hanning_win, Nol);
	N_r_m = 0;
	N_i_m = 0;
	for (i = 1; i <= Buffer_Size / 2 + 1; i++)
	{
		N_r_m += N[i + i - 1][1];
		N_i_m += N[i + i][1];
	}
	N_r_m /= 257.0;
	N_i_m /= 257.0;
	Nvar = 0;
	for (i = 1; i <= Buffer_Size / 2 + 1; i++)
		Nvar += pow(N[i + i - 1][1] - N_r_m, 2) + pow(N[i + i][1] - N_i_m, 2);
	Nvar /= 257.0;

	if (count == 0)
		Nvar_mean = Nvar;
	else
		Nvar_mean = (Nvar_mean + Nvar) / 2;

	count++;
	free2(N);
	//printf("%f\n", Nvar_mean);

}

int oldVAD::Calculate_VAD(double *input)
{
	int i;
	double gamma[257], PSY_temp[257] = { 0, };
	double LRT = 0;


	PSY = 1;

	X = stft(input, Buffer_Size, hanning_win, Nol); // 입력값은 이것 역시 ANC 출력 신호
	for (i = 1; i <= Buffer_Size / 2 + 1; i++)
	{
		//2017.12.29===========Start==============================================
		gamma[i - 1] = (pow(X[i + i - 1][1], 2) + pow(X[i + i][1], 2)) / Nvar_mean;
		PSY_temp[i - 1] = gamma[i - 1] - log(gamma[i - 1]) - 1.0;
		PSY = PSY + PSY_temp[i - 1];
		//2017.12.29===========End================================================
	}
	free2(X);
	PSY = PSY / 257.0;
	if (PSY > VAD_THRESHOLD + 100)
		PSY = VAD_THRESHOLD + 100;



	EPDcase = 0;
	// EPD using VAD hangover scheme
	if (EPDstate == 0)
	{
		if (PSY > VAD_THRESHOLD)
		{
			EPDstate = 2;
			CNThangup = 1;
			EPDcase = 1;
		}
	}
	else if (EPDstate == 2)
	{
		if (PSY <= VAD_THRESHOLD)
		{
			if (CNThangup > FramesHangUp)
			{
				EPDstate = 1;
				CNThangdown = 1;
			}
			else
			{
				EPDstate = 0;
				EPDcase = 2;
			}
		}
		else
			CNThangup++;
	}
	else if (EPDstate == 1)
	{
		if (PSY <= VAD_THRESHOLD)
		{
			CNThangdown++;
			if (CNThangdown > FramesHangDn)
			{
				EPDstate = 0;
				EPDcase = 3;
			}
		}
		else
			EPDstate = 2;
	}



	//printf("%d\n", EPDstate);
	return EPDcase;

}

int oldVAD::EPD_print()
{
	return EPDstate;

}

int oldVAD::power_test(double *input)
{
	int i;
	double gamma[257], PSY_temp[257] = { 0, };
	double LRT = 0;
	int ok;

	PSY = 1;

	X = stft(input, Buffer_Size, hanning_win, Nol); // 입력값은 이것 역시 ANC 출력 신호
	for (i = 1; i <= Buffer_Size / 2 + 1; i++)
	{
		/*gamma[i - 1] = (pow(X[i + i - 1][1], 2) + pow(X[i + i][1], 2)) / Nvar_mean;
		PSY_temp[i - 1] = exp(gamma[i - 1] - log(gamma[i - 1]) - 1.0);
		PSY = PSY * PSY_temp[i - 1];*/

		gamma[i - 1] = (pow(X[i + i - 1][1], 2) + pow(X[i + i][1], 2)) / Nvar_mean;
		PSY_temp[i - 1] = gamma[i - 1] - log(gamma[i - 1]) - 1.0;
		PSY = PSY + PSY_temp[i - 1];
	}
	free2(X);
	//PSY = log(pow(PSY, 1.0 / 257.0));
	PSY = PSY / 257.0;
	if (PSY > VAD_THRESHOLD + 100)
		PSY = VAD_THRESHOLD + 100;
	ok = 0;
	if (PSY > VAD_THRESHOLD)
		ok = 1;

	return ok;

}