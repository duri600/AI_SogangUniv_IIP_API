#define _CRT_SECURE_NO_WARNINGS
#define MAKE_FILE	1		//1 : make

#include <stdio.h>
#include "ProcBuffers.h"
#include "VAD.h"
#include "RTIVA.h"
#include "SNMF.h"
#include "sigproc.h"
#include <io.h>
#include <time.h>
#include <iostream>
#include "header.h"

oldVAD *my_VAD;
realtime_IVA *rtIVA;
SNMF *m_snmf;
SNMF *f_snmf;
#if MAKE_FILE == 1
FILE **nmf, **iva;
double **out_buff;
short **IVA_out, **NMF_out;
#endif
double **Input;
double **basis;
double **basis_a, **basis_b, **basis_an, **basis_bn, **basis_n;
double **output, **output_temp, **input_temp;
double *test_buffer1;
double *test_buffer2;
double **xx_lp, **InitCond, *XX_LP, *XX, **x, *LPF;
double **vad_input_tmp;
double lpf[256] = { 2.08798576826208e-06, 2.29954344986799e-05, -8.35216069314825e-06, -1.01421599408546e-05, 1.75825071131748e-05, 6.68506278162495e-07, -2.37466843099458e-05, 1.65292262827601e-05, 1.86603000155727e-05, -3.47276454619664e-05, 1.39791928852911e-06, 4.27828841291473e-05,	-3.24679407522093e-05, -2.99725427486013e-05, 6.15260248626211e-05, -6.94391496384368e-06, -7.03040741114150e-05, 5.82408009886648e-05, 4.39850560826189e-05, -0.000101231819760140, 1.84188639093049e-05,
0.000108128664125297, -9.77414224812796e-05, -5.99649810253215e-05, 0.000157527792910078, -3.91650080197067e-05, -0.000157914699984545, 0.000155735031049291, 7.64078744251690e-05, -0.000234514782604306, 7.36122297632366e-05, 0.000221009284717929, -0.000237857827130511, -9.08819603229543e-05, 0.000336524435543129, -0.000127344522646467, -0.000298213623651684, 0.000350655612780098, 9.97613501424823e-05, -0.000468059888805500, 0.000207292215334781, 0.000389557076800700, -0.000501506101607508,
-9.80403878092228e-05, 0.000633560527987403, -0.000321828720343297, -0.000494095184236539, 0.000698643273967565, 7.90873763046442e-05, -0.000837285503953346, 0.000480844420057879, 0.000609619959544147, -0.000951134221010459, -3.44598790385002e-05, 0.00108314551209364, -0.000696005652956184, -0.000732455819102594, 0.00126890531102301, -4.64362037536257e-05, -0.00137463837825748, 0.000980942207594894, 0.000857217253013655, -0.00166293374478135, 0.000176679007271611, 0.00171482336592332, -0.00135171031885412,
-0.000976555872744570, 0.00214562278599450, -0.000372365974724178, -0.00210647594676115, 0.00182754320491227, 0.00108084917623704, -0.00273156581117595, 0.000653525908926587, 0.00255245356503915, -0.00243222943518740, -0.00115775411275296, 0.00343894766908327, -0.00104561120424174, -0.00305643069554936, 0.00319651199968515, 0.00119146660897295, -0.00429207443129812, 0.00158219963360274, 0.00362422367640721, -0.00416242631985649, -0.00116137980775673, 0.00532598247979084, -0.00230995159742169,
-0.00426620546541984, 0.00539142505709083, 0.00103945666930688, -0.00659509904844149, 0.00329821462490733, 0.00500185874499057, -0.00698024820349760, -0.000784864830481446, 0.00819058918837528, -0.00465870264975467, -0.00586892292500419, 0.00909441684396274, 0.000332170044372207, -0.0102780677249779, 0.00658974395002571, 0.00694399559445460, -0.0120466066251864, 0.000437482118368211, 0.0131904224264900, -0.00948907025662504, -0.00839595749341352, 0.0165095446428109, -0.00177826835482988,
-0.0176999110017040, 0.0143004638937376, 0.0106573513166513, -0.0242397828308403, 0.00437616097914038, 0.0260668422700029, -0.0239555296472238, -0.0151991935796053, 0.0416959813296214, -0.0110150817261730, -0.0488064516064210, 0.0544275764049411, 0.0319437385822245, -0.127073765350224, 0.0615707177784359, 0.517443329375340, 0.517443329375340, 0.0615707177784359, -0.127073765350224, 0.0319437385822245, 0.0544275764049411, -0.0488064516064210, -0.0110150817261730, 0.0416959813296214,
-0.0151991935796053, -0.0239555296472238, 0.0260668422700029, 0.00437616097914038, -0.0242397828308403, 0.0106573513166513, 0.0143004638937376, -0.0176999110017040, -0.00177826835482988, 0.0165095446428109, -0.00839595749341352, -0.00948907025662504, 0.0131904224264900, 0.000437482118368211, -0.0120466066251864, 0.00694399559445460, 0.00658974395002571, -0.0102780677249779, 0.000332170044372207, 0.00909441684396274, -0.00586892292500419, -0.00465870264975467, 0.00819058918837528,
-0.000784864830481446, -0.00698024820349760, 0.00500185874499057, 0.00329821462490733, -0.00659509904844149, 0.00103945666930688, 0.00539142505709083, -0.00426620546541984, -0.00230995159742169, 0.00532598247979084, -0.00116137980775673, -0.00416242631985649, 0.00362422367640721, 0.00158219963360274, -0.00429207443129812, 0.00119146660897295, 0.00319651199968515, -0.00305643069554936, -0.00104561120424174, 0.00343894766908327, -0.00115775411275296, -0.00243222943518740, 0.00255245356503915,
0.000653525908926587, -0.00273156581117595, 0.00108084917623704, 0.00182754320491227, -0.00210647594676115, -0.000372365974724178, 0.00214562278599450, -0.000976555872744570, -0.00135171031885412, 0.00171482336592332, 0.000176679007271611, -0.00166293374478135, 0.000857217253013655, 0.000980942207594894, -0.00137463837825748, -4.64362037536257e-05, 0.00126890531102301, -0.000732455819102594, -0.000696005652956184, 0.00108314551209364, -3.44598790385002e-05, -0.000951134221010459, 0.000609619959544147,
0.000480844420057879, -0.000837285503953346, 7.90873763046442e-05, 0.000698643273967565, -0.000494095184236539, -0.000321828720343297, 0.000633560527987403, -9.80403878092228e-05, -0.000501506101607508, 0.000389557076800700, 0.000207292215334781, -0.000468059888805500, 9.97613501424823e-05, 0.000350655612780098, -0.000298213623651684, -0.000127344522646467, 0.000336524435543129, -9.08819603229543e-05, -0.000237857827130511, 0.000221009284717929, 7.36122297632366e-05, -0.000234514782604306, 7.64078744251690e-05,
0.000155735031049291, -0.000157914699984545, -3.91650080197067e-05, 0.000157527792910078, -5.99649810253215e-05, -9.77414224812796e-05, 0.000108128664125297, 1.84188639093049e-05, -0.000101231819760140, 4.39850560826189e-05, 5.82408009886648e-05, -7.03040741114150e-05, -6.94391496384368e-06, 6.15260248626211e-05, -2.99725427486013e-05, -3.24679407522093e-05, 4.27828841291473e-05, 1.39791928852911e-06, -3.47276454619664e-05, 1.86603000155727e-05, 1.65292262827601e-05, -2.37466843099458e-05, 6.68506278162495e-07,
1.75825071131748e-05, -1.01421599408546e-05, -8.35216069314825e-06, 2.29954344986799e-05, 2.08798576826208e-06 };


ProcBuffers::ProcBuffers()
{
	int ch, i, j;
	my_VAD = new oldVAD(BufferSize, Nch, SamplingFreq);
	rtIVA = new realtime_IVA();
	m_snmf = new SNMF();
	f_snmf = new SNMF();
	int total_stack = tr_sec * SamplingFreq;

	Input = new double *[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		Input[ch] = new double[4 * BufferSize];
	}
	basis = new double*[size_basis];
	basis_a = new double*[size_basis];
	basis_b = new double*[size_basis];
	basis_an = new double*[size_basis];
	basis_bn = new double*[size_basis];
	for (i = 0; i < size_basis; i++)
	{
		basis[i] = new double[2 * R_x + 3 * R_d / 2];
		basis_a[i] = new double[R_x];
		basis_b[i] = new double[R_x];
		basis_an[i] = new double[R_d];
		basis_bn[i] = new double[R_d];
	}
	output = new double*[Nch];
	output_temp = new double*[Nch];
	input_temp = new double*[Nch];
#if MAKE_FILE == 1
	out_buff = new double*[Nch];
	IVA_out = new short*[Nch];
	NMF_out = new short*[Nch];
#endif
	for (ch = 0; ch < Nch; ch++)
	{
		output[ch] = new double[BufferSize];
		output_temp[ch] = new double[Nwin];
		input_temp[ch] = new double[Nwin];
#if MAKE_FILE == 1
		out_buff[ch] = new double[BufferSize];
		IVA_out[ch] = new short[frame_shift];
		NMF_out[ch] = new short[frame_shift];
#endif
	}
	test_buffer1 = new double[frame_shift];
	test_buffer2 = new double[frame_shift];

	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < Nwin; i++)
		{
			output_temp[ch][i] = 0.0;
		}
	}

#if MAKE_FILE == 1
	char file_name[2][500];
	nmf = new FILE*[Nch];
	iva = new FILE*[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		sprintf(file_name[0], ".\\output\\NMF_ch%d.pcm", ch + 1);
		nmf[ch] = fopen(file_name[0], "wb");
		sprintf(file_name[0], ".\\output\\IVA_ch%d.pcm", ch + 1);
		iva[ch] = fopen(file_name[0], "wb");
	}
#endif
	//down sampling & LPF
	InitCond = new double *[Nch];
	x = new double *[Nch];
	for (i = 0; i < Nch; i++)
	{
		InitCond[i] = new double[BufferSize];
		x[i] = new double[BufferSize];
		for (j = 0; j < BufferSize; j++)
		{
			InitCond[i][j] = 0;
			x[i][j] = 0;
		}
	}
	xx_lp = new double *[Nch];
	for (i = 0; i < Nch; i++)
	{
		xx_lp[i] = new double[3 * BufferSize];
		for (j = 0; j < 3 * BufferSize; j++)
			xx_lp[i][j] = 0;
	}

	XX = new double[BufferSize * 2 + 2];
	XX_LP = new double[BufferSize * 2 + 2];
	LPF = new double[BufferSize * 2 + 2];
	for (i = 0; i < BufferSize * 2 + 2; i++)
	{
		XX[i] = 0;
		XX_LP[i] = 0;
		LPF[i] = 0;
	}
	for (i = 0; i < BufferSize; i++)
		LPF[i] = lpf[i];
	hfft1_a(LPF, BufferSize * 2, 1);

	vad_input_tmp = zeros2(Nch, Nwin);
}

ProcBuffers::~ProcBuffers()
{
	int ch, i;
	delete my_VAD;
	delete rtIVA;
	delete m_snmf;
	delete f_snmf;

	for (ch = 0; ch < Nch; ch++)
	{
		delete[] Input[ch];
		delete[] output[ch];
		delete[] output_temp[ch];
		delete[] input_temp[ch];
#if MAKE_FILE == 1
		delete[] out_buff[ch];
		delete[] IVA_out[ch];
		delete[] NMF_out[ch];
#endif
	}
#if MAKE_FILE == 1
	delete[] IVA_out;
	delete[] NMF_out;
	delete[] out_buff;
#endif
	delete[] Input;
	delete[] output;
	delete[] output_temp;
	delete[] input_temp;

#if MAKE_FILE == 1
	char file_name[2][500];
	for (ch = 0; ch < Nch; ch++)
	{
		fclose(nmf[ch]);
		sprintf(file_name[0], ".\\output\\NMF_ch%d.pcm", ch + 1);
		sprintf(file_name[1], ".\\output\\NMF_ch%d.wav", ch + 1);
		pcm2wav(file_name[0], file_name[1], (long)(SamplingFreq));
		remove(file_name[0]);
		fclose(iva[ch]);
		sprintf(file_name[0], ".\\output\\IVA_ch%d.pcm", ch + 1);
		sprintf(file_name[1], ".\\output\\IVA_ch%d.wav", ch + 1);
		pcm2wav(file_name[0], file_name[1], (long)(SamplingFreq));
		remove(file_name[0]);
	}
#endif
	for (i = 0; i < size_basis; i++)
	{
		delete[] basis[i];
		delete[] basis_a[i];
		delete[] basis_b[i];
		delete[] basis_an[i];
		delete[] basis_bn[i];
	}
	delete[] basis;
	delete[] basis_a;
	delete[] basis_b;
	delete[] basis_an;
	delete[] basis_bn;
	delete[] test_buffer1;
	delete[] test_buffer2;
#if MAKE_FILE == 1
	delete[] nmf;
	delete[] iva;
#endif
	//Down Sampling & LPF
	for (i = 0; i < Nch; i++)
	{
		delete[] xx_lp[i];
		delete[] x[i];
		delete[] InitCond[i];
	}

	delete[] xx_lp;
	delete[] x;
	delete[] InitCond;
	delete[] LPF;
	delete[] XX_LP;
	delete[] XX;
	free2(vad_input_tmp);

}

void ProcBuffers::Process(double**input)
{
	FILE *nFile;
	static int Process_count = 0;
	static int VAD_count = 0;
	static int BuffCnt = 0, isNew16k = 0;
	int EPDcase, ch, i, j;
	int total_stack = tr_sec * SamplingFreq;

	if (Process_count == 0)
	{

		nFile = fopen("basis_mat.txt", "rt");
		if (nFile != NULL)
		{
			for (i = 0; i < size_basis; i++)
			{
				for (j = 0; j < 2 * R_x + 3 * R_d / 2; j++)
				{
					fscanf(nFile, "%lf", &basis[i][j]);
				}
			}
		}

		for (i = 0; i < fftlen / 2 + 1; i++)
		{
			for (j = 0; j < R_x; j++)
			{
				basis_a[i][j] = basis[i][j];
			}
			for (; j < 2 * R_x; j++)
			{
				basis_b[i][j - R_x] = basis[i][j];
			}
			for (; j < 2 * R_x + R_d; j++) //female noise = noise 50rank female_noise 50rank
			{
				basis_an[i][j - 2 * R_x] = basis[i][j];
			}
			for (j = 2 * R_x; j < 2 * R_x + R_d / 2; j++)//male noise
			{
				basis_bn[i][j - 2 * R_x] = basis[i][j];
			}
			for (j = 2 * R_x + R_d; j < 2 * R_x + 3 * R_d / 2; j++) //male noise
			{
				basis_bn[i][j - (2 * R_x + R_d / 2)] = basis[i][j];
			}

		}
	}
	for (ch = 0; ch < Nch; ch++) //NMF input
	{
		for (i = 0; i < 3 * BufferSize; i++)
		{
			input_temp[ch][i] = input_temp[ch][BufferSize + i];
		}
		for (i = 0; i < BufferSize; i++)
		{
			input_temp[ch][3 * BufferSize + i] = input[ch][i];
		}
	}
	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < Nwin; i++)
		{
			vad_input_tmp[ch + 1][i + 1] = input_temp[ch][i];
		}
	}

	Process_count++;
	if (Process_count >= 1 && Process_count <= 10)
	{
		my_VAD->Estimate_Noise(vad_input_tmp[1]);//10번째 프레임까지 노이즈 추측
	}
	else
	{
		EPDcase = my_VAD->Calculate_VAD(vad_input_tmp[1]);
		if (my_VAD->EPD_print() == 2 || my_VAD->EPD_print() == 1)
		{

			rtIVA->IVAprocessing(input_temp, output);

#if MAKE_FILE == 1
			for (i = 0; i < Nch; i++)
			{
				for (j = 0; j < BufferSize; j++)
				{
					out_buff[i][j] = output[i][j] * 32768.0;
					IVA_out[i][j] = (short)(out_buff[i][j]);
				}
				fwrite(IVA_out[i], sizeof(short), BufferSize, iva[i]);
			}
#endif
			for (ch = 0; ch < Nch; ch++) //NMF input
			{
				for (i = 0; i < 3 * BufferSize; i++)
				{
					output_temp[ch][i] = output_temp[ch][BufferSize + i];
				}
				for (i = 0; i < BufferSize; i++)
				{
					output_temp[ch][3 * BufferSize + i] = output[ch][i] * 32768.0;
				}
			}

			f_snmf->SNMF_test(output_temp[0], test_buffer1, basis_a, basis_bn);

			m_snmf->SNMF_test(output_temp[1], test_buffer2, basis_b, basis_an);

#if MAKE_FILE == 1
			for (i = 0; i < frame_shift; i++)
			{
				NMF_out[0][i] = (short)(test_buffer1[i]);
				NMF_out[1][i] = (short)(test_buffer2[i]);
			}
			fwrite(NMF_out[0], sizeof(short), frame_shift, nmf[0]);
			fwrite(NMF_out[1], sizeof(short), frame_shift, nmf[1]);
#endif
		}
	}
}

void ProcBuffers::Training(int type)
{
	int i, j, ch;
	double *a_stack, *b_stack, *n_stack;
	double **a_tr_out, **b_tr_out, **n_tr_out, **an_tr_out, **bn_tr_out;
	double **total_basis_x;
	int total_stack = tr_sec * SamplingFreq;

	a_stack = new double[total_stack];
	b_stack = new double[total_stack];
	n_stack = new double[total_stack];

	a_tr_out = new double*[size_basis];
	b_tr_out = new double*[size_basis];
	n_tr_out = new double*[size_basis];
	an_tr_out = new double*[size_basis];
	bn_tr_out = new double*[size_basis];

	total_basis_x = new double*[size_basis];
	for (i = 0; i < size_basis; i++)
	{
		a_tr_out[i] = new double[R_x];
		b_tr_out[i] = new double[R_x];
		n_tr_out[i] = new double[R_d / 2];
		an_tr_out[i] = new double[R_d / 2];
		bn_tr_out[i] = new double[R_d / 2];
		total_basis_x[i] = new double[2 * R_x + 3 * R_d / 2];
	}

	char mat_name[50] = "basis_mat_type";
	char type_char[5];
	
	sprintf(type_char, "%d", type);
	strcat(mat_name, type_char);
	strcat(mat_name, ".txt");

	if (_access(mat_name, 0) != 0)
	{
		m_snmf->basis_stack(total_stack, a_stack, b_stack, n_stack, type);

		m_snmf->SNMF_training(a_stack, total_stack, R_x, a_tr_out); //F03 or M01 data training

		std::cout << "a_tr_out fin!" << std::endl;

		m_snmf->SNMF_training(b_stack, total_stack, R_x, b_tr_out); //M02 data training

		std::cout << "b_tr_out fin!" << std::endl;

		m_snmf->SNMF_training(n_stack, total_stack, R_d / 2, n_tr_out);

		std::cout << "n_tr_out fin!" << std::endl;

		m_snmf->SNMF_training(b_stack, total_stack, R_d / 2, bn_tr_out);

		std::cout << "bn_tr_out fin!" << std::endl;

		m_snmf->SNMF_training(a_stack, total_stack, R_d / 2, an_tr_out);

		std::cout << "an_tr_out fin!" << std::endl;

		for (i = 0; i < size_basis; i++)
		{
			for (j = 0; j < R_x; j++)//F03 or M01_basis 100rank
			{
				total_basis_x[i][j] = a_tr_out[i][j];
			}
			for (; j < 2 * R_x; j++)//M02_basis 100rank
			{
				total_basis_x[i][j] = b_tr_out[i][j - R_x];
			}
			for (; j < 2 * R_x + R_d / 2; j++)//noise_basis 50rank
			{
				total_basis_x[i][j] = n_tr_out[i][j - 2 * R_x];
			}
			for (; j < 2 * R_x + R_d; j++)//female_noise_basis 50rank
			{
				total_basis_x[i][j] = an_tr_out[i][j - (2 * R_x + R_d / 2)];
			}
			for (; j < 2 * R_x + 3 * R_d / 2; j++)//male_noise_basis 50rank
			{
				total_basis_x[i][j] = bn_tr_out[i][j - (2 * R_x + R_d)];
			}
		}

		/*basis 저장하기*/
		FILE *pFile;
		pFile = fopen(mat_name, "wt");
		if (pFile != NULL)
		{
			for (i = 0; i < size_basis; i++)
			{
				for (j = 0; j < 2 * R_x + 3 * R_d / 2; j++)
				{
					fprintf(pFile, "%E\t", total_basis_x[i][j]);
				}
			}
		}
		fclose(pFile);
	}

	for (i = 0; i < size_basis; i++)
	{
		delete[] a_tr_out[i];
		delete[] b_tr_out[i];
		delete[] n_tr_out[i];
		delete[] an_tr_out[i];
		delete[] bn_tr_out[i];
		delete[] total_basis_x[i];
	}
	delete[] a_tr_out;
	delete[] b_tr_out;
	delete[] n_tr_out;
	delete[] an_tr_out;
	delete[] bn_tr_out;

	delete[] total_basis_x;

	delete[] a_stack;
	delete[] b_stack;
	delete[] n_stack;

}

