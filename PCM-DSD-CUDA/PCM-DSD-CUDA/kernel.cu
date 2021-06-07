#include <fstream>
#include <string>
#include <omp.h>
#include <windows.h>
#include <iostream>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cufft.h>

using namespace std;
void Check(cudaError_t status)
{

	if (status != cudaSuccess)
	{
		cout << "行号:" << __LINE__ << endl;
		cout << "错误:" << cudaGetErrorString(status) << endl;
	}
}
template<class T>
T reverse_endian(T value)
{
	char* first = reinterpret_cast<char*>(&value);
	char* last = first + sizeof(T);
	std::reverse(first, last);
	return value;
}

bool RequireWriteData(const wchar_t* filepath, const wchar_t* flag, const wchar_t* FileMode, FILE** WriteData) {
	wchar_t DriveName[3];
	wchar_t DirectoryName[256];
	wchar_t FileName[256];
	_wsplitpath_s(filepath, DriveName, 3, DirectoryName, 256, FileName, 256, NULL, 0);
	wchar_t WritePath[260];
	wcscpy_s(WritePath, 260, DriveName);
	wcscat_s(WritePath, 260, DirectoryName);
	wcscat_s(WritePath, 260, FileName);
	wcscat_s(WritePath, 260, flag);

	errno_t error;

	if ((error = _wfopen_s(WriteData, WritePath, FileMode)) != 0) {
		return false;
	}
	return true;
}

bool TrushFile(const wchar_t* filepath, const wchar_t* flag) {
	wchar_t DriveName[3];
	wchar_t DirectoryName[256];
	wchar_t FileName[256];
	_wsplitpath_s(filepath, DriveName, 3, DirectoryName, 256, FileName, 256, NULL, 0);
	wchar_t DeletePath[260];
	wcscpy_s(DeletePath, 260, DriveName);
	wcscat_s(DeletePath, 260, DirectoryName);
	wcscat_s(DeletePath, 260, FileName);
	wcscat_s(DeletePath, 260, flag);

	if (!DeleteFile(DeletePath)) {
		return false;
	}
	return true;
}

bool WAV_Metadata(const wchar_t* filepath, unsigned __int32& samplingrate, unsigned short& bitdepth, unsigned __int32& samplesize)
{
	FILE* fprwav;
	errno_t error;
	if ((error = _wfopen_s(&fprwav, filepath, L"rb")) != 0) {
		return false;
	}

	unsigned short fmtID;
	unsigned short chnum;

	_fseeki64(fprwav, 20, SEEK_CUR);
	fread(&fmtID, 2, 1, fprwav);
	if (fmtID != 1) {
		fclose(fprwav);
		return false;
	}

	fread(&chnum, 2, 1, fprwav);
	if (chnum != 2) {
		fclose(fprwav);
		return false;
	}

	fread(&samplingrate, 4, 1, fprwav);
	if (samplingrate == 44100 || samplingrate == 44100 * 2 || samplingrate == 44100 * 4 || samplingrate == 44100 * 8 || samplingrate == 44100 * 16 ||
		samplingrate == 48000 || samplingrate == 48000 * 2 || samplingrate == 48000 * 4 || samplingrate == 48000 * 8 || samplingrate == 48000 * 16) {
	}
	else {
		fclose(fprwav);
		return false;
	}

	_fseeki64(fprwav, 6, SEEK_CUR);
	fread(&bitdepth, 2, 1, fprwav);
	if (bitdepth == 16 || bitdepth == 24 || bitdepth == 32 || bitdepth == 64) {
	}
	else {
		fclose(fprwav);
		return false;
	}
	_fseeki64(fprwav, 4, SEEK_CUR);
	fread(&samplesize, 4, 1, fprwav);

	fclose(fprwav);
	return true;
}

bool TmpWriteData(const wchar_t* filepath, unsigned short bitdepth, unsigned __int32 samplesize, FILE* tmpl, FILE* tmpr, int Times) {
	FILE* wavread;
	errno_t error;

	if ((error = _wfopen_s(&wavread, filepath, L"rb")) != 0) {
		return false;
	}

	_fseeki64(wavread, 44, SEEK_CUR);

	__int64 buffer_int = 0;
	double buffer_double = 0;
	double bit = pow(2, bitdepth - 1);
	unsigned __int64 writelength = samplesize / (bitdepth / 4);
	__int64 fillsize = 4096 * Times - writelength % (4096 * Times);
	for (int i = 0; i < fillsize; i++) {
		fwrite(&buffer_double, 8, 1, tmpl);
		fwrite(&buffer_double, 8, 1, tmpr);
	}


	for (int i = 0; i < writelength; i++) {
		fread(&buffer_int, bitdepth / 8, 1, wavread);
		buffer_int = buffer_int << (64 - bitdepth);
		buffer_int = buffer_int >> (64 - bitdepth);
		buffer_double = buffer_int / bit;
		fwrite(&buffer_double, 8, 1, tmpl);

		fread(&buffer_int, bitdepth / 8, 1, wavread);
		buffer_int = buffer_int << (64 - bitdepth);
		buffer_int = buffer_int >> (64 - bitdepth);
		buffer_double = buffer_int / bit;
		fwrite(&buffer_double, 8, 1, tmpr);
	}

	_fseeki64(tmpl, 0, SEEK_SET);
	_fseeki64(tmpr, 0, SEEK_SET);
	fclose(wavread);
	return true;
}


bool WAV_Filter_Renew(FILE* UpSampleData, FILE* OrigData, unsigned int Times, omp_lock_t* myLock) {
	omp_set_lock(myLock);

	ifstream ifs(".\\FIRFilter.dat");
	string str;
	unsigned __int64 samplesize;
	if (ifs.fail())
	{
		return false;
	}
	double* firfilter_table = new double[4095];
	unsigned __int32 i = 0;
	while (getline(ifs, str))
	{
		firfilter_table[i] = atof(str.c_str());
		i++;
	}

	ifs.close();


	ifstream ifsNS(".\\NoiseShapingCoeff.dat");
	if (ifsNS.fail())
	{
		return false;
	}
	i = 0; int s = 0;
	getline(ifsNS, str);
	unsigned int order = atoi(str.c_str());

	double** NS = new double* [2];
	NS[0] = new double[order];
	NS[1] = new double[order];
	while (getline(ifsNS, str))
	{
		if (str != "0") {
			if (i == 0)
				NS[i][s] = atof(str.c_str());
			else {
				NS[i][order - s - 1] = atof(str.c_str());
			}
			s++;
		}
		else {
			s = 0;
			i++;
		}
	}
	ifsNS.close();
	for (i = 0; i < order; i++) {
		NS[0][i] = NS[0][i] - NS[1][i];
	}

	_fseeki64(OrigData, 0, SEEK_END);
	samplesize = _ftelli64(OrigData);
	_fseeki64(OrigData, 0, SEEK_SET);
	samplesize = samplesize / 8;


	const unsigned int logtimes = unsigned int(log(Times) / log(2));
	const unsigned int fftsize = 4096 * Times;
	const unsigned int datasize = fftsize / 2;
	unsigned int* nowfftsize = new  unsigned int[logtimes];
	unsigned int* zerosize = new  unsigned int[logtimes];
	unsigned int* puddingsize = new  unsigned int[logtimes];
	unsigned int* realfftsize = new unsigned int[logtimes];
	unsigned int* addsize = new  unsigned int[logtimes];
	double** prebuffer = new double* [logtimes];
	double gain = 1;

	double* buffer = new double[fftsize];
	unsigned char* out = new unsigned char[datasize];
	for (i = 0; i < datasize; i++) {
		out[i] = 0;
	}

	double* deltabuffer = new double[order + 1];
	for (i = 0; i < order+1; i++) {
		deltabuffer[i] = 0;
	}

	double x_in = 0;
	double error_y = 0;
	double deltagain = 0.5;
	cufftDoubleReal** fftin;
	cufftDoubleComplex** fftout;
	cufftDoubleReal** ifftout;
	cufftDoubleComplex** ifftin;
	cufftDoubleComplex** firfilter_table_fft;
	Check(cudaMallocManaged((void**)&fftin, logtimes * sizeof(cufftDoubleReal)));
	Check(cudaMallocManaged((void**)&fftout, logtimes * sizeof(cufftDoubleComplex)));
	Check(cudaMallocManaged((void**)&ifftout, logtimes * sizeof(cufftDoubleReal)));
	Check(cudaMallocManaged((void**)&ifftin, logtimes * sizeof(cufftDoubleComplex)));
	Check(cudaMallocManaged((void**)&firfilter_table_fft, logtimes * sizeof(cufftDoubleComplex)));

	cufftHandle* FFT, * iFFT;
	Check(cudaMallocManaged((void**)&FFT, logtimes * sizeof(cufftHandle)));
	Check(cudaMallocManaged((void**)&iFFT, logtimes * sizeof(cufftHandle)));

	unsigned int p = 0;
	unsigned int k = 0;
	unsigned int t = 0;
	unsigned int q = 0;
	for (i = 1; i < Times; i = i * 2) {
		nowfftsize[p] = 4096 * (i * 2);
		realfftsize[p] = nowfftsize[p] / 2 + 1;
		zerosize[p] = nowfftsize[p] / 4;
		puddingsize[p] = realfftsize[p] - 1;
		gain = gain * (2.0 / nowfftsize[p]);
		prebuffer[p] = new double[fftsize];

		Check(cudaMallocManaged((void**)&firfilter_table_fft[logtimes - p - 1], unsigned int(fftsize/i/2+1) * sizeof(cufftDoubleComplex)));
		Check(cudaMallocManaged((void**)&fftin[logtimes - p - 1], unsigned int(fftsize/ i) * sizeof(cufftDoubleReal)));
		Check(cudaMallocManaged((void**)&fftout[logtimes - p - 1], unsigned int(fftsize/i/2+1) * sizeof(cufftDoubleComplex)));
		Check(cudaMallocManaged((void**)&ifftout[logtimes - p - 1], unsigned int(fftsize / i) * sizeof(cufftDoubleReal)));
		Check(cudaMallocManaged((void**)&ifftin[logtimes - p - 1], unsigned int(fftsize/i/2+1) * sizeof(cufftDoubleComplex)));

		for (k = 0; k < fftsize / i; k++) {
			fftin[logtimes - p - 1][k] = 0;
			ifftout[logtimes - p - 1][k] = 0;
		}
		for (k = 0; k < fftsize / i / 2 + 1; k++) {
			ifftin[logtimes - p - 1][k].x = 0;
			ifftin[logtimes - p - 1][k].y = 0;
			fftout[logtimes - p - 1][k].x = 0;
			fftout[logtimes - p - 1][k].y = 0;

		}
		for (k = 0; k < fftsize; k++) {
			prebuffer[p][k] = 0;
		}
		p++;
	}
	
	p = 0;
	for (i = 1; i < Times; i = i * 2) {
		cufftPlan1d(&FFT[logtimes - p - 1], fftsize/i, CUFFT_D2Z, 1);
		cudaDeviceSynchronize();
		cufftPlan1d(&iFFT[logtimes - p - 1], fftsize/i, CUFFT_Z2D, 1);
		cudaDeviceSynchronize();
		p++;
	}

	for (k = 0; k < logtimes; k++) {
		for (i = 0; i < 4095; i++) {
			fftin[k][i] = firfilter_table[i];
		}
		for (i = 4095; i < nowfftsize[k]; i++) {
			fftin[k][i] = 0;
		}
	}

	for (k = 0;k < logtimes;k++) {
		cufftExecD2Z(FFT[logtimes - k - 1], fftin[logtimes - k - 1], fftout[logtimes - k - 1]);
		cudaDeviceSynchronize();
		for (p = 0; p < realfftsize[logtimes - k - 1]; p++) {
			firfilter_table_fft[logtimes - k - 1][p].x = fftout[logtimes - k - 1][p].x;
			firfilter_table_fft[logtimes - k - 1][p].y = fftout[logtimes - k - 1][p].y;
		}
	}

	unsigned __int64 SplitNum = unsigned __int64((samplesize / datasize) * Times);
	deltagain = gain * deltagain;
	omp_unset_lock(myLock);
	for (k = 0; k < SplitNum; k++) {
		fread(buffer, 8, datasize / Times, OrigData);
		for (t = 0;t < logtimes;t++) {
			q = 0;
			for (p = 0; p < zerosize[t]; p++) {
				fftin[t][q] = buffer[p];
				q++;
				fftin[t][q] = 0;
				q++;
			}
			memset(fftin[t] + q, 0, 8 * (nowfftsize[t] - q));
			cufftExecD2Z(FFT[t], fftin[t], fftout[t]);
			cudaDeviceSynchronize();
			for (p = 0; p < realfftsize[t]; p++) {
				ifftin[t][p].x = fftout[t][p].x * firfilter_table_fft[t][p].x - fftout[t][p].y * firfilter_table_fft[t][p].y;
				ifftin[t][p].y = fftout[t][p].x * firfilter_table_fft[t][p].y + firfilter_table_fft[t][p].x * fftout[t][p].y;
			}
			cufftExecZ2D(iFFT[t], ifftin[t], ifftout[t]);
			cudaDeviceSynchronize();
			for (p = 0; p < realfftsize[t]; p++) {
				ifftout[t][p] = ifftout[t][p] / (fftsize / i);
			}
			for (p = 0; p < puddingsize[t]; p++) {
				buffer[p] = prebuffer[t][p] + ifftout[t][p];
			}
			q = 0;
			for (p = puddingsize[t]; p < nowfftsize[t]; p++) {
				prebuffer[t][q] = ifftout[t][p];
				q++;
			}
		}

		for (q = 0; q < datasize; q++) {

			x_in = buffer[q] * deltagain;

			for (t = 0; t < order; t++) {
				x_in += NS[0][t] * deltabuffer[t];
			}

			if (x_in >= 0.0) {
				out[q] = 1;
				error_y = -1.0;
			}
			else {
				out[q] = 0;
				error_y = 1.0;
			}
			for (t = order; t > 0; t--) {
				deltabuffer[t] = deltabuffer[t - 1];
			}

			deltabuffer[0] = x_in + error_y;

			for (t = 0; t < order; t++) {
				deltabuffer[0] += NS[1][t] * deltabuffer[t + 1];
			}
		}
		fwrite(out, 1, datasize, UpSampleData);
	}


	for (i = 0; i < logtimes; i++) {
		delete[] prebuffer[i];
	}

	delete[] NS[0];
	delete[] NS[1];
	delete[] NS;
	delete[] nowfftsize;
	delete[] zerosize;
	delete[] puddingsize;
	delete[] realfftsize;
	delete[] out;
	delete[] prebuffer;
	delete[] buffer;
	delete[] deltabuffer;
	delete[] firfilter_table;
	return true;
}

bool DSD_Write(FILE* LData, FILE* RData, FILE* WriteData, unsigned int DSDSamplingRate, unsigned short bitdepth, unsigned __int32 samplesize, unsigned int Times) {
	unsigned __int64 writelength = samplesize / (bitdepth / 4);
	unsigned __int64 DSD_SampleSize = writelength * Times;
	unsigned __int64 DSD_DataSize = DSD_SampleSize / 4;
	_fseeki64(LData, 0, SEEK_END);
	_fseeki64(RData, 0, SEEK_END);
	_fseeki64(LData, _ftelli64(LData) - DSD_SampleSize, SEEK_SET);
	_fseeki64(RData, _ftelli64(RData) - DSD_SampleSize, SEEK_SET);

	fwrite("FRM8", 4, 1, WriteData);
	unsigned __int64 binary = 0;
	unsigned short ushort = 0;
	unsigned char uchar = 0;
	unsigned __int64 ulong = 0;
	binary = reverse_endian(DSD_DataSize + 152);
	fwrite(&binary, 8, 1, WriteData);
	fwrite("DSD ", 4, 1, WriteData);
	fwrite("FVER", 4, 1, WriteData);
	binary = 0;
	fwrite(&binary, 4, 1, WriteData);
	binary = reverse_endian(4);
	fwrite(&binary, 4, 1, WriteData);


	binary = 1;
	fwrite(&binary, 1, 1, WriteData);
	binary = 5;
	fwrite(&binary, 1, 1, WriteData);
	binary = 0;
	fwrite(&binary, 1, 1, WriteData);
	binary = 0;
	fwrite(&binary, 1, 1, WriteData);

	fwrite("PROP", 4, 1, WriteData);
	binary = 0;
	fwrite(&binary, 4, 1, WriteData);
	binary = reverse_endian(108);
	fwrite(&binary, 4, 1, WriteData);
	fwrite("SND ", 4, 1, WriteData);

	fwrite("FS  ", 4, 1, WriteData);
	binary = 0;
	fwrite(&binary, 4, 1, WriteData);
	binary = reverse_endian(4);
	fwrite(&binary, 4, 1, WriteData);
	unsigned __int32 binary1;
	binary1 = reverse_endian(DSDSamplingRate);
	fwrite(&binary1, 4, 1, WriteData);

	fwrite("CHNL", 4, 1, WriteData);
	binary = 0;
	fwrite(&binary, 4, 1, WriteData);
	binary = reverse_endian(10);
	fwrite(&binary, 4, 1, WriteData);
	binary = 0;
	fwrite(&binary, 1, 1, WriteData);
	binary = 2;
	fwrite(&binary, 1, 1, WriteData);
	fwrite("SLFT", 4, 1, WriteData);
	fwrite("SRGT", 4, 1, WriteData);

	fwrite("CMPR", 4, 1, WriteData);
	binary = 0;
	fwrite(&binary, 4, 1, WriteData);
	binary = reverse_endian(20);
	fwrite(&binary, 4, 1, WriteData);

	fwrite("DSD ", 4, 1, WriteData);
	binary = 14;
	fwrite(&binary, 1, 1, WriteData);
	fwrite("not compressed ", 15, 1, WriteData);

	fwrite("ABSS", 4, 1, WriteData);
	binary = 0;
	fwrite(&binary, 4, 1, WriteData);
	binary = reverse_endian(8);
	fwrite(&binary, 4, 1, WriteData);
	fwrite(&ushort, 2, 1, WriteData);
	fwrite(&uchar, 1, 1, WriteData);
	fwrite(&uchar, 1, 1, WriteData);
	fwrite(&ulong, 4, 1, WriteData);

	fwrite("LSCO", 4, 1, WriteData);
	binary = 0;
	fwrite(&binary, 4, 1, WriteData);
	binary = reverse_endian(2);
	fwrite(&binary, 4, 1, WriteData);
	fwrite(&ushort, 2, 1, WriteData);

	fwrite("DSD ", 4, 1, WriteData);
	binary = reverse_endian(DSD_DataSize);
	fwrite(&binary, 8, 1, WriteData);
	unsigned __int64 i = 0;
	unsigned char* onebyte = new unsigned char[2];
	unsigned char* tmpdataL = new unsigned char[8];
	unsigned char* tmpdataR = new unsigned char[8];

	for (i = 0; i < DSD_SampleSize / 8; i++) {
		fread(tmpdataL, 1, 8, LData);
		fread(tmpdataR, 1, 8, RData);
		onebyte[0] = tmpdataL[0] << 7;
		onebyte[0] += tmpdataL[1] << 6;
		onebyte[0] += tmpdataL[2] << 5;
		onebyte[0] += tmpdataL[3] << 4;
		onebyte[0] += tmpdataL[4] << 3;
		onebyte[0] += tmpdataL[5] << 2;
		onebyte[0] += tmpdataL[6] << 1;
		onebyte[0] += tmpdataL[7] << 0;
		onebyte[1] = tmpdataR[0] << 7;
		onebyte[1] += tmpdataR[1] << 6;
		onebyte[1] += tmpdataR[2] << 5;
		onebyte[1] += tmpdataR[3] << 4;
		onebyte[1] += tmpdataR[4] << 3;
		onebyte[1] += tmpdataR[5] << 2;
		onebyte[1] += tmpdataR[6] << 1;
		onebyte[1] += tmpdataR[7] << 0;
		fwrite(onebyte, 1, 2, WriteData);
	}
	delete[] onebyte;
	delete[] tmpdataL;
	delete[] tmpdataR;
	return true;
}

bool WAV_Convert(const wchar_t* filepath, unsigned int DSD_Times) {
	unsigned __int32 samplingrate;
	unsigned short bitdepth;
	unsigned __int32 samplesize;
	if (!WAV_Metadata(filepath, samplingrate, bitdepth, samplesize))
		return false;
	unsigned int Times;
	unsigned int DSDSamplingRate;
	if (0 == samplingrate % 44100) {
		Times = DSD_Times / (samplingrate / 44100);
		DSDSamplingRate = samplingrate * Times;
	}
	else {
		Times = DSD_Times / (samplingrate / 48000);
		DSDSamplingRate = samplingrate * Times;
	}
	bool flag = true;
	bool flagl = true;
	bool flagr = true;
	FILE* tmpl;
	FILE* tmpr;

	if (!RequireWriteData(filepath, L"_tmpL0", L"wb", &tmpl)) {
		flagl = false;
		flag = false;
	}
	if (!RequireWriteData(filepath, L"_tmpR0", L"wb", &tmpr)) {
		flagr = false;
		flag = false;
	}
	if (flag)if (!TmpWriteData(filepath, bitdepth, samplesize, tmpl, tmpr, Times)) {
		flag = false;
	}
	if (flagl) {
		fclose(tmpl);
	}
	if (flagr) {
		fclose(tmpr);
	}
	omp_lock_t myLock;
	omp_init_lock(&myLock);
#pragma omp parallel
#pragma omp sections
	{
#pragma omp section
		{

			if (flag) {
				bool flagUpl = true;
				bool flagOrigl = true;
				FILE* tmpl;
				FILE* UpsampleDataL;
				omp_set_lock(&myLock);
				if (!RequireWriteData(filepath, L"_tmpL0", L"rb", &tmpl)) {
					flagOrigl = false;
					flag = false;
				}
				if (!RequireWriteData(filepath, L"_tmpLDSD", L"wb", &UpsampleDataL)) {
					flagUpl = false;
					flag = false;
				}
				omp_unset_lock(&myLock);
				if (flag)if (!WAV_Filter_Renew(UpsampleDataL, tmpl, Times, &myLock)) {
					flag = false;
				}
				if (flagUpl) {
					fclose(UpsampleDataL);
				}
				if (flagOrigl) {
					fclose(tmpl);
				}
			}
		}
#pragma omp section  
		{

			if (flag) {
				bool flagUpr = true;
				bool flagOrigr = true;
				FILE* tmpr;
				FILE* UpsampleDataR;
				omp_set_lock(&myLock);
				if (!RequireWriteData(filepath, L"_tmpR0", L"rb", &tmpr)) {
					flagOrigr = false;
					flag = false;
				}
				if (!RequireWriteData(filepath, L"_tmpRDSD", L"wb", &UpsampleDataR)) {
					flagUpr = false;
					flag = false;
				}
				omp_unset_lock(&myLock);
				if (flag)if (!WAV_Filter_Renew(UpsampleDataR, tmpr, Times, &myLock)) {
					flag = false;
				}
				if (flagUpr) {
					fclose(UpsampleDataR);
				}
				if (flagOrigr) {
					fclose(tmpr);
				}
			}
		}
	}
	omp_destroy_lock(&myLock);

	if (flag) {
		FILE* tmpDSD;
		bool flagOrigl = true;
		bool flagOrigr = true;
		bool flagDSD = true;
		if (!RequireWriteData(filepath, L"_tmpLDSD", L"rb", &tmpl)) {
			TrushFile(filepath, L"_tmpLDSD");
			flagOrigl = false;
			flag = false;
		}
		if (!RequireWriteData(filepath, L"_tmpRDSD", L"rb", &tmpr)) {
			TrushFile(filepath, L"_tmpRDSD");
			flagOrigr = false;
			flag = false;
		}
		if (!RequireWriteData(filepath, L".dff", L"wb", &tmpDSD)) {
			TrushFile(filepath, L"_tmpLDSD");
			TrushFile(filepath, L"_tmpRDSD");
			flagDSD = false;
			flag = false;
		}
		if (flag)if (!DSD_Write(tmpl, tmpr, tmpDSD, DSDSamplingRate, bitdepth, samplesize, Times)) {
			TrushFile(filepath, L"_tmpLDSD");
			TrushFile(filepath, L"_tmpRDSD");
		}
		if (flagOrigl) {
			fclose(tmpl);
		}
		if (flagOrigr) {
			fclose(tmpr);
		}
		if (flagDSD) {
			fclose(tmpDSD);
		}
	}

	TrushFile(filepath, L"_tmpL0");
	TrushFile(filepath, L"_tmpR0");
	TrushFile(filepath, L"_tmpLDSD");
	TrushFile(filepath, L"_tmpRDSD");

	if (!flag) {
		return false;
	}
	return true;
}

int wmain(int argc, wchar_t* argv[])
{
	if (argc != 3)
		return -1;
	unsigned int DSD_Times = _wtoi(argv[2]);
	wchar_t cmdflac[260] = LR"(flac -d -f ")";
	wcscat_s(cmdflac, 260, argv[1]);
	wcscat_s(cmdflac, 260, LR"(")");
	STARTUPINFO si;
	PROCESS_INFORMATION pi;
	ZeroMemory(&si, sizeof(si));
	ZeroMemory(&pi, sizeof(pi));
	DWORD flacExitCode;
	DWORD wvExitCode;

	if (CreateProcess(NULL, (LPWSTR)cmdflac, NULL, NULL, false, 0, NULL, NULL, &si, &pi))
	{
		CloseHandle(pi.hThread);
		WaitForSingleObject(pi.hProcess, INFINITE);
		GetExitCodeProcess(pi.hProcess, &flacExitCode);
		CloseHandle(pi.hProcess);
	}
	ZeroMemory(&si, sizeof(si));
	ZeroMemory(&pi, sizeof(pi));
	if (!flacExitCode)
	{
		wchar_t DriveName[3];
		wchar_t DirectoryName[256];
		wchar_t FileName[256];
		_wsplitpath_s(argv[1], DriveName, 3, DirectoryName, 256, FileName, 256, NULL, 0);
		wchar_t ReadPath[260];
		wcscpy_s(ReadPath, 260, DriveName);
		wcscat_s(ReadPath, 260, DirectoryName);
		wcscat_s(ReadPath, 260, FileName);
		wchar_t wavPath[260];
		wcscpy_s(wavPath, 260, ReadPath);
		wcscat_s(wavPath, 260, LR"(.wav)");
		if (WAV_Convert(wavPath, DSD_Times)) {
			DeleteFile(wavPath);
			wchar_t cmdwv[260] = LR"(wavpack -d -h ")";
			wcscat_s(cmdwv, 260, ReadPath);
			wcscat_s(cmdwv, 260, LR"(.dff")");
			if (CreateProcess(NULL, (LPWSTR)cmdwv, NULL, NULL, false, 0, NULL, NULL, &si, &pi))
			{
				CloseHandle(pi.hThread);
				WaitForSingleObject(pi.hProcess, INFINITE);
				GetExitCodeProcess(pi.hProcess, &wvExitCode);
				CloseHandle(pi.hProcess);
			}
			ZeroMemory(&si, sizeof(si));
			ZeroMemory(&pi, sizeof(pi));
			if (!wvExitCode)
				return 0;
		}
	}
}
