//Hao Zhi Rao 30031059
//
//compiled with: 
// no complier optimizations
// g++ -o convolve convolve.cpp -pg -no-pie
// with compiler optimize
// g++ -o convolve convolve.cpp -O3 -pg -no-pie

//to time results
//(time ./convolve guitar.wav bighall.wav out.wav ) &> result.txt

//profiler used, gprof
//gprof -b convolve > gprof.txt

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>

using namespace std;

//#define DEBUG
#define FAST 1
#define MAX_SHORT 32767
#define MONO 1
#define FFT 1
#define IFFT -1

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE 16
/*  Standard sample rate in Hz  */
#define SAMPLE_RATE 44100.0
/*  Standard sample size in bytes  */
#define BYTES_PER_SAMPLE (BITS_PER_SAMPLE / 8)

//for four1
#define SIZE 8
#define PI 3.141592653589793
#define TWO_PI (2.0 * PI)
#define SWAP(a, b) tempr = (a);(a) = (b);(b) = tempr

void convolve(double x[], int N, double h[], int M, double y[], int P);
int checkHeader(FILE *wavFile);
void writeWaveFileHeader(int channels, int numberSamples,double outputRate, FILE *outputFile);
size_t fwriteIntLSB(int data, FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);
void four1(double data[], int nn, int isign);

double* fftconvolve(double x[], int N, double h[], int M);
int powerTwo(int num);
double *padArray(double x[], int N, int N2);

int main(int argc, char *argv[])
{
    char *inputFile = NULL, *irFile = NULL, *outFile = NULL;

    if (argc != 4)
    {
        cout << "Usage: " << argv[0] << " inputfile IRfile outputfile" << endl;
        return -1;
    }
    else
    {
        inputFile = argv[1];
        irFile = argv[2];
        outFile = argv[3];
    }
    FILE *input = fopen(inputFile, "rb");
    FILE *ir = fopen(irFile, "rb");
    FILE *output = fopen(outFile, "wb");

    int numSamplesInput = checkHeader(input);

    int numSamplesIR = checkHeader(ir);

    int ySize = numSamplesInput + numSamplesIR - 1;

    double *x = new double[numSamplesInput];

    short readbits;
    //partial unrolling
    int i;
    for (i = 0; i < numSamplesInput-1; i+=2)
    {
        fread(&readbits, 2, 1, input);
        x[i] = readbits;
        fread(&readbits, 2, 1, input);
        x[i+1] = readbits;
    }
    if(i == numSamplesInput-1){
        fread(&readbits, 2, 1, input);
        x[numSamplesInput-1] = readbits;
    }

    double *h = new double[numSamplesIR];

    //partial unrolling again
    for (i = 0; i < numSamplesIR-1; i+=2)
    {
        fread(&readbits, 2, 1, ir);
        h[i] = readbits;
        fread(&readbits, 2, 1, ir);
        h[i+1] = readbits;
    }
    if(i == numSamplesIR -1){
        fread(&readbits, 2, 1, ir);
        h[numSamplesIR -1] = readbits;
    }

    double *y;
    if(FAST)
        y = fftconvolve(x,numSamplesInput, h, numSamplesIR);
    else{
    y = new double[ySize];
     convolve(x, numSamplesInput, h, numSamplesIR, y, ySize);
}
    writeWaveFileHeader(MONO, ySize, SAMPLE_RATE, output);

    double minVal = 0.0;
    double maxVal = 0.0;
    double tempVal;
    for (int i = 0; i < ySize; i++)
    {
        tempVal = y[i];

        if (tempVal > maxVal)
            maxVal = tempVal;
        if (tempVal < minVal)
            minVal = tempVal;
    }
   
    //Jamming
    //combine normalize into writing
    //write back as approp 16 bit value
    short sampleVal;
    for (i = 0; i < ySize-1; i+=1)
    {// plus another partial unrolling 
        sampleVal = (short)rint(((y[i] - minVal) / (maxVal - minVal) * MAX_SHORT));
        fwriteShortLSB(sampleVal, output);
    }
    if(i == ySize -1){
        sampleVal = (short)rint(((y[ySize -1] - minVal) / (maxVal - minVal) * MAX_SHORT));
        fwriteShortLSB(sampleVal, output);
    }

    fclose(output);
    fclose(input);
    fclose(ir);
}

//strength reduction
//changed mul by 2 and div by 2 to respective shifts
double* fftconvolve(double x[], int N, double h[], int M)
{
    //get new indices that are pow2 + imaginary part
    // x[] MUST be the input ie longest
    int length =  powerTwo(N) << 1;

    double *xf = padArray(x, N, length);
    double *hf = padArray(h, M, length);

    four1(xf - 1, length >> 1, FFT);
    four1(hf - 1, length >> 1, FFT);

    //complex multiplication to get result

    double *yf = new double[length];

    for (int i = 0; i < length; i += 2)
    {
        yf[i] =     xf[i] * hf[i]     - xf[i + 1] * hf[i + 1];
        yf[i + 1] = xf[i + 1] * hf[i] + xf[i] * hf[i + 1];
    }

    //get inverse fft
     four1(yf - 1, length >> 1, IFFT);

    double * ret = new double[length];
    for(int i = 0; i < length; i++){
        ret[i] = yf[i << 1];
    }


    return ret;

}

//round up power two
//http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
int powerTwo(int num)
{
    unsigned int v = num;
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}
//minimzing work done in array
//reorders data, real is signal, imag is zeroed
double* padArray(double x[], int N, int N2)
{
    double *padded = new double[N2];
    memset(padded, 0, sizeof(padded));
    for (int k = 0, i = 0; k < N; k++, i += 2)
    {
        if (k < N)
            padded[i] = x[k];
        else
            padded[i] = 0;
        padded[i + 1] = 0;
    }
    return padded;
}
/*
* taken from the fft demo, test.c
*/

//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below).

void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2)
    {
        if (j > i)
        {
            SWAP(data[j], data[i]);
            SWAP(data[j + 1], data[i + 1]);
        }
        m = nn;
        while (m >= 2 && j > m)
        {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    mmax = 2;
    while (n > mmax)
    {
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;   
        wi = 0.0;
        for (m = 1; m < mmax; m += 2)
        {
            for (i = m; i <= n; i += istep)
            {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j + 1];
                tempi = wr * data[j + 1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j + 1] = data[i + 1] - tempi;
                data[i] += tempr;
                data[i + 1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}

/* 
* taken from the convolution demo, convolve.c 
*/

void convolve(double x[], int N, double h[], int M, double y[], int P)
{
    int n, m;
    /*  Make sure the output buffer is the right size: P = N + M - 1  */
    if (P != (N + M - 1))
    {
        printf("Output signal vector is the wrong size\n");
        printf("It is %-d, but should be %-d\n", P, (N + M - 1));
        printf("Aborting convolution\n");
        return;
    }
    /*  Clear the output buffer y[] to all zero values  */
    for (n = 0; n < P; n++)
        y[n] = 0.0;

    /*  Do the convolution  */
    /*  Outer loop:  process each input value x[n] in turn  */
    for (n = 0; n < N; n++)
    {
        /*  Inner loop:  process x[n] with each sample of h[]  */
        for (m = 0; m < M; m++)
            y[n + m] += x[n] * h[m];
    }
}

/*
*  methods taken from testtone.c
* Author:        Leonard Manzara
* Date:          November 21, 2009
* writeWaveFileHeader() , fwriteIntLSB, fwriteShortLSB 
*/

void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, FILE *outputFile)
{
    /*  Calculate the total number of bytes for the data chunk  */
    int dataChunkSize = channels * numberSamples * BYTES_PER_SAMPLE;

    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;

    /*  Calculate the total number of bytes per frame  */
    short int frameSize = channels * BYTES_PER_SAMPLE;

    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);

    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);

    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);

    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);

    /*  Number of channels  */
    fwriteShortLSB((short)channels, outputFile);

    /*  Output Sample Rate  */
    fwriteIntLSB((int)outputRate, outputFile);

    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    fwriteShortLSB(BITS_PER_SAMPLE, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
}

/******************************************************************************
*
*       function:       fwriteIntLSB
*
*       purpose:        Writes a 4-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*
*       internal
*       functions:      none
*
*       library
*       functions:      fwrite
*
******************************************************************************/

size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];
    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}

/******************************************************************************
*
*       function:       fwriteShortLSB
*
*       purpose:        Writes a 2-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*
*       internal
*       functions:      none
*
*       library
*       functions:      fwrite
*
******************************************************************************/

size_t fwriteShortLSB(short int data, FILE *stream)
{
    unsigned char array[2];
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}

/*
* checks file header, currently does nothing really
*/
int checkHeader(FILE *wavFile)
{
    char string[5];
    string[4] = 0;

    fread(string, 4, 1, wavFile);
#ifdef DEBUG
    cout << string << endl;
#endif

    int chunkSize;
    fread(&chunkSize, 4, 1, wavFile);

    fread(string, 4, 1, wavFile);
#ifdef DEBUG
    cout << string << endl;
#endif
    fread(string, 4, 1, wavFile);
#ifdef DEBUG
    cout << string << endl;
#endif

    int subChunk1Size;
    fread(&subChunk1Size, 4, 1, wavFile);

    short int audioFormat;
    fread(&audioFormat, 2, 1, wavFile);

    short int numChannels;
    fread(&numChannels, 2, 1, wavFile);

    int sampleRate;
    fread(&sampleRate, 4, 1, wavFile);

    int byteRate;
    fread(&byteRate, 4, 1, wavFile);

    short int blockAlign;
    fread(&blockAlign, 2, 1, wavFile);

    short int bitsPerSample;
    fread(&bitsPerSample, 2, 1, wavFile);

    if (subChunk1Size != 16)
    {
        short temp;
        fread(&temp, 2, 1, wavFile);
        cout << "extra params " << temp << endl;
        fread(NULL, temp, 1, wavFile);
    }

    fread(string, 4, 1, wavFile);
#ifdef DEBUG
    cout << string << endl;
#endif

    int subChunk2Size;

    fread(&subChunk2Size, 4, 1, wavFile);
#ifdef DEBUG
    cout << chunkSize << endl;
    cout << subChunk1Size << endl;
    cout << audioFormat << endl;
    cout << numChannels << endl;
    cout << sampleRate << endl;
    cout << byteRate << endl;
    cout << blockAlign << endl;
    cout << bitsPerSample << endl;
    cout << subChunk2Size << endl;
#endif
    // x # bytes, but samples are shorts, so return x/2 # of shorts
    return subChunk2Size / 2;
}
