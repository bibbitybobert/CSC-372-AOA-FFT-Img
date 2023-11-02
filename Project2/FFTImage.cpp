#include "FFT.h"
#include "BMPImage.h"
#include "PNGImage.h"
#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <cmath>

using namespace std;

const bool VERBOSE = false;

typedef vector<complex<double>> complex1D;
typedef vector<vector<complex<double>>> complex2D;


/*
Helper function to convert a 1D vector of double type, complex values to a w by h 2D
vector of double type, complex values
*/
complex2D makeComplex2D(complex1D v, unsigned int w, unsigned int h)
{
    complex2D newVector(w, complex1D(h));

    for (unsigned int i = 0; i < w; i++)
    {
        for (unsigned int j = 0; j < h; j++)
        {
            newVector[i][j] = v[i * w + j];
        }
    }
    return newVector;
}

/*
Helper function to convert a 2D vector of double type, complex values 1D
vector of double type, complex values in row-column order
*/
complex1D makeComplex1D(complex2D v)
{
    size_t w = v.size();
    size_t h = v[0].size();

    complex1D  newV(w * h);

    for (unsigned int i = 0; i < w; i++)
    {
        for (unsigned int j = 0; j < h; j++)
        {
            newV[i * h + j] = v[i][j];
        }
    }
    return newV;
}

/*
Helper function to convert a 1D vector of unsigned chars (pixel values) to
a 1D vector of double type, complex values
*/
complex1D convertToComplex(vector<unsigned char> v)
{
    complex1D newV(v.size());

    for (unsigned int i = 0; i < v.size(); i++)
    {
        newV[i] = complex<double>(v[i]);
    }
    return newV;
}

/*
Helper function to convert a 1D vector of double type, complex values  to
a 1D vector of unsigned chars (pixel values). This pulls only the "magnitude"
values from the complex number. It has minor bounds checking where magnitudes
over 255 are clamped to 255. If unnormalized complex arrays are sent in,
the resulting chars can all be 255.
*/
vector<unsigned char> convertToUnsignedChar(complex1D v, int n)
{
    vector<unsigned char> newV(v.size());

    double magnitude;
    for (unsigned int i = 0; i < v.size(); i++)
    {
        magnitude = sqrt(v[i].real() * v[i].real() + v[i].imag() * v[i].imag());

        if (magnitude > 255)
            magnitude = 255;

        newV[i] = (unsigned char)magnitude;
    }

    return newV;
}





int main(int argc, char* argv[])
{

    vector<unsigned char> origImage;
    unsigned int width;
    unsigned int height;

    //example aliasing use
    complex1D data;
    complex1D result(8); //make it size 8
    complex2D data2D;
    complex2D result2D = complex2D(data.size(), complex1D(data.size()));
  
    
    //example FFT calls
    origImage = {1, 2, 3, 4, 5, 6, 7, 8};
    data = convertToComplex(origImage);
    // Do forward FFT
    FFT::fft(data, result, FFT::FORWARD);

    // do inverse FFT
    FFT::fft(result, data, FFT::INVERSE);


    return 0;
}
