////////////////////////////////////////////////////////////////////
//           FFT based on code from CLRS Algorithms Text          //
////////////////////////////////////////////////////////////////////

#define _USE_MATH_DEFINES

#include "FFT.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <algorithm>
#include <complex>
using namespace std;

const double EPSILON = 1.0e-12;
const std::complex<double> I{ 0, 1 };

////////////////////////////////////////////////////////////////////
// complex_round - just set to zero values that are very small    //
//             makes output easier to read                        //
////////////////////////////////////////////////////////////////////
void FFT::complex_round(vector<complex <double>> a)
{
    int n = (int)a.size();

    for (int i = 0; i < n; i++)
    {
        double x = a[i].real();
        double y = a[i].imag();
        if (fabs(x) < EPSILON)
            x = 0;
        if (fabs(y) < EPSILON)
            y = 0;
        a[i] = complex<double>(x, y);
    }
}


////////////////////////////////////////////////////////////////////
//                        FFT based on CLRS                       //
////////////////////////////////////////////////////////////////////
void FFT::fft(vector<complex <double>> a, vector<complex <double>>& y,
    direction dir)
{
    fftCalculate(a, y, dir);

    if (dir == INVERSE)
    {
        int n = (int)a.size();
        for (int i = 0; i < n; i++)
        {
            y[i] /= n;
        }
    }
}

/*
primary FFT function. Finds the FFt for vector a and puts it in y.
*/
void FFT::fftCalculate(vector<complex <double>> a, vector<complex <double>>& y, direction dir)
{
    int n = (int)a.size();
    int n2 = n / 2;
    vector<complex  <double>> even(n / 2);
    vector<complex  <double>> even_fft(n / 2);
    vector<complex  <double>> odd(n / 2);
    vector<complex  <double>> odd_fft(n / 2);
    complex  <double> omega;
    complex  <double> omega_power;

    //base case
    if (n == 1)
    {
        y[0] = a[0];
        return;
    }

    //split the vector into even and odd
    int j = 0;
    for (int i = 0; i < n; i += 2)
    {
        even[j] = a[i];
        odd[j] = a[i + 1];
        j++;
    }

    //recursion
    fftCalculate(even, even_fft, dir);
    fftCalculate(odd, odd_fft, dir);

    //determine which rotation to use depending on direction
    if (dir != FORWARD)
        omega = cos(-2.0 * M_PI / n) + I * sin(-2.0 * M_PI / n);
    else
        omega = cos(2.0 * M_PI / n) + I * sin(2.0 * M_PI / n);

    //recombine even and odd
    omega_power = 1;
    for (int i = 0; i < n2; i++)
    {
        y[i] = (even_fft[i] + omega_power * odd_fft[i]);
        y[i + n2] = (even_fft[i] - omega_power * odd_fft[i]);
        omega_power = omega * omega_power;
    }
}

