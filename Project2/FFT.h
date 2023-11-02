////////////////////////////////////////////////////////////////////
//           FFT based on code from CLRS Algorithms Text          //
////////////////////////////////////////////////////////////////////


#include <complex>
#include <vector>


class FFT
{
public:
    enum  direction { FORWARD, INVERSE };
    static void complex_round(std::vector<std::complex <double>> a);
    static void fft(std::vector<std::complex <double>> a, std::vector<std::complex <double>>& y,
        direction dir);
private:
    static void fftCalculate(std::vector<std::complex <double>> a, std::vector<std::complex <double>>& y,
        direction dir);
};
