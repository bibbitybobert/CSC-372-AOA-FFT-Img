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
vector<unsigned char> convertToUnsignedChar(complex1D v)
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

ifstream openInFile(string file_name) {
    ifstream inFile;
    inFile.open(file_name, ios::in | ios::binary);
    if (!inFile) {
        cout << "Error opening input file: " << file_name;
        cout << "quitting program\n";
        exit(1);
    }
    return inFile;
}

ofstream openOutFile(string file_name) {
    ofstream outFile;
    outFile.open(file_name, ios::out | ios::binary);
    if (!outFile) {
        cout << "Error opening output file: " << file_name;
        cout << "quitting program\n";
        exit(1);
    }
    return outFile;
}

vector<unsigned char> loadImage(string fileName, unsigned int& width, unsigned int& height) {
    vector<unsigned char> outVec;
    string image_type;
    ifstream in_file;

    image_type = fileName.substr(fileName.size() - 3, 3);
    in_file = openInFile(fileName);

    if (image_type == "png") {
        outVec = loadPNG(fileName, width, height);
    }
    else if (image_type == "bmp") {
        outVec = loadBMP(fileName, width, height);
    }

    return outVec;
}

void saveImage(string fileName, unsigned int width, unsigned int height, vector<unsigned char> image) {
    string image_type;
    ofstream out_file;

    image_type = fileName.substr(fileName.size() - 3, 3);
    out_file = openOutFile(fileName);
    
    if (image_type == "png") {
        savePNG(fileName, width, height, image);
    }
    else if (image_type == "bmp") {
        saveBMP(fileName, width, height, image);
    }
}

void applyBands(complex2D& image, int band_size) {
    complex<double> zero = complex<double>(0);

    int start_loc = (image.size() / 2) - (band_size / 2); //from index
    int end_loc = image.size() - start_loc; //to index (exclusive)

    vector<complex<double>> temp;
    temp.resize(image.size(), zero);
    for (int i = 0; i < image.size(); i++) {
        if (i >= start_loc && i < end_loc) { //convert to all black rows
            image[i] = temp;
        }
        else {
            for (int j = start_loc; j < end_loc; j++) {
                image[i][j] = zero;
            }
        }
    }
}

vector<unsigned char> normalize(complex1D image) {
    double max = 0.0;
    vector<unsigned char> normalizedImage(image.size());
    for (int i = 1; i < image.size(); i++) {
        double mag = sqrt(image[i].real() * image[i].real() + image[i].imag() * image[i].imag());
        if (mag > max)
            max = mag;
    }

    for (int i = 0; i < image.size(); i++) {
        double mag = sqrt(image[i].real() * image[i].real() + image[i].imag() * image[i].imag());
        unsigned char normalizedVal = (mag / max) * 255;
        normalizedImage[i] = normalizedVal;
    }
    return normalizedImage;
}

int main(int argc, char* argv[])
{

    vector<unsigned char> origImage;
    vector<unsigned char> outImage;
    vector<unsigned char> interImage;
    unsigned int width = 0;
    unsigned int height = 0;
    

    //example aliasing use
    complex1D data;
    complex1D result(8); //make it size 8
    complex1D interData;

    complex2D data2D;
    complex2D result2D; 
  
    
    //load image
    //origImage = {1, 2, 3, 4, 5, 6, 7, 8};
    origImage = loadImage(argv[1], width, height);

    data = convertToComplex(origImage);
    result.resize(data.size());
    

    data2D = makeComplex2D(data, width, height);
    result2D = complex2D(data2D.size(), complex1D(data2D.size()));

    // Do forward FFT
    FFT::fft_2D(data2D, result2D, FFT::FORWARD);

    //apply bands 
    applyBands(result2D, stoi(argv[4]));
    interData = makeComplex1D(result2D);
    interImage = normalize(interData);
    
    saveImage(argv[2], width, height, interImage);


    // do inverse FFT
    FFT::fft_2D(result2D, data2D, FFT::INVERSE);

    data = makeComplex1D(data2D);
    outImage = convertToUnsignedChar(data);

    saveImage(argv[3], width, height, outImage);


    return 0;
}

