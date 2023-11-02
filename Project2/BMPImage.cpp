#include "BMPImage.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

constexpr int HEADER_LENGTH = 54;
constexpr int DIB_HEADER_LENGTH = 40;
constexpr bool DEBUG = false;

/*
load method for 24 bit BMP images (no other are supported)
It will return a vector of unsigned chars of the colors in the BMP image converted to grayscale..
width and hieght will be updated to the images width and height.
*/
vector<unsigned char> loadBMP(string filename, unsigned int& width, unsigned int& height)
{


    ifstream f(filename, std::ios::binary);

    if (f.fail())
    {
        throw "Argument Exception";
        return{ };
    }

    char info[HEADER_LENGTH];
    f.read(info, HEADER_LENGTH);

    //check some header info to make sure it is probably valid
    char bChar = *(char*)&info[0];
    char mChar = *(char*)&info[1];
    width = *(int*)&info[18];
    height = *(int*)&info[22];
    short depth = *(short*)&info[28];

    if ((bChar != 'B') || (mChar != 'M'))
    {
        throw "Invalid header in file";
        return{ };
    }

    if (depth != 24)
    {
        throw "Only 24 bit encodings are supported";
        return{ };
    }

    if (DEBUG)
    {
        cout << endl;
        cout << "  Name: " << filename << endl;
        cout << " Width: " << width << endl;
        cout << "Height: " << height << endl;
    }

    //bmp file MUST have rows evenly divisable by 4
    int row_padded = width * 3;
    if (row_padded % 4 != 0)
        row_padded += 4 - row_padded % 4;


    //read in data one row at a time
    char* data = new char[row_padded];
    vector<unsigned char> img(width * height);
    unsigned char tmp;
    for (int i = height - 1; i >= 0; i--)
    {
        //grab entire row and then parse it
        f.read(data, row_padded);
        for (unsigned int j = 0; j < width * 3; j += 3)
        {
            // Convert (B, G, R) to (R, G, B)
            tmp = data[j];
            int R = (unsigned char)data[j + 2];
            int G = (unsigned char)data[j + 1];
            int B = tmp;

            //put grayscale value into image index
            int index = i * width + j / 3;
            img[index] = (int)(0.3 * R + 0.59 * G + 0.11 * B);

            if (DEBUG)
            {
                cout << "(" << i << ", " << j / 3 << ") R: " << R << " G: " << G << " B: " << B << endl;
            }
        }
    }

    f.close();

    delete[] data;
    return img;
}


/*
save method for a BMP. It will save a BMP image of size width by height using the vector
of grey scale values given. The resulting image will have all three channels (RGB)
*/
void saveBMP(string filename, unsigned int width, unsigned int height, vector<unsigned char> img)
{
    ofstream f(filename, ios::out | ios::binary);

    // Each pixel is stored as three bytes.  One for each of R, G, and B.
    int lineLength = width * 3;

    //  The BMP format specifies that the number of bytes to represent a
    //  line must be a multiple of four bytes.
    if (lineLength % 4 != 0)
        lineLength += 4 - lineLength % 4;

    f.write("B", 1);
    f.write("M", 1);

    // Determine the file size, which is HEADER_LENGTH bytes plus the
    // amount of data in the image, and thenInsert the file size (as an int)
    int fileSize = HEADER_LENGTH + (lineLength * height);
    f.write((char*)&fileSize, 4);

    //Insert a 32 bit int with a value of 0. (The value is not used)
    int pad = 0;
    f.write((char*)&pad, 4);

    //Insert the offset where the pixel array starts (as an int).
    f.write((char*)&HEADER_LENGTH, 4);

    // Insert the size of the DIB header (as an int)
    f.write((char*)&DIB_HEADER_LENGTH, 4);

    //  Insert the width and height of the image (as ints).
    f.write((char*)&width, 4);
    f.write((char*)&height, 4);

    // Insert a 1 for 1 color plane (as a short)
    short plane = 1;
    f.write((char*)&plane, 2);

    // Insert a 24 for 24 bits per pixel. (as a short)
    short bitSize = 24;
    f.write((char*)&bitSize, 2);

    // Insert a zero to indicate no pixel array compression (as an int)
    f.write((char*)&pad, 4);

    // Insert the size of the raw bitmap data (as an int)
    int rawSize = lineLength * height;
    f.write((char*)&rawSize, 4);

    // Insert four 32 bit 0 values for unused fields (as ints)
    f.write((char*)&pad, 4);
    f.write((char*)&pad, 4);
    f.write((char*)&pad, 4);
    f.write((char*)&pad, 4);

    //Insert the pixel data into the byte buffer
    //    a) Loop over all rows
    for (int r = height - 1; r >= 0; r--)
    {


        // b) Loop over each pixel in the row
        for (unsigned int c = 0; c < width; c++)
        {
            // Get pixel from the image
            unsigned char temp = img[r * width + c];

            //Normally, insert blue, green, and then red bytes (as bytes) (but this is grayscale)
            f.write((char*)&temp, 1);
            f.write((char*)&temp, 1);
            f.write((char*)&temp, 1);
        }


        // c) Add padding so that the number of bytes for each row is
        //    evenly divisible by 4.  Add additional bytes with the value
        //    of 0 for the padding bytes.
        if (width * 3 % 4 != 0)
        {
            for (unsigned int pad = 0; pad < (4 - width * 3 % 4); pad++)
            {
                f.write((char*)&pad, 1);
            }
        }
    }
    f.close();
}
