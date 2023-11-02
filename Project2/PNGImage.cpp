#include "PNGImage.h"
#include "lodepng.h"
#include <vector>
#include <iostream>

using namespace std;

/*
load method for PNG image that automatticaly converts to greyscale
It will return a vector of unsigned chars of the colors in the PNG image convrted to grayscale..
width and height will be updated to the images width and height.
*/
vector<unsigned char> loadPNG(string filename, unsigned int& width, unsigned int& height)
{
    vector<unsigned char> image;

    //decode
    unsigned error = lodepng::decode(image, width, height, filename);

    if (error)
    {
        cout << "decoder error while reading file " << filename << endl;;
        cout << "error code" << error << ": " << lodepng_error_text(error) << "\n";
        exit(2);
        return{};
    }

    vector<unsigned char> greyImage(width * height);
    for (unsigned int i = 0; i < width * height * 4; i += 4)
    {
        greyImage[i / 4] = (int)(0.3 * image[i] + 0.59 * image[i + 1] + 0.11 * image[i + 2]);
    }

    return greyImage;
}

/*
save method for a PNG. It will save a PNG image of size width by height using the vector
of grey scale values given. The resulting image will have all three channels (RGB)
*/
void savePNG(string filename, unsigned int width, unsigned int height, vector<unsigned char> image)
{
    //make it 3 channels again
    vector<unsigned char> threeChannel(width * height * 4);
    for (unsigned int i = 0; i < image.size(); i++)
    {
        threeChannel[i * 4] = image[i];
        threeChannel[i * 4 + 1] = image[i];
        threeChannel[i * 4 + 2] = image[i];
        threeChannel[i * 4 + 3] = 255;
    }

    //Encode the image
    unsigned error = lodepng::encode(filename, threeChannel, width, height);

    //if there's an error, display it
    if (error)
    {
        cout << "encoder error " << error << ": " << lodepng_error_text(error) << endl;
    }
}
