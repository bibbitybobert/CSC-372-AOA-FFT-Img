#ifndef _loadpng_h
#define _loadpng_h

#include <vector>
#include <string>

std::vector<unsigned char> loadPNG(std::string filename, unsigned int& width, unsigned int& height);
void savePNG(std::string filename, unsigned int width, unsigned int height, std::vector<unsigned char> img);

#endif
