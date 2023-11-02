#pragma once

#include <vector>
#include <string>

std::vector<unsigned char> loadBMP(std::string filename, unsigned int& width, unsigned int& height);
void saveBMP(std::string filename, unsigned int width, unsigned int height, std::vector<unsigned char> img);
