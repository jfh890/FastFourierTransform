/*
 * FFT.h
 *
 *  Created on: Feb 18, 2011
 *      Author: dr
 */

#ifndef FFT_H_
#define FFT_H_

// Global includes
#include <string>
#include <math.h>
#include <vector>
#include <complex>

// Local includes
#include "Image.h"

using namespace std;

class FFT
{

private:

  // Run 2D FFT in the direction specified
  void Run_2D(Image<complexd> &image, int dir);

  // Run 1D FFT in the direction specified (uses a recursive implementation)
  vector<complexd> Run_1D_Recursive(const vector<complexd> &v, int n, int dir);

  // Run 1D FFT in the direction specified (uses an iterative implementation)
  vector<complexd> Run_1D_Iterative(const vector<complexd> &v, int n, int dir);

  // Extracts even or odd memebers from vector
  vector<complexd> getMembers(const vector<complexd> &v, bool even);

  // Compute magnitude of complex image
  Image<double> Magnitude(Image<complexd> c);

  // Find maximum value in image (for logOp)
  complexd findMax(Image<complexd> input);

  // checks whether specified point is within dist from any of the corners
  bool withinDistanceToCorner(int row, int col, int maxRow, int maxCol, double dist);

  // Creates bit reverse copy ordering of the vector
  void BitReverseCopy(const vector<complexd> &v, vector<complexd> &bitRevCopy, int n);

  // Reverses bits in integer
  unsigned int revBits(unsigned int n, int size);

  // Converts image of doubles to image of complex<double>
  Image<complexd> imgToComplex(Image<double> image);

public:

  // Constructor
  FFT();

  // logarithmic operator to display image with applied fourier transform
  Image<double> logOp(Image<complexd> image);

  // Filter an image using low pass, high pass and band pass filters. The
  // parameters max and min specify the frequencies for the filters
  void Filter(Image<complexd> &image, int type, double max, double min);

  // Run forward 2D FFT on a 2D image of type double. Returns an image of
  // complex doubles
  Image<complexd> ForwardTransform(const Image<double> image);

  // Run reverse 2D FFT on a 2D image of type complex double. Returns an image
  // doubles
  Image<double> ReverseTransform(const Image<complexd> image);

};

// Constructor with not much to do
inline FFT::FFT()
{}

#endif /* FFT_H_ */
