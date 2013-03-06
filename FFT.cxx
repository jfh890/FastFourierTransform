/*
 * FFT.cxx
 *
 *  Created on: Feb 18, 2011
 *      Author: dr
 */

// Global includes
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <complex>
#include <algorithm>
// Local includes
#include "FFT.h"

using namespace std;

#define ITERATIVE

Image<double> FFT::Magnitude(Image<complexd> input) {
	Image<double> out(input.numRows(), input.numCols());
	int size = input.size();
	for( int i = 0; i < size; i++ ) {
		out[i] = abs(input[i]);
	}
	return out;
}


/* Used for displaying image with applied FT */
Image<double> FFT::logOp(Image<complexd> input) {
	Image<double> out(input.numRows(), input.numCols());
	int size     = input.size();
	complexd max = findMax(input);
	int c        = 255/(log(1 + abs(max)));
	for( int i = 0; i < size; i++ ) {
		out[i] = c*log(1+abs(input[i]));
	}
	return out;
}

bool magnitude (const complexd & lhs, const complexd & rhs) {
   return abs(lhs) < abs(rhs);
}

complexd FFT::findMax(Image<complexd> input) {
	return *max_element(input.begin(), input.end(), &magnitude);
}

void FFT::Filter(Image<complexd> &image, int type, double maxcut, double mincut) {
	int rows = image.numRows();
	int cols = image.numCols();
	double cornerToCentre = sqrt(pow(cols/2,2)+pow(rows/2,2));

	if (type == 1) {
		cout << " low-pass ..." << endl;
	} else if (type == 2) {
		cout << " high-pass ..." << endl;
	} else if (type == 3) {
		cout << " band-pass ..." << endl;
	}

	double minDistance = mincut * cornerToCentre;
	double maxDistance = maxcut * cornerToCentre;
	for(int i = 0; i < rows; i++) {
		for(int k = 0; k < cols; k++) {
			if(( type == 3 || type == 2 ) && withinDistanceToCorner(i, k, rows, cols, minDistance)) {
				image.setData(i, k, complexd(0,0));
			} else if(( type == 3 || type == 1 ) && !withinDistanceToCorner(i, k, rows, cols, maxDistance)) {
				image.setData(i, k, complexd(0,0));
			}
		}
	}
}

bool FFT::withinDistanceToCorner(int row, int col, int maxRow, int maxCol, double dist) {
	bool withinCorner1 = sqrt(pow(row,2)+pow(col,2)) <= dist;
	bool withinCorner2 = sqrt(pow(row,2)+pow(col-maxCol,2)) <= dist;
	bool withinCorner3 = sqrt(pow(row-maxRow,2)+pow(col,2)) <= dist;
	bool withinCorner4 = sqrt(pow(row-maxRow,2)+pow(col-maxCol,2)) <= dist;
	return withinCorner1 || withinCorner2 || withinCorner3 || withinCorner4;
}

Image<complexd> FFT::imgToComplex(Image<double> image) {
	Image<complexd> compImage(image.numRows(), image.numCols());
	int size = image.size();
	for(int i = 0; i < size; i++) {
		compImage[i] = complexd(image[i], 0);
	}
	return compImage;
}


Image<complexd> FFT::ForwardTransform(Image<double> image) {
	Image<complexd> compImage = imgToComplex(image);
	Run_2D(compImage, 1);
	return compImage;
}

Image<double> FFT::ReverseTransform(Image<complexd> image) {
	Run_2D(image, -1);
	Image<double> img = Magnitude(image);
	return img;
}

void FFT::Run_2D(Image<complexd> &image, int dir) {
	int rows = image.numRows();
	int cols = image.numCols();

	for(int i = 0; i < cols; i++) {
		vector<complexd> col = image.extractCol(i);
		#ifdef ITERATIVE
			vector<complexd> colFT = Run_1D_Iterative(col, rows, dir);
		#else
			vector<complexd> colFT = Run_1D_Recursive(col, rows, dir);
		#endif
		image.fillCol(i, colFT);
	}

	for(int i = 0; i < rows; i++) {
		vector<complexd> row = image.extractRow(i);
		#ifdef ITERATIVE
			vector<complexd> rowFT = Run_1D_Iterative(row, cols, dir);
		#else
			vector<complexd> rowFT = Run_1D_Recursive(row, cols, dir);
		#endif
		image.fillRow(i, rowFT);
	}
}

vector<complexd> FFT::Run_1D_Recursive(const vector<complexd> &v, int n, int dir) {
	vector<complexd> retY(v);
	if (n == 1) {
		return retY;
	}

	complexd Wn     = exp(complexd(0, 2 * M_PI / n));
	complexd dirMod(1,0);
	complexd w(1,0);
	if( dir == -1 ) {
		Wn     = complexd(1,0) / Wn;
		dirMod = complexd(2,0);
	}

	vector<complexd> vEven    = getMembers(v, true);
	vector<complexd> vOdd     = getMembers(v, false);
	vector<complexd> retYEven = Run_1D_Recursive(vEven, n/2, dir);
	vector<complexd> retYOdd  = Run_1D_Recursive(vOdd, n/2, dir);

	for( int k = 0; k < n/2; k++ ) {
		complexd t    = w * retYOdd[k];
		complexd u    = retYEven[k];
		retY[k]       = (u + t) / dirMod;
		retY[k+(n/2)] = (u - t) / dirMod;
		w             *= Wn;
	}
	return retY;
}

vector<complexd> FFT::Run_1D_Iterative(const vector<complexd> &v, int n, int dir) {
	vector<complexd> bitRevCopy(n);
	BitReverseCopy(v, bitRevCopy, n);
	double upperLoopBound = log2(n);
	for (int s = 1; s <= upperLoopBound; s++){
		int m = pow(2,s);
		complexd Wm     = exp(complexd(0, 2 * M_PI / m));
		complexd dirMod(1,0);
		if (dir == -1) {
			Wm     = complexd(1,0) / Wm;
			dirMod = complexd(2,0);
		}
		for (int k = 0; k < n; k += m){
			complexd w(1,0);
			for (int i = 0; i < m/2; i++){
				complexd t          = w * bitRevCopy[k+i+m/2];
				complexd u          = bitRevCopy[k+i];
				bitRevCopy[k+i]     = (u + t) / dirMod;
				bitRevCopy[k+i+m/2] = (u - t) / dirMod;
				w                   *= Wm;
			}
		}
	}
	return bitRevCopy;
}

void FFT::BitReverseCopy(const vector<complexd> &v, vector<complexd> &bitRevCopy, int n) {
	unsigned int maxbits = (unsigned int) log2(n);
	for(int i = 0; i < n; i++) {
		bitRevCopy[revBits(i,maxbits)] = v[i];
	}
}

unsigned int FFT::revBits(unsigned int n, int size) {
	unsigned int result = 0;

	for (unsigned int i = size; i; --i) {
		result <<= 1;
		if (n & 1) {
			result |= 1;
		}
		n >>= 1;
	}
	return result;
}

vector<complexd> FFT::getMembers(const vector<complexd> &v, bool even) {
	int size = v.size()/2;
	vector<complexd> retV(size);
	int modifier = even ? 0 : 1;
	for(int i = 0; i < size; i++) {
		retV[i] = v[i*2+modifier];
	}
	return retV;
}