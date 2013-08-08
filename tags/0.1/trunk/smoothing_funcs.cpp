/*

smoothing_funcs.cpp provides functions to perform non-local means filtering following Baudes et al (2005)

Copyright (C) 2013 Martin D. Hurst

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Martin D. Hurst
mhurst@bgs.ac.uk

*/

//INCLUDE MODULES
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "TNT/tnt.h"

//SET NAME SPACE
using namespace std;
using namespace TNT;


void make_kernel(Array2D<float>& kernel, int t, int f, int h) {

  /*
  gaussian weigthed kernel must be pre-declared of size (2*f+1,2*f+1) and consist of zeros:
  Array2D<double> kernel(2*f+1, 2*f+1) = 0.0;
	*/

	double value;
	double bottom;
	double wgt = 0;

    for (int d=0; d<f; ++d) {

		bottom = 2.0*(d+1) + 1;
    value = 1.0 / bottom*bottom;

    for (int i=-d; i<d+1; ++i) {
        for (int j=-d; j<d+1; ++j) {
            kernel[f-i][f-j]= kernel[f-i][f-j] + value;
            wgt += value;
			}
		}
	}

	//may need to declare f as an array here for division:
	int kernel_size = 2*f+1;
	Array2D<float> denominator(kernel_size,kernel_size);
	for (int i=0; i<kernel_size; ++i) {
	   	for (int j=0; j<kernel_size; ++j) {
	  			kernel[i][j] = kernel[i][j]/wgt;
		}
	}
}



void make_gaussian_kernel(Array2D<float>& kernel, float sigma, int kernel_size) {

	/*

	Generate gaussian weighted kernel
	kernel array must be predeclared of size kernel_size and consist of zeros:
	Array2D<double> kernel(kernel_size, kernel_size) = 0.0;

	Kernel generated using:
  G(x,y) = (1/2*pi*sigma^2) exp ((-x^2+y^2)/(2*sigma^2))

  Martin Hurst, Feb 2012

  */

  float pi = 3.1416;
  double left_side = 1/(2*pi*sigma*sigma);
  double twosigma2 = 2.0*sigma*sigma;
  double right_side;
  double wgt = 0;
  double value;

	//calculate kernel values
  for (int i=0;i<2*kernel_size+1;++i) {
    for (int j=0;j<2*kernel_size+1;++j) {
    	right_side = -(((j-kernel_size)*(j-kernel_size) + (i-kernel_size)*(i-kernel_size))/twosigma2);
            //cout << i << " " << j << " " << right_side << endl;
      right_side = exp(right_side);
      value = left_side*right_side;
      kernel[i][j] = value;
      wgt += value;
    }
	}

	//scale to sum to 1
	for (int i=0;i<2*kernel_size+1;++i) {
    for (int j=0;j<2*kernel_size+1;++j) {
		  kernel[i][j] = kernel[i][j]/wgt;
		}
	}
}


void pad_array_symettric(Array2D<float>& array, Array2D<float>& padded_array, int nrows, int ncols, int pad_size) {

	/*

	Creates a buffer around an array (of size pad_size) and gives the new border
	mirror symmetric values of	the original array reflected across the boundary.
	pad_size should be the size of the window if filtering

	New array has size nrows + 2*pad_size x ncols + 2*pad_size

	Martin Hurst, Feb 2012

	*/

	int padded_nrows = nrows + 2*pad_size;
	int padded_ncols = ncols + 2*pad_size;

	int minus_i;
	int minus_j;

	for (int i=0; i<padded_nrows; ++i) {
		for (int j=0; j<padded_ncols; ++j) {

			//reverse of i and j
			minus_i = padded_nrows-1-i;
			minus_j = padded_ncols-1-j;

			//cout << i << endl;

			//top edge
			if (i<pad_size) {
				if (j<pad_size) {
					padded_array[i][j] = array[pad_size-i][pad_size-j];
				}
				else if (j>(ncols-1+pad_size)) {
					padded_array[i][j] = array[pad_size-i][j-pad_size-2*(pad_size-minus_j)];
					//cout << i << " " << j << "  " << padded_array[i][j] << endl;
				}
				else {
					padded_array[i][j] = array[pad_size-i][j-pad_size];
				}
			}
			//bottom edge
			else if (i>nrows-1+pad_size) {
				if (j<pad_size) {
					padded_array[i][j] = array[i-pad_size-2*(pad_size-minus_i)][pad_size-j];
				}
				else if (j>ncols+pad_size) {
					padded_array[i][j] = array[i-pad_size-2*(pad_size-minus_i)][j-pad_size-2*(pad_size-minus_j)];
				}
				else {
					padded_array[i][j] = array[i-pad_size-2*(pad_size-minus_i)][j-pad_size];
				}
			}
			//left side
			else if (j<pad_size) {
				padded_array[i][j] = array[i-pad_size][pad_size-j];
			}
			//right side
			else if (j>ncols-1+pad_size) {
				padded_array[i][j] = array[i-pad_size][j-pad_size-2*(pad_size-minus_j)];
			}
			//copy rest of array
			else {
				padded_array[i][j] = array[i-pad_size][j-pad_size];
			}
		}
	}
}


void NLocal_meansfilter(Array2D<float>& zeta, Array2D<float>& zeta_filt, Array2D<float>& noise, int nrows, int ncols, int t, int f, int h)
{
	/*

  zeta: elevation array
  t: radius of search window
  f: radius of similarity window
  h: degree of filtering

  t has to be <= f ?

  Adapted from a matlab script by:
  Author: Jose Vicente Manjon Herrera & Antoni Buades
  Date: 09-03-2006

  Implementation of the Non local filter proposed for A. Buades, B. Coll and J.M. Morel in
  "A non-local algorithm for image denoising"

  Martin Hurst, Feb 2012

  **Added soft threshold optimal correction - David Milodowski, 05/2012
  */

  //initiated arrays
  //Array2D<float> updated_zeta(nrows,ncols);
  Array2D<float> padded_zeta(nrows+2*f, ncols+2*f,0.0);
  //Created padded array with replicated boundaries
  pad_array_symettric(zeta, padded_zeta, nrows, ncols, f);

	//initiate kernel
	Array2D<float> kernel(2*f+1,2*f+1,0.0);
  //Generate kernel
  make_gaussian_kernel(kernel, 1.25, f);

  //initiate temporary arrays
  Array2D<float> W1(2*f+1,2*f+1);
  Array2D<float> W2(2*f+1,2*f+1);

	//initiate temp variables
	float w, wmax, average, sweight, d;

	//loop through DEM
	int i1, j1, rowmin, rowmax, colmin, colmax;

  for (int i=0; i<nrows; ++i) {

  i1 = i+f;
    for (int j=0; j<ncols; ++j) {
      j1 = j+f;
      //Get DEM sample  with size f centred on cell of interest
      for (int a=0; a<(2*f+1); ++a) {
  		  for (int b=0; b<(2*f+1); ++b) {
          W1[a][b] = padded_zeta[i1-f+a][j1-f+b];
				}
      }

      wmax=0;
      average=0;
  		sweight=0;

      //get bounding conditions
      rowmin = max(i1-t,f);
      rowmax = min(i1+t,nrows+f-1);
      colmin = max(j1-t,f);
      colmax = min(j1+t,ncols+f-1);

      //loop to calculate weigths for each cell
      for (int row=rowmin; row<rowmax+1; ++row) {
        for (int col=colmin; col<colmax+1; ++col) {
          d=0;
          //If centre cell do nothing
          if (row==i1 && col==j1) { }
          //Otherwise do stuff
          else {
  				  //Extract DEM centred around each point in kernel
            for (int a=0; a<(2*f+1); ++a) {
              for (int b=0; b<(2*f+1); ++b) {
						    W2[a][b] = padded_zeta[row+a-f][col+b-f];
	 					    d += kernel[a][b]*(W1[a][b]-W2[a][b])*(W1[a][b]-W2[a][b]);

              }
				    }

				    w = exp(-d/(h*h));

            if (w>wmax) wmax=w;
            sweight += w;
            average += w*padded_zeta[row][col];
          }
				}
      }

      average += wmax*padded_zeta[i1][j1];
      sweight += wmax;

      //if (sweight > 0) zeta_filt[i][j] = average/sweight;
      if (sweight > 0) {
        zeta_filt[i][j] = average/sweight;
      }
      else zeta_filt[i][j] = zeta[i][j];
		
		if (zeta_filt[i][j] < 0) zeta_filt[i][j] = -9999;
		
      // Also extract a record of the noise
      noise[i][j]=zeta[i][j]-zeta_filt[i][j];
    }
  }
}

