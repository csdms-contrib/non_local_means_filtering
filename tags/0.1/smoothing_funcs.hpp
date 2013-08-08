/*

smoothing_funcs.h provides functions to perform non-local means filtering following Baudes et al (2005)

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

#ifndef SMOOTH_FUNC
#define SMOOTH_FUNC

//INCLUDE MODULES
#include <iostream>
#include <fstream>
#include <iomanip>
#include "math.h"
#include "TNT/tnt.h"

//SET NAME SPACE
using namespace std;
using namespace TNT;


void make_kernel(Array2D<float>& kernel, int t, int f, int h);

    /*
    gaussian weigthed kernel must be pre-declared of size (2*f+1,2*f+1) and consist of zeros:
    Array2D<double> kernel(2*f+1, 2*f+1) = 0.0;
	*/


void make_gaussian_kernel(Array2D<float>& kernel, float sigma, int kernel_size);

	/*

	Generate gaussian weighted kernel
	kernel array must be predeclared of size kernel_size and consist of zeros:
	Array2D<double> kernel(kernel_size, kernel_size) = 0.0;

	Kernel generated using:
    G(x,y) = (1/2*pi*sigma^2) exp ((-x^2+y^2)/(2*sigma^2))

    Martin Hurst, Feb 2012

    */


void pad_array_symettric(Array2D<float> array, Array2D<float>& padded_array, int nrows, int ncols, int pad_size);

	/*

	Creates a buffer around an array (of size pad_size) and gives the new border
	mirror symmetric values of	the original array reflected across the boundary.
	pad_size should be the size of the window if filtering

	New array has size nrows + 2*pad_size x ncols + 2*pad_size

	Martin Hurst, Feb 2012

	*/


void NLocal_meansfilter(Array2D<float>& zeta, Array2D<float>& zeta_filt, Array2D<float>& noise, int nrows, int ncols, int t, int f, int h);

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

    */

#endif // #ifndef SMOOTH_FUNC
