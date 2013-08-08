/*

NL_means_filter.cpp driver file to perform non-local means filtering following Baudes et al (2005)
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
#include <cstring>
#include "math.h"
#include "TNT/tnt.h"
#include "read_dem_v1.3.hpp"
#include "smoothing_funcs.hpp"

//SET NAME SPACE
using namespace std;
using namespace TNT;

//------------------------------------------------------------------------------
//
// Perform Non-local means filtering on a DEM following Baude et al. [2005]
//
// Inputs required: input dem file name, output smoothed_dem file name, the
// search window radius, the similarity window and the degree of filtering
//
// Martin Hurst, February, 2012
// Modified by David Milodowski, May 2012- generates grid of recording filtered
// noise
//
//------------------------------------------------------------------------------

int main(int nNumberofArgs,char *argv[]) {

	//Check arguments
	if (nNumberofArgs!=5) {
		cout << "FATAL ERROR: You must provide four input arguments:\n\t"
				" - an input file descriptor\n\t - the search window radius (try 2) \n\t"
				" - the similarity window radius (try 2) \n\t - the degree of filtering (try 2)\n" << endl;
		exit(EXIT_FAILURE);
	}

	//get the file names as arguments
	char *dem_in = argv[1];
	int t = atoi(argv[2]);
	int f = atoi(argv[3]);
	int h = atoi(argv[4]);

	//generate file names
	int len = strlen(dem_in);
	char *dem_header = new char[len+8];
	char *dem_file = new char[len+8];
	char *filt_file = new char[len+9];
	char *filt_header = new char[len+9];
	char *noise_file = new char[len+10];
	char *noise_header = new char[len+10];

	strcpy(dem_header,dem_in);
	strcat(dem_header,"_dem.hdr");
	strcpy(dem_file,dem_in);
	strcat(dem_file,"_dem.flt");
	strcpy(filt_header,dem_in);
	strcat(filt_header,"_filt.hdr");
	strcpy(filt_file,dem_in);
	strcat(filt_file,"_filt.flt");
	strcpy(noise_header,dem_in);
	strcat(noise_header,"_noise.hdr");
	strcpy(noise_file,dem_in);
	strcat(noise_file,"_noise.flt");

	//Declare parameters
	int ncols, nrows;
	double xmin, ymin;
	float dem_res;
	int ndv;

	//READ IN DATA
	//read in header
	read_raster_hdr(ncols, nrows, xmin, ymin, dem_res, ndv, dem_header);

	//Declare Arrays
	Array2D<float> zeta(nrows,ncols);
	Array2D<float> zeta_filt(nrows,ncols);
	Array2D<float> noise(nrows,ncols);
	//read in arrays
	read_raster_flt(ncols, nrows, zeta, dem_file);

	cout << "\n*** Non-local means filtering: smooths a DEM ***\n" << endl;

	Array2D<float> kernel(2*f+1,2*f+1,0.0);
	NLocal_meansfilter(zeta, zeta_filt, noise, nrows, ncols, t, f, h);

	cout << "\n\tWriting filtered DEM to " << filt_file << endl;
	write_raster_hdr(ncols,nrows,xmin,ymin,dem_res,ndv,filt_header);
	write_raster_flt(ncols,nrows,zeta_filt,filt_file);
	cout << "\n\tWriting noise to " << noise_file << endl;
	write_raster_hdr(ncols,nrows,xmin,ymin,dem_res,ndv,noise_header);
	write_raster_flt(ncols,nrows,noise,noise_file);
	cout << "\tDone!\n" << endl;

}

