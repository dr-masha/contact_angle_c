/******************************************************************************
 *
 *   Author:   Masa Prodanovic
 *   Copyright (c) 2009, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <./segfl_util.h>

/*  Usage:
    convert_ubc_3dma fname_in_ubc nx ny nz fname_out_3dma
    
    fname_in_ubc - filename (including the path to it) of the unsigned char 
                   binary array; it can be zipped (.gz)
    nx,ny,nz     - number of voxels in x,y,z directions
                   The data is assumed to be packed as follows
		   packing_index = k*nx*ny + j*ny + i
		   (store data z-plane by z-plane; within each plane,
		    y direction is stored first then x
		    This is C/C++ convention.
		    Fortran and Matlab store data as 		   
		   packing_index = k*nx*ny + i*nx + j
    fname_out_3dma - filename for 3DMA file (it will be zipped  by default)	
*/

int main(int argc, char **argv) 
{
      char fname_in[256], fname_out[256], command[256];
      int  n[3], zs, ze;
      unsigned char *data;

      if( argc >= 6)
      {
        sprintf(fname_in,"%s",argv[1]);
        n[0] = atoi(argv[2]);
        n[1] = atoi(argv[3]);
        n[2] = atoi(argv[4]);
        sprintf(fname_out,"%s",argv[5]);
      }
      else
      {
        printf("\nUsage: convert_ubc_3dma fname_in_ubc nx ny nz fname_out_3dma\n");	
	return 1;
      }
      
      /** read ubc binary file **/
      data = read_segfl_ubc(fname_in,n);
      
      /*** write 3DMA data ***/
      zs = 1;
      ze = n[2];
      unchar2bit_out_vol(fname_out,data,n[0],n[1],zs,ze);
      zipFile(fname_out,GZIP);
      
      free(data);
      
      return 0;
}
