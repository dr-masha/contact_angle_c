/******************************************************************************
 *
 *   Author:   Masa Prodanovic
 *   Copyright (c) 2009, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/

#include <arpa/inet.h>

#define LENGTH 8
#define CSZ   sizeof(char)
#define ISZ   sizeof(int)
#define UCSZ  sizeof(unsigned char)

#define   ZIP(fn) "gzip -nf %s",fn
#define UNZIP(raw_fn,uncmp_fn) "gunzip -cf %s > %s",raw_fn,uncmp_fn

#define NO_ZIP 0
#define GZIP   1
#define BZIP2  2


#define Test_endian(test)\
{\
	int	i = 1;\
	test = ( i == htonl(i) ) ? 'E' : 'e';\
}

#define Write_endian(fb)\
{\
	char test;\
\
	Test_endian(test)\
	fwrite(&test,CSZ,1,fb);\
}


/*! 
 * checkUnzipFile() checks if file has .gz or .bz2 extention and 
 *   uncompresses the file.
 *
 *   Arguments:
 *    - file_name(in):     name of the file
 *    - pzip_status(out):  integer pointer that points to status of the
 *                         file (NO_ZIP, GZIP, BZIP2)
 *    - pfile_base (out):  pointer to string containing the file name base 
 *                         (file_name of the uncompressed file, same as the 
 *                         file name if there was no compression to begin with).
 *
 *  Return value:   none
 *
 *  Notes: 
 *     - There is no attempt to open the file.
 *     - file_base is allocated to 256 bytes within the function, and needs 
 *       to be freed later on
 *
 */        
void checkUnzipFile(char *file_name,int *pzip_status,char **pfile_base);

/*! zipFile() compresses the file according to its status and frees
 *   file_base pointer.
 *
 *   Arguments:
 *    - file_base(in):     name of the uncompressed file
 *    - zip_status(in):    integer compression status of the file 
 *                        (NO_ZIP,GZIP,BZIP2)
 *
 *  Return value:   none
 * 
 *  Notes: In case the status is NO_ZIP, function doesn't do anything.
 */        
void   zipFile(char *file_base,int zip_status);




/*!	Bitpacked segmented volume file syntax
*	endian  - binary character, 'e': file written on little endian hardware
*				    'E': file written on  big   endian hardware
*	nx
*	ny
*	zs      - starting z-slice (convention is to start from 1)
*	ze      -  ending  z-slice
*	b_nxyz  - number of full bytes to store the volume
*	n_resid - number of bits used in the last byte of storage required
*	b_nxyz + 1 unsigned char containing the bit-packed segmentation
*		   information for the volume
*/


/* utility functions for reading bitpacked segmented (volume) file
   Each voxel information is stored in one bit (NOT byte) of data
  */
void	bit2unchar(unsigned char *,unsigned char **,int,int,int *);
void	bit2unchar_in_vol1(char *,unsigned char **,int *,int *,int *,int *,
			                                        int *,int *);
void	bit_in_vol1(char *,int *,int *,int *,int *,int *,int *,unsigned char **);
void	bit2unchar(unsigned char *,unsigned char **,int,int,int *);
void	set_uncompressed_filename(int, char *,char *);
void	cleanup_tmp_files(int,char *);

/* utility functions for writing bitpacked segmented (volume) file */
void	bit_out_vol(char *,int,int,int,int,int,int,unsigned char *);
void	unchar2bit(unsigned char *,unsigned char **,int,int *,int *);
void	unchar2bit_out_vol(char *,unsigned char *,int,int,int,int);


/*!	Unsigned array segmented volume file syntax
	nx      - binary integer (number of voxels in x dir)
	ny      - binary integer (number of voxels in y dir)
	nz      - binary integer (number of voxels in z dir)
        array of nx*ny*nz unsigned characters (0 and 1)
        
	Packing convention:
         for k=1:nz
          for j=1:ny
	    for i=1:nx
	       store (i,j,k) position information 
*/

/* read/write functions for unsigned char segmented (volume) arrays
   Note: unlike the segmented files above, these use 1 byte of data
   per voxel stored */
   
void write_segfl_ubc(
     unsigned char *seg_data,
     int *n,
     char *file_name,
     int zip_status);
     
unsigned char *read_segfl_ubc(
     char     *file_name,
     int      *grid_dims_ghostbox);
   
