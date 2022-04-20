/******************************************************************************
 *
 *   Author:   Masa Prodanovic
 *   Copyright (c) 2009, The University of Texas at Austin. All rights reserved.
 *
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "segfl_util.h"

//utility functions for zipping/unzipping file
void checkUnzipFile(char *file_name,int *pzip_status, char **pfile_base)
{
     int     zip_status = (int)NO_ZIP, length_base;
     char    command[256], *file_base, *gptr, *bptr;
     
     file_base = (char *)malloc(256*sizeof(char));
     
     gptr = strstr(file_name,".gz");
     bptr = strstr(file_name,".bz2");
     if( gptr != (char *)NULL )
     {
       zip_status = GZIP;
       sprintf(command,"gunzip -f %s",file_name);
       system(command);
       //printf("\n%s",command);
       
       length_base = (strlen(file_name)-3);       
       file_base = strncpy(file_base,file_name,length_base);
       file_base[length_base] = '\0';
       //printf("\n%s",file_base); fflush(stdout); fflush(stderr);
     }
     else if( bptr!= (char *)NULL )
     {
       zip_status = BZIP2;
       sprintf(command,"bunzip2 -f %s",file_name);
       system(command);
       length_base =strlen(file_name)-4;
       file_base = strncpy(file_base,file_name,length_base);
       file_base[length_base] = '\0';
     }
     else 
       file_base = strcpy(file_base,file_name);
     
    *pzip_status = zip_status;
    *pfile_base = file_base;
}

void  zipFile(char *file_base,int zip_status)
{
     char command[256];
     
     if( zip_status == GZIP ) 
     {
        sprintf(command,"gzip -f %s",file_base);
        system(command);
     }
     else if( zip_status == BZIP2 ) 
     {
        sprintf(command,"bzip2 -f %s",file_base);
        system(command);
     }
}



//utility functions for reading segmented file

void	bit2unchar(
        unsigned char *bvect,
        unsigned char **evect,
	int           bnxy,
        int           resid,
        int           *pnxy)
{
	unsigned char *pb, *pe;
	unsigned char temp;

	int	nxy;
	int	i, j;

	nxy = bnxy*LENGTH + resid;	*pnxy = nxy;

	*evect = (unsigned char *)malloc(nxy*UCSZ);

	pe = *evect;	pb = bvect;

	for(  i = 0;  i < bnxy;  i++, pb++ )
	{
		temp = *pb;
		for( j = 0;  j < LENGTH-1;  j++, pe++ )
		{
			*pe = 0;
			*pe = (temp & 0200) ? 1 : 0;
			temp = temp << 1;
		}
		*pe = 0;
		*pe = (temp & 0200) ? 1 : 0;
		pe++;
	}

	temp = *pb;
	for( j = 0;  j < resid-1;  j++, pe++ )
	{
		*pe = 0;
		*pe = (temp & 0200) ? 1 : 0;
		temp = temp << 1;
	}
	*pe = 0;
	*pe = (temp & 0200) ? 1 : 0;

	temp = temp << LENGTH - resid;
}


void	bit_in_vol1(
		char		*fn,
		int		*nx,
		int		*ny,
		int		*zs,
		int		*ze,
		int		*bnxyz,
		int		*resid,
		unsigned char	**bvect)
{
	FILE	*fb;
	char	endian;
	int	nin, bnxyz1;


	if( (fb = fopen(fn,"r")) == NULL )
	{
		fprintf(stderr,"Unable to open file %s\n",fn);
	}

	fread(&endian,CSZ,1,fb);
	fread(nx     ,ISZ,1,fb);
	fread(ny     ,ISZ,1,fb);
	fread(zs     ,ISZ,1,fb);
	fread(ze     ,ISZ,1,fb);
	fread(bnxyz  ,ISZ,1,fb);
	fread(resid  ,ISZ,1,fb);

	bnxyz1 = (*bnxyz)+1;
	*bvect = (unsigned char *)malloc(bnxyz1*UCSZ);

	nin = fread(*bvect,UCSZ,bnxyz1,fb);

	if( nin != (*bnxyz)+1 )
	{
		fprintf(stderr,"Error in reading %s.\n",fn);
		fprintf(stderr,"incorrect item count\n");
	}
	fclose(fb);
}


void	bit2unchar_in_vol1(
		char		*fn,
		unsigned char	       **evect,
		int		*nx,
		int		*ny,
		int		*zs,
		int		*ze,
		int		*pbnxyz,
		int		*presid)
{
	unsigned char	*bvect;

	int	nxyz, bnxyz, resid, renxyz;

	bit_in_vol1(fn,nx,ny,zs,ze,&bnxyz,&resid,&bvect); /* BIT INPUT */

	bit2unchar(bvect,evect,bnxyz,resid,&renxyz);	 /* BIT TO UNCHAR */

	free(bvect);

	*pbnxyz = bnxyz;
	*presid = resid;

	nxyz = (*nx)*(*ny)*( (*ze)-(*zs)+1 );
	if( renxyz != nxyz )
	{
		printf("bit2unchar_in_vol() conversion error\n");
	}
}


void	set_uncompressed_filename(
		int		is_cmprss,
		char		*raw_fn,
		char		*uncmp_fn)
{
	char	cmnd[256];
	char	*after;

	if( is_cmprss )
	{
		after = strrchr(raw_fn,'/');
		if( after == NULL )
		{
			strcpy(uncmp_fn,"/tmp/");
			strcat(uncmp_fn,raw_fn);
		}
		else
		{
			strcpy(uncmp_fn,"/tmp");
			strcat(uncmp_fn,after);
		}

		sprintf(cmnd,UNZIP(raw_fn,uncmp_fn));
		system(cmnd);
	}
	else strcpy(uncmp_fn,raw_fn);
}

void	cleanup_tmp_files(
		int		in_cmprss,
		char		*uncmp_fn)
{
	char	cmnd[256];
	if( in_cmprss )
	{
		sprintf(cmnd,"/bin/rm -f %s",uncmp_fn);
		system(cmnd);
	}
}

//utility functions for writing a segmented file

void	unchar2bit(
	unsigned char	*cvect,
	unsigned char	**bvect,
	int		nxy,
	int		*pbnxy,
	int		*presid)
{
	unsigned char *pc, *pb;
	unsigned char temp;

	int	bnxy, resid;
	int	i, j;

	bnxy = nxy/LENGTH;	resid = nxy - bnxy*LENGTH;

	if( resid == 0 )  { bnxy--;	resid = LENGTH; }

	*pbnxy  = bnxy;		*presid = resid;

	*bvect = (unsigned char *)malloc((bnxy+1)*UCSZ);

	pc = cvect;	pb = *bvect;

	for(  i = 0;  i < bnxy;  i++, pb++ )
	{
		temp = 0;
		for( j = 0;  j < LENGTH-1;  j++, pc++ )
		{
			if( *pc == 1 ) temp |= 01;
			temp = temp << 1;
		}
		if( *pc == 1 ) temp |= 01;
		pc++;
		*pb = temp;
	}

	temp = 0;
	for( j = 0;  j < resid-1;  j++, pc++ )
	{
		if( *pc == 1 ) temp |= 01;
		temp = temp << 1;
	}
	if( *pc == 1 ) temp |= 01;
	temp = temp << (LENGTH - resid);

	*pb = temp;
}

void	bit_out_vol(
	char		*fn,
	int		nx,
	int		ny,
	int		zs,
	int		ze,
	int		bnxyz,
	int		resid,
	unsigned char	*bvect)
{

	FILE	*fb;
	int	nout;

	if( (fb = fopen(fn,"w")) == NULL )
	{
		fprintf(stderr,"Unable to open file %s\n",fn);
	}

	Write_endian(fb)
	fwrite(&nx   ,ISZ,1,fb);
	fwrite(&ny   ,ISZ,1,fb);
	fwrite(&zs   ,ISZ,1,fb);
	fwrite(&ze   ,ISZ,1,fb);
	fwrite(&bnxyz,ISZ,1,fb);
	fwrite(&resid,ISZ,1,fb);

	nout = fwrite(bvect,UCSZ,bnxyz+1,fb);
	if( nout != bnxyz+1 )
	{
		fprintf(stderr,"Error in writing %s\n",fn);
		fprintf(stderr,"incorrect item count\n");
	}
	fclose(fb);
}


void	unchar2bit_out_vol(
	char		*fn,
	unsigned char	*cvect,
	int		nx,
	int		ny,
	int		zs,
	int		ze)
{

	unsigned char	*bvect;

	int	bnxyz, resid, nz;

	nz = ze - zs + 1;
	unchar2bit(cvect,&bvect,nx*ny*nz,&bnxyz,&resid);    /* UNCHAR TO BIT */
	bit_out_vol(fn,nx,ny,zs,ze,bnxyz,resid,bvect);	    /* BIT OUTPUT */
	fprintf(stderr,"\nnxyz %d bnxyz %d resid %d\n",nx*ny*nz,bnxyz,resid);
}


/* reads binary image stored in segmented 3DMA volume file format 
*/
void    read_segfl(
        char           *segfl_name,
	unsigned char  **pdat,
	int            *nx,
	int            *ny,
	int            *nz)
{
   char	uncmp_fn[256], *zip_end = ".gz";
   int	rdzs, rdze, bnxyz, resid, zip;

   unsigned char  *udat;
  
   
   if( strstr(segfl_name,zip_end) ) zip = 1;
  
   set_uncompressed_filename(zip,segfl_name,uncmp_fn);
   bit2unchar_in_vol1(uncmp_fn,&udat,nx,ny,&rdzs,&rdze,&bnxyz,&resid);
   *nz = rdze - rdzs + 1;

   printf("\nnx = %d,  ny = %d,  nz =%d\n",*nx,*ny,*nz);

   *pdat  = udat;
}
 
/*
*	Bitpacked segmented volume file syntax
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



void write_segfl_ubc(
     unsigned char *seg_data,
     int *n,
     char *file_name,
     int zip_status)
{
   FILE *fp;
   
   fp = fopen(file_name,"w");   

   /* write grid dimensions */
   //fwrite(n, sizeof(int), 3, fp); 

   /* write data array */
   fwrite(seg_data, sizeof(unsigned char), n[0]*n[1]*n[2], fp);

   fclose(fp);
   zipFile(file_name,zip_status);
}

unsigned char    *read_segfl_ubc(
        char     *file_name,
	int      *grid_dims_ghostbox)
{
   FILE    *fp;
   int     zip_status;
   int     num_gridpts;
   unsigned char    *data = NULL;
   char             *file_base;
   
   checkUnzipFile(file_name,&zip_status,&file_base);
   
   fp = fopen(file_base,"r");

   if( fp != NULL)
   {
     /* read grid dimensions; if commented they need to be provided
        in grid_dims_ghostbox */
     //fread(grid_dims_ghostbox, sizeof(int), 3, fp); 
  
     /* allocate memory for data array */ 
     num_gridpts = grid_dims_ghostbox[0] * grid_dims_ghostbox[1]
                 * grid_dims_ghostbox[2];
     data = (unsigned char *) malloc(num_gridpts*sizeof(unsigned char));

     /* read data array */ 
     fread(data, sizeof(unsigned char), num_gridpts, fp);

     fclose(fp);
     
     zipFile(file_base,zip_status);
   }
   else
   {
      printf("\nCould not open file %s",file_name);
   }
   free(file_base);
   return data;
}
