#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>

#include "segfl_util.h"
#include "mem.h"

//utility functions for reading segmented file

void	bit2unchar(
        guchar *bvect,
        guchar **evect,
	gint           bnxy,
        gint           resid,
        gint           *pnxy)
{
	guchar *pb, *pe;
	guchar temp;

	gint	nxy;
	gint	i, j;

	nxy = bnxy*LENGTH + resid;	*pnxy = nxy;

	*evect = (guchar *)MALLOC(nxy*UCSZ);

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
		gchar		*fn,
		gint		*nx,
		gint		*ny,
		gint		*zs,
		gint		*ze,
		gint		*bnxyz,
		gint		*resid,
		guchar	**bvect)
{
	FILE	*fb;
	gchar	endian;
	gint	nin, bnxyz1;


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
	*bvect = (guchar *)MALLOC(bnxyz1*UCSZ);

	nin = fread(*bvect,UCSZ,bnxyz1,fb);

	if( nin != (*bnxyz)+1 )
	{
		fprintf(stderr,"Error in reading %s.\n",fn);
		fprintf(stderr,"incorrect item count\n");
	}
	fclose(fb);
}


void	bit2unchar_in_vol1(
		gchar		*fn,
		guchar	**evect,
		gint		*nx,
		gint		*ny,
		gint		*zs,
		gint		*ze,
		gint		*pbnxyz,
		gint		*presid)
{
	guchar	*bvect;

	gint	nxyz, bnxyz, resid, renxyz;

	bit_in_vol1(fn,nx,ny,zs,ze,&bnxyz,&resid,&bvect); /* BIT INPUT */

	bit2unchar(bvect,evect,bnxyz,resid,&renxyz);	 /* BIT TO UNCHAR */

	FREE(bvect);

	*pbnxyz = bnxyz;
	*presid = resid;

	nxyz = (*nx)*(*ny)*( (*ze)-(*zs)+1 );
	if( renxyz != nxyz )
	{
		printf("bit2unchar_in_vol() conversion error\n");
	}
}


void	set_uncompressed_filename(
		gint		is_cmprss,
		gchar		*raw_fn,
		gchar		*uncmp_fn)
{
	gchar	cmnd[256];
	gchar	*after;

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
		gint		in_cmprss,
		gchar		*uncmp_fn)
{
	gchar	cmnd[256];
	if( in_cmprss )
	{
		sprintf(cmnd,"/bin/rm -f %s",uncmp_fn);
		system(cmnd);
	}
}

void	read_segfl(
	gchar           *segfl_base,
        gint              zip,
	Cube_Info	*cube,
	guchar	**dat)
{
	guchar	*pdat;
        gchar	uncmp_fn[256];
	gint	rdx, rdy, rdz, rdzs, rdze, bnxyz, resid;

        set_uncompressed_filename(zip,segfl_base,uncmp_fn);
	bit2unchar_in_vol1(uncmp_fn,&pdat,&rdx,&rdy,&rdzs,&rdze,&bnxyz,
                                                                       &resid);
	rdz = rdze - rdzs + 1;

        cube->nx = rdx;  cube->ny = rdy;  cube->nz = rdz;
        cube->nxy = rdx * rdy;            cube->nxyz = rdx*rdy*rdz;
	cube->zs = rdzs; cube->ze = rdze;
        cube->xs = cube->ys = 0;
        cube->xe = rdx-1;    cube->ye = rdy-1;

	*dat = pdat;
        cleanup_tmp_files(zip,uncmp_fn);
}



void    read_segfl_top(
        gchar *segfl_name,
	guchar **pdat,
	Cube_Info *pcube )
{
   gchar       *segfl_base;
   gboolean    zip, file_exists;

   Cube_Info   cube;
   guchar  *dat; 
   // FILE   *fp;
   
   file_exists = g_file_test(segfl_name,G_FILE_TEST_EXISTS);

   if( file_exists )
   {
       zip = g_str_has_suffix(segfl_name, ".gz");       //test if file is zipped

       if( zip )
       { 
	   segfl_base = g_strndup(segfl_name, strlen(segfl_name)-3);
       }
       else segfl_base = g_strndup(segfl_name, strlen(segfl_name));

       read_segfl(segfl_base,zip,&cube,&dat); //change into read size info only

       g_message("x = %d..%d,  y = %d..%d,  z =%d..%d",
		     cube.xs, cube.xe, cube.ys, cube.ye,cube.zs, cube.ze);

       *pcube = cube;
       *pdat  = dat;
   }
   else
   {
       g_error("\nFile %s doesn't exist.",segfl_name);

   }
   //print temporary segfile data onto disc
   //fp = fopen("/tmp/tmp_segfl","w");
   //fwrite(&cube,1,sizeof(Cube_Info),fp);
   //fwrite(dat,cube.nxyz,UCSZ,fp);
   //fclose(fp);
}
 
