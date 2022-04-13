#ifndef INCLUDE_SEG_UTIL_H
#define INCLUDE_SEG_UTIL_H


/*
*	Copyrighted, Research Foundation of SUNY, 2005
*/

//from 3DMA
typedef struct _Cube_Info
{
	gint	xs, xe, ys, ye, zs, ze;
	gint	nx, ny, nz, nxy, nxyz;
} Cube_Info;

#define LENGTH 8
#define CSZ   sizeof(gchar)
#define ISZ   sizeof(gint)
#define UCSZ  sizeof(guchar)

#define UNZIP(raw_fn,uncmp_fn) "gunzip -cf %s > %s",raw_fn,uncmp_fn

//utility functions for reading 3DMA segmented file
void    read_segfl_top(gchar *, guchar **,Cube_Info *);
void	read_segfl(gchar *,gint,Cube_Info *,guchar **);

void	bit2unchar(guchar *,guchar **,gint,gint,gint *);
void	bit2unchar_in_vol1(gchar *,guchar **,gint *,gint *,gint *,gint *,
			                                        gint *,gint *);
void	bit_in_vol1(gchar *,gint *,gint *,gint *,gint *,gint *,gint *,
                                                                    guchar **);
void	bit2unchar(guchar *,guchar **,gint,gint,gint *);
void	set_uncompressed_filename(gint, gchar *,gchar *);
void	cleanup_tmp_files(gint,gchar *);

#endif
