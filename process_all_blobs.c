#include <stdio.h>
#include <stdlib.h>

/***
'process_all_blobs basename phase1 phase2 blob_id_list out_fname'

basename - blob fluid segmented files basename w/o blob id or fluid
           identifier; files assumed gzipped
           -- basically a basename given to blob segfiles in case7.14,
              3DMA-Rock
blob_id_list - files storing id's of blobs to process
phase1,phase2 - phase identifiers e.g. OIL,H2O
         - if processing oil blobs put OIL as first identifier and vice versa
out_fname - basename for areas/contact angles filename
***/

int main(int argc, char * argv[])
{
    char *basename, *phase1, *phase2, *blob_id_list, *out_fname;
    char fname1[256], fname2[256], cmnd[256];
    int  n, i, num, *blob_id;
    FILE *fp;

    basename = argv[1];
    phase1 = argv[2];
    phase2 = argv[3];
    blob_id_list = argv[4];
    out_fname = argv[5];
    
    //fprintf(stderr,"\n%s %s %s %s %s",argv[1],argv[2],argv[3],argv[4],argv[5]); 
    //fflush(stdout);

    //read blob id's from blob_id_list
    fp = fopen(blob_id_list,"r");
    n = 0;
    while( fscanf(fp,"%d",&num) != EOF ) n++; //count number of elements

    blob_id =(int *)calloc(n,sizeof(int));

    n = 0;  rewind(fp);
    while( fscanf(fp,"%d",blob_id+n) != EOF ) n++; //store numbers
    fclose(fp);

    //run contact angles/areas code for all blobs
    for(i=0; i < n;i++)
    {
       fprintf(stderr,"\nBlob %d",blob_id[i]);
       sprintf(fname1,"%s_%d_%s.gz",basename,blob_id[i],phase1);
       sprintf(fname2,"%s_%d_%s.gz",basename,blob_id[i],phase2);
       sprintf(cmnd,"contact_angle %s %s %s %d", fname1,fname2,out_fname,blob_id[i]);
       printf("\n%s\n",cmnd); fflush(stdout);
       system(cmnd);
    }
}
