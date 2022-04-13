void MEM_START(void);
void MEM_STATS(void);

gpointer MALLOC(gulong);
gpointer CALLOC(gulong,gulong);
gpointer REALLOC(gpointer,gulong);
void  FREE(gpointer);
