#ifndef INCLUDE_DOUBLE_DATA_UTIL_H
#define INCLUDE_DOUBLE_DATA_UTIL_H

void read_gdouble_array(gdouble **p_phi,gint *n,gchar *fname,gint);
void read_gfloat_array_into_gdouble(gdouble **p_phi,gint *n,gchar *fname,gint);
void write_gdouble_array(gdouble *phi,gint *n,gchar *fname);

#endif
