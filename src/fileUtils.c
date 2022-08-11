/** An interface to missing file operations in Fortran **/

#include <stdlib.h>
#include <stdio.h>

#ifdef ENABLE_GZ
#include <zlib.h>
#define MYFILE_T gzFile
#define fopen(fname, fmode) gzopen(fname, fmode)
#define fread(the_ptr, the_size, the_count, the_file)                          \
  (gzread(the_file, the_ptr, (the_size * the_count)) / (the_size))
#define fseek(the_file, the_value, the_mode)                                   \
  gzseek(the_file, the_value, the_mode)
#define fclose(the_file) gzclose(the_file)
#else
#define MYFILE_T FILE *
#endif


int64_t openIt(char *name) {
  MYFILE_T fp;
  fp = fopen(name, "rb");
  return (int64_t) fp;
}

void seekIt(int64_t fp, int64_t pos) {
  (void) fseek((MYFILE_T ) fp, pos, SEEK_SET);
}
  
void closeIt(int64_t fp) {
  fclose((MYFILE_T) fp);
}

void readIt(void *ptr, int64_t fp, int64_t offset, int64_t isize) {
  seekIt(fp, offset);
  fread(ptr, 1, isize, (MYFILE_T )fp);
}

