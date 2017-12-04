#include <stdio.h>
void auto_show_flags(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"CPP = /usr/bin/cpp\n");
   fprintf (fp,"CC  = /usr/bin/gcc\n");
   fprintf (fp,"CXX = /usr/bin/g++\n");
   fprintf (fp,"FC  = /usr/local/bin/gfortran\n");
   fprintf (fp,"F90 = /usr/local/bin/gfortran\n");
   fprintf (fp,"LD  = /usr/bin/gcc\n");
   fprintf (fp,"LIBTOOL = /usr/local/bin/glibtool\n");
   fprintf (fp,"\n");
   fprintf (fp,"DEFINES = -DLINUX -DH5_USE_16_API   -DCONFIG_BFLOAT_8\n");
   fprintf (fp,"\n");
   fprintf (fp,"INCLUDES = -I/usr/local/include    -I.\n");
   fprintf (fp,"\n");
   fprintf (fp,"CPPFLAGS = -P -traditional \n");
   fprintf (fp,"CFLAGS   =  -O2 \n");
   fprintf (fp,"CXXFLAGS =  -O2 \n");
   fprintf (fp,"FFLAGS   = -fno-second-underscore -m64 -O2 \n");
   fprintf (fp,"F90FLAGS = -fno-second-underscore -m64 -O2 \n");
   fprintf (fp,"LDFLAGS  =  \n");
   fprintf (fp,"\n");
   fprintf (fp,"LIBS     = -L/usr/local/lib -lhdf5 -L/usr/local/lib/gcc/7 -lgfortran  \n");
   fprintf (fp,"\n");
}
