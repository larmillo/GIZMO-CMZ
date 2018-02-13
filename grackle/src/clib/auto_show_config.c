#include <stdio.h>
void auto_show_config(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"   MACHINE: Use apt-get to install csh libhdf5-serial-dev gfortran\n");
   fprintf (fp,"   MACHINE-NAME: avatar\n");
   fprintf (fp,"\n");
   fprintf (fp,"   CONFIG_PRECISION  [precision-{32,64}]                     : 64\n");
   fprintf (fp,"   CONFIG_OPT  [opt-{warn,debug,high,aggressive}]            : high\n");
   fprintf (fp,"   CONFIG_OMP  [omp-{on,off}]                                : off\n");
   fprintf (fp,"\n");
}
