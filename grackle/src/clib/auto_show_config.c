#include <stdio.h>
void auto_show_config(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"   MACHINE: Darwin (OSX) \n");
   fprintf (fp,"   MACHINE-NAME: darwin\n");
   fprintf (fp,"\n");
   fprintf (fp,"   CONFIG_PRECISION  [precision-{32,64}]                     : 64\n");
   fprintf (fp,"   CONFIG_OPT  [opt-{warn,debug,high,aggressive}]            : high\n");
   fprintf (fp,"\n");
}
