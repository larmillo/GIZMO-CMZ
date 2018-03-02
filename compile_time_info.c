#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        HYDRO_MESHLESS_FINITE_MASS\n"
"        EOS_GAMMA=(5.0/3.0)\n"
"        NOGRAVITY\n"
"        ANALYTIC_GRAVITY\n"
"        HAVE_HDF5\n"
"        OUTPUT_ADDITIONAL_RUNINFO\n"
"\n");
}
