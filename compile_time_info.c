#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        HYDRO_MESHLESS_FINITE_MASS\n"
"        EOS_GAMMA=(5.0/3.0)\n"
"        COOLING\n"
"        GRACKLE\n"
"        GRACKLE_CHEMISTRY=0\n"
"        GRACKLE_OPTS\n"
"        MULTIPLEDOMAINS=64\n"
"        ADAPTIVE_GRAVSOFT_FORGAS\n"
"        HAVE_HDF5\n"
"        OUTPUT_ADDITIONAL_RUNINFO\n"
"        SLUG\n"
"        STAR_FORMATION\n"
"\n");
}
