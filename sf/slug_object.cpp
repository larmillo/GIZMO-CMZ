#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <slug.H>
#include <slug_cluster.H>
#include <slug_predefined.H>
#include "slug_object.h"

void slug_object::construct_cluster(double particle_mass) {
  cluster = new slug_cluster(1, particle_mass, 0.0, 
                             slug_predef.imf("chabrier", 9.0, POISSON), 
                             slug_predef.tracks("modp020.dat"), 
                             slug_predef.specsyn("sb99",
                                                 slug_predef.tracks("modp020.dat"),
                                                 slug_predef.imf("chabrier", 9.0, POISSON)), 
                             slug_predef.filter_set("QH0"), 
                             nullptr, nullptr, 
                             slug_predef.yields("SNII_Sukhbold16_nodecay"),
							 nullptr,
                             slug_predef.ostreams, nullptr, true);

}

// Method to reconstruct the slug_cluster object from a serialized buffer
void slug_object::reconstruct_cluster(char *buf) {

  cluster = new slug_cluster((slug_cluster_buffer *) buf,
                             slug_predef.imf("chabrier", 9.0, POISSON),
                             slug_predef.tracks("modp020.dat"),
                             slug_predef.specsyn("sb99",
                                                 slug_predef.tracks("modp020.dat"),
                                                 slug_predef.imf("chabrier", 9.0, POISSON)),
                             slug_predef.filter_set("QH0"), 
                             nullptr, nullptr, 
                             slug_predef.yields("SNII_Sukhbold16_nodecay"),nullptr,
                             slug_predef.ostreams, nullptr, true);

}

// Method to return the member functions of the slug_cluster
void slug_object::pack_buffer(char *buf) { cluster->pack_buffer((slug_cluster_buffer *) buf); }
int slug_object::buffer_size() { return cluster->buffer_size(); }
void slug_object::advance(double particle_age) { cluster->advance(particle_age); }
int slug_object::get_stoch_sn() { return cluster->get_stoch_sn(); }
double slug_object::get_birth_mass() { return cluster->get_birth_mass(); }
double slug_object::get_stellar_mass() { return cluster->get_stellar_mass(); }
double slug_object::get_photometry_QH0() { return cluster->get_photometry()[0]; }


extern "C" {
	slug_object *slug_object_new() {return new slug_object();}
	void   slug_object_delete(slug_object *S) {delete S;}
	void   slug_construct_cluster(slug_object *S, double particle_mass) {return S->construct_cluster(particle_mass);}	
    void   slug_reconstruct_cluster(slug_object *S, char *buf) {return S->reconstruct_cluster(buf);}
    void   slug_pack_buffer(slug_object *S, char *buf) {return S->pack_buffer(buf);}
    int    slug_buffer_size(slug_object *S) {return S->buffer_size();}
    void   slug_advance(slug_object *S, double particle_age) { return S->advance(particle_age);}
    int    slug_get_stoch_sn(slug_object *S){return S->get_stoch_sn();}
	double slug_get_birth_mass(slug_object *S){return S->get_birth_mass();}
	double slug_get_stellar_mass(slug_object *S){return S->get_stellar_mass();}
    double slug_get_photometry_QH0(slug_object *S){return S->get_photometry_QH0();}
}	