#include <PlasmaInjector.H>

using namespace amrex;

bool
InjectorPosition::useRandom () const noexcept
{
    switch (type)
    {
    case InjectorPositionType::random:
    {
        return true;
    }
    default:
    {
        return false;
    }
    };
}

RandomPosition::RandomPosition(int num_particles_per_cell):
  _num_particles_per_cell(num_particles_per_cell)
{}

void RandomPosition::getPositionUnitBox(vec3& r, int i_part, int ref_fac){
    r[0] = amrex::Random();
    r[1] = amrex::Random();
    r[2] = amrex::Random();
}

RegularPosition::RegularPosition(const amrex::Vector<int>& num_particles_per_cell_each_dim)
    : _num_particles_per_cell_each_dim(num_particles_per_cell_each_dim)
{}

void RegularPosition::getPositionUnitBox(vec3& r, int i_part, int ref_fac)
{
  int nx = ref_fac*_num_particles_per_cell_each_dim[0];
  int ny = ref_fac*_num_particles_per_cell_each_dim[1];
#if AMREX_SPACEDIM == 3
  int nz = ref_fac*_num_particles_per_cell_each_dim[2];
#else
  int nz = 1;
#endif
  
  int ix_part = i_part/(ny * nz);
  int iy_part = (i_part % (ny * nz)) % ny;
  int iz_part = (i_part % (ny * nz)) / ny;

  r[0] = (0.5+ix_part)/nx;
  r[1] = (0.5+iy_part)/ny;
  r[2] = (0.5+iz_part)/nz;
}
