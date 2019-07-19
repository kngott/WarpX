#ifndef INJECTOR_POSITION_H_
#define INJECTOR_POSITION_H_

#include <AMReX_Gpu.H>
#include <AMReX_Dim3.H>
#include <AMReX_Utility.H>

struct InjectorPositionRandom
{
    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getPositionUnitBox (int i_part, int ref_fac=1) const noexcept
    {
        return amrex::XDim3{amrex::Random(), amrex::Random(), amrex::Random()};
    }
};

struct InjectorPositionRegular
{
    InjectorPositionRegular (amrex::Dim3 const& a_ppc) noexcept : ppc(a_ppc) {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getPositionUnitBox (int i_part, int ref_fac=1) const noexcept
    {
        int nx = ref_fac*ppc.x;
        int ny = ref_fac*ppc.y;
#if (AMREX_SPACEDIM == 3)
        int nz = ref_fac*ppc.z;
#else
        int nz = 1;
#endif
        int ix_part = i_part/(ny*nz);  // written this way backward compatibility
        int iz_part = (i_part-ix_part*(ny*nz)) / ny;
        int iy_part = (i_part-ix_part*(ny*nz)) - ny*iz_part;
        return amrex::XDim3{(0.5+ix_part)/nx, (0.5+iy_part)/ny, (0.5+iz_part) / nz};
    }
private:
    amrex::Dim3 ppc;
};

struct InjectorPosition
    : public amrex::Gpu::Managed
{
    InjectorPosition (InjectorPositionRandom* t,
                      amrex::Real a_xmin, amrex::Real a_xmax,
                      amrex::Real a_ymin, amrex::Real a_ymax,
                      amrex::Real a_zmin, amrex::Real a_zmax)
        : type(Type::random),
          object(t),
          xmin(a_xmin), xmax(a_xmax),
          ymin(a_ymin), ymax(a_ymax),
          zmin(a_zmin), zmax(a_zmax)
    { }

    InjectorPosition (InjectorPositionRegular* t,
                      amrex::Real a_xmin, amrex::Real a_xmax,
                      amrex::Real a_ymin, amrex::Real a_ymax,
                      amrex::Real a_zmin, amrex::Real a_zmax,
                      amrex::Dim3 const& a_ppc)
        : type(Type::regular),
          object(t, a_ppc),
          xmin(a_xmin), xmax(a_xmax),
          ymin(a_ymin), ymax(a_ymax),
          zmin(a_zmin), zmax(a_zmax)
    { }

    InjectorPosition (InjectorPosition const&) = delete;
    InjectorPosition (InjectorPosition&&) = delete;
    void operator= (InjectorPosition const&) = delete;
    void operator= (InjectorPosition &&) = delete;

    std::size_t sharedMemoryNeeded () const noexcept { return 0; }

    bool useRandom () const noexcept;

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getPositionUnitBox (int i_part, int ref_fac=1) const noexcept
    {
        switch (type)
        {
        case Type::regular:
        {
            return object.regular.getPositionUnitBox(i_part, ref_fac);
        }
        default:
        {
            return object.random.getPositionUnitBox(i_part, ref_fac);
        }
        };
    }

    AMREX_GPU_HOST_DEVICE
    bool
    insideBounds (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        return (x < xmax and x >= xmin and
                y < ymax and y >= ymin and
                z < zmax and z >= zmin);
    }

private:
    enum struct Type { random, regular };
    Type type;

    union Object {
        Object (InjectorPositionRandom*) noexcept : random() {}
        Object (InjectorPositionRegular*, amrex::Dim3 const& a_ppc) noexcept
            : regular(a_ppc) {}
        InjectorPositionRandom random;
        InjectorPositionRegular regular;
    };
    Object object;

    amrex::Real xmin, xmax;
    amrex::Real ymin, ymax;
    amrex::Real zmin, zmax;
};

///
/// PlasmaParticlePosition describes how particles are initialized
/// into each cell box. Subclasses must define a
/// getPositionUnitBox function that returns the position of
/// particle number i_part in a unitary box.
///
class PlasmaParticlePosition{
public:
  using vec3 = std::array<amrex::Real, 3>;
  virtual ~PlasmaParticlePosition() {};
    virtual void getPositionUnitBox(vec3& r, int i_part, int ref_fac=1) = 0;
};

///
/// Particles are initialized with a random uniform
/// distribution inside each cell
///
class RandomPosition : public PlasmaParticlePosition{
public:
    RandomPosition(int num_particles_per_cell);
    virtual void getPositionUnitBox(vec3& r, int i_part, int ref_fac=1) override;
//private:
    amrex::Real _x;
    amrex::Real _y;
    amrex::Real _z;
    int _num_particles_per_cell;
};

///
/// Particles are regularly distributed inside each cell. The user provides
/// a 3d (resp. 2d) vector num_particles_per_cell_each_dim that contains
/// the number of particles per cell along each dimension.
///
class RegularPosition : public PlasmaParticlePosition{
public:
  RegularPosition(const amrex::Vector<int>& num_particles_per_cell_each_dim);
    virtual void getPositionUnitBox(vec3& r, int i_part, int ref_fac=1) override;
//private:
  amrex::Real _x;
  amrex::Real _y;
  amrex::Real _z;
  amrex::Vector<int> _num_particles_per_cell_each_dim;
};

#endif
