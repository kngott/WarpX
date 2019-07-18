#ifndef INJECTOR_MOMENTUM_H_
#define INJECTOR_MOMENTUM_H_

#include <AMReX_Gpu.H>
#include <AMReX_Dim3.H>
#include <GpuParser.H>

enum struct InjectorMomentumType { constant, custom, gaussian, RadialExpansion, parser };

struct InjectorMomentumConstant
{
    InjectorMomentumConstant (amrex::Real a_ux, amrex::Real a_uy, amrex::Real a_uz) noexcept
        : m_ux(a_ux), m_uy(a_uy), m_uz(a_uz) {}

    AMREX_GPU_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real, amrex::Real, amrex::Real) const noexcept
    {
        return amrex::XDim3{m_ux,m_uy,m_uz};
    }
private:
    amrex::Real m_ux, m_uy, m_uz;
};

struct InjectorMomentumGaussian
{
    InjectorMomentumGaussian (amrex::Real a_ux_m, amrex::Real a_uy_m,
                              amrex::Real a_uz_m, amrex::Real a_ux_th,
                              amrex::Real a_uy_th, amrex::Real a_uz_th) noexcept
        : m_ux_m(a_ux_m), m_uy_m(a_uy_m), m_uz_m(a_uz_m),
          m_ux_th(a_ux_th), m_uy_th(a_uy_th), m_uz_th(a_uz_th)
        {}

    AMREX_GPU_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        return amrex::XDim3{amrex::RandomNormal(m_ux_m, m_ux_th),
                            amrex::RandomNormal(m_uy_m, m_uy_th),
                            amrex::RandomNormal(m_uz_m, m_uz_th)};
    }
private:
    amrex::Real m_ux_m, m_uy_m, m_uz_m;
    amrex::Real m_ux_th, m_uy_th, m_uz_th;
};

struct InjectorMomentumParser
{
    InjectorMomentumParser (WarpXParser const& a_ux_parser,
                            WarpXParser const& a_uy_parser,
                            WarpXParser const& a_uz_parser) noexcept
        : m_ux_parser(a_ux_parser), m_uy_parser(a_uy_parser),
          m_uz_parser(a_uz_parser) {}

    AMREX_GPU_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        return amrex::XDim3{m_ux_parser(x,y,z),m_uy_parser(x,y,z),m_uz_parser(x,y,z)};
    }

    GpuParser m_ux_parser, m_uy_parser, m_uz_parser;
};

struct InjectorMomentum
    : public amrex::Gpu::Managed
{
    InjectorMomentum (InjectorMomentumType a_type, amrex::Real a_ux,
                      amrex::Real a_uy, amrex::Real a_uz)
        : object(a_ux, a_uy, a_uz),
          type(a_type)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_type == InjectorMomentumType::constant,
                                         "InjectorMomentum: must be type constant");
    }

    InjectorMomentum (InjectorMomentumType a_type, WarpXParser const& a_ux_parser,
                      WarpXParser const& a_uy_parser, WarpXParser const& a_uz_parser)
        : object(a_ux_parser, a_uy_parser, a_uz_parser),
          type(a_type)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_type == InjectorMomentumType::parser,
                                         "InjectorMomentum: must be type parser");
    }

    InjectorMomentum (InjectorMomentumType a_type, amrex::Real a_ux_m,
                      amrex::Real a_uy_m, amrex::Real a_uz_m, amrex::Real a_ux_th,
                      amrex::Real a_uy_th, amrex::Real a_uz_th)
        : object(a_ux_m,a_uy_m,a_uz_m,a_ux_th,a_uy_th,a_uz_th),
          type(a_type)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_type == InjectorMomentumType::gaussian,
                                         "InjectorMomentum: must be type gaussian");
    }

    ~InjectorMomentum ();

    std::size_t sharedMemoryNeeded () const noexcept;

    bool useRandom () const noexcept;

    AMREX_GPU_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        switch (type)
        {
        case InjectorMomentumType::parser:
        {
            return object.parser.getMomentum(x,y,z);
        }
        case InjectorMomentumType::gaussian:
        {
            return object.gaussian.getMomentum(x,y,z);
        }
        case InjectorMomentumType::constant:
        {
            return object.constant.getMomentum(x,y,z);
        }
        default:
        {
            amrex::Abort("InjectorMomentum: unknown type");
            return {0.0,0.0,0.0};
        }
        }
    }

private:
    union Object {
        Object (amrex::Real a_ux, amrex::Real a_uy, amrex::Real a_uz) noexcept
            : constant(a_ux,a_uy,a_uz) {}
        Object (amrex::Real a_ux_m, amrex::Real a_uy_m,
                amrex::Real a_uz_m, amrex::Real a_ux_th,
                amrex::Real a_uy_th, amrex::Real a_uz_th) noexcept
            : gaussian(a_ux_m,a_uy_m,a_uz_m,a_ux_th,a_uy_th,a_uz_th) {}
        Object (WarpXParser const& a_ux_parser, WarpXParser const& a_uy_parser,
                WarpXParser const& a_uz_parser) noexcept
            : parser(a_ux_parser, a_uy_parser, a_uz_parser) {}
        InjectorMomentumConstant constant;
        InjectorMomentumGaussian gaussian;
        InjectorMomentumParser   parser;
    };
    Object object;
    InjectorMomentumType type;
};


///
/// PlasmaMomentumDistribution describes how the particle momenta
/// are set. Subclasses must define a getMomentum method that fills
/// a u with the 3 components of the particle momentum
///
class PlasmaMomentumDistribution
{
public:
    using vec3 = std::array<amrex::Real, 3>;
    virtual ~PlasmaMomentumDistribution() {};
    virtual void getMomentum(vec3& u, amrex::Real x, amrex::Real y, amrex::Real z) = 0;
};

///
/// This is a constant momentum distribution - all particles will
/// have the same ux, uy, and uz
///
class ConstantMomentumDistribution : public PlasmaMomentumDistribution
{
public:
    ConstantMomentumDistribution(amrex::Real ux,
                                 amrex::Real uy,
                                 amrex::Real uz);
    virtual void getMomentum(vec3& u, amrex::Real x, amrex::Real y, amrex::Real z) override;

//private:
    amrex::Real _ux;
    amrex::Real _uy;
    amrex::Real _uz;
};

///
/// This describes a custom momentum distribution. Users can supply
/// in their problem directory.
///
///
class CustomMomentumDistribution : public PlasmaMomentumDistribution
{
public:
    CustomMomentumDistribution(const std::string& species_name);
    virtual void getMomentum(vec3& u, amrex::Real x, amrex::Real y, amrex::Real z) override;

private:
    amrex::Vector<amrex::Real> params;
};


///
/// This is a Gaussian Random momentum distribution.
/// Particles will get random momenta, drawn from a normal.
/// ux_m, ux_y, and ux_z describe the mean components in the x, y, and z
/// directions, while u_th is the standard deviation of the random
/// component.
///
class GaussianRandomMomentumDistribution : public PlasmaMomentumDistribution
{
public:
    GaussianRandomMomentumDistribution(amrex::Real ux_m,
                                       amrex::Real uy_m,
                                       amrex::Real uz_m,
                                       amrex::Real ux_th,
                                       amrex::Real uy_th,
                                       amrex::Real uz_th);
    virtual void getMomentum(vec3& u, amrex::Real x, amrex::Real y, amrex::Real z) override;
//private:
    amrex::Real _ux_m;
    amrex::Real _uy_m;
    amrex::Real _uz_m;
    amrex::Real _ux_th;
    amrex::Real _uy_th;
    amrex::Real _uz_th;
};

///
/// This is a radially expanding momentum distribution
/// Particles will have a radial momentum proportional to their 
/// radius, with proportionality constant u_over_r
class RadialExpansionMomentumDistribution : public PlasmaMomentumDistribution
{
public:
  RadialExpansionMomentumDistribution( amrex::Real u_over_r );
  virtual void getMomentum(vec3& u, amrex::Real x, amrex::Real y, amrex::Real z) override;
private:
    amrex::Real _u_over_r;
};

///
/// This describes a momentum distribution function parsed in the input file. 
///
class ParseMomentumFunction : public PlasmaMomentumDistribution
{
public:
    ParseMomentumFunction(const std::string _parse_momentum_function_ux,
                          const std::string _parse_momentum_function_uy,
                          const std::string _parse_momentum_function_uz);
    virtual void getMomentum(vec3& u, 
                             amrex::Real x,
                             amrex::Real y,
                             amrex::Real z) override;
//private:
    std::string _parse_momentum_function_ux;
    std::string _parse_momentum_function_uy;
    std::string _parse_momentum_function_uz;
    WarpXParser parser_ux;
    WarpXParser parser_uy;
    WarpXParser parser_uz;
};

#endif
