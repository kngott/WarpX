#ifndef INJECTOR_DENSITY_H_
#define INJECTOR_DENSITY_H_

#include <AMReX_Gpu.H>
#include <AMReX_Dim3.H>
#include <GpuParser.H>

enum struct InjectorDensityType { constant, custom, predefined, parser };

struct InjectorDensityConstant
{
    InjectorDensityConstant (amrex::Real a_rho) noexcept : m_rho(a_rho) {}

    AMREX_GPU_DEVICE
    amrex::Real
    getDensity (amrex::Real, amrex::Real, amrex::Real) const noexcept
    {
        return m_rho;
    }

private:
    amrex::Real m_rho;
};

struct InjectorDensityParser
{
    InjectorDensityParser (WarpXParser const& a_parser) noexcept
        : m_parser(a_parser) {}

    AMREX_GPU_DEVICE
    amrex::Real
    getDensity (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        return m_parser(x,y,z);
    }

    GpuParser m_parser;
};

struct InjectorDensity
    : public amrex::Gpu::Managed
{
    InjectorDensity (InjectorDensityType a_type, amrex::Real a_rho)
        : object(a_rho),
          type(a_type)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_type == InjectorDensityType::constant,
                                         "InjectorDensity: must be type constant");
    }

    InjectorDensity (InjectorDensityType a_type, WarpXParser const& a_parser)
        : object(a_parser),
          type(a_type)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_type == InjectorDensityType::parser,
                                         "InjectorDensity: must be type parser");
    }

    ~InjectorDensity ();

    std::size_t sharedMemoryNeeded () const noexcept;

    AMREX_GPU_DEVICE
    amrex::Real
    getDensity (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        switch (type)
        {
        case InjectorDensityType::parser:
        {
            return object.parser.getDensity(x,y,z);
        }
        case InjectorDensityType::constant:
        {
            return object.constant.getDensity(x,y,z);
        }
        default:
        {
            amrex::Abort("InjectorDensity: unknown type");
            return 0.0;
        }
        }
    }

private:
    union Object {
        Object (amrex::Real a_rho) noexcept : constant(a_rho) {}
        Object (WarpXParser const& a_parser) noexcept : parser(a_parser) {}
        InjectorDensityConstant constant;
        InjectorDensityParser   parser;
    };
    Object object;
    InjectorDensityType type;
};


enum class predefined_profile_flag { Null, parabolic_channel };

///
/// PlasmaDensityProfile describes how the charge density
/// is set in particle initialization. Subclasses must define a
/// getDensity function that describes the charge density as a
/// function of x, y, and z.
///
class PlasmaDensityProfile
{
public:
    virtual ~PlasmaDensityProfile() {};
    virtual amrex::Real getDensity(amrex::Real x,
                                   amrex::Real y,
                                   amrex::Real z) const = 0;
protected:
    std::string _species_name;
};

///
/// This describes a constant density distribution.
///
class ConstantDensityProfile : public PlasmaDensityProfile
{
public:
    ConstantDensityProfile(amrex::Real _density);
    virtual amrex::Real getDensity(amrex::Real x,
                                   amrex::Real y,
                                   amrex::Real z) const override;

//private:
    amrex::Real _density;
};

///
/// This describes a custom density distribution. Users can supply
/// in their problem directory.
///
///
class CustomDensityProfile : public PlasmaDensityProfile
{
public:
    CustomDensityProfile(const std::string& species_name);
    virtual amrex::Real getDensity(amrex::Real x,
                                   amrex::Real y,
                                   amrex::Real z) const override;
private:
    amrex::Vector<amrex::Real> params;
};

///
/// This describes predefined density distributions.
///
class PredefinedDensityProfile : public PlasmaDensityProfile
{
public:
    PredefinedDensityProfile(const std::string& species_name);
    virtual amrex::Real getDensity(amrex::Real x,
                                   amrex::Real y,
                                   amrex::Real z) const override;
    amrex::Real ParabolicChannel(amrex::Real x,
                                 amrex::Real y,
                                 amrex::Real z) const;
private:
    predefined_profile_flag which_profile = predefined_profile_flag::Null;
    amrex::Vector<amrex::Real> params;
};

///
/// This describes a density function parsed in the input file. 
///
class ParseDensityProfile : public PlasmaDensityProfile
{
public:
    ParseDensityProfile(const std::string _parse_density_function);
    virtual amrex::Real getDensity(amrex::Real x,
                                   amrex::Real y,
                                   amrex::Real z) const override;
//private:
    std::string _parse_density_function;
    WarpXParser parser_density;
};

#endif
