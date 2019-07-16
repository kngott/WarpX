#ifndef INJECTOR_DENSITY_H_
#define INJECTOR_DENSITY_H_

#include <AMReX_Gpu.H>
#include <AMReX_Dim3.H>
#include <WarpXParser.H>

enum struct InjectorDensityType { constant, custom, predefined, parser };



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

private:
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
private:
    std::string _parse_density_function;
    WarpXParser parser_density;
};

#endif
