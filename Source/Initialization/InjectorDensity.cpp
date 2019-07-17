#include <PlasmaInjector.H>

using namespace amrex;

InjectorDensity::~InjectorDensity ()
{
    switch (type)
    {
    case InjectorDensityType::parser:
    {
        object.parser.m_parser.clear();
    }
    }
}

std::size_t
InjectorDensity::sharedMemoryNeeded () const noexcept
{
    switch (type)
    {
    case InjectorDensityType::parser:
    {
        return amrex::Gpu::numThreadsPerBlockParallelFor() * sizeof(double);
    }
    default:
        return 0;
    }
}

ConstantDensityProfile::ConstantDensityProfile(Real density)
    : _density(density)
{}

Real ConstantDensityProfile::getDensity(Real x, Real y, Real z) const
{
    return _density;
}

CustomDensityProfile::CustomDensityProfile(const std::string& species_name)
{
    ParmParse pp(species_name);
    pp.getarr("custom_profile_params", params);
}

PredefinedDensityProfile::PredefinedDensityProfile(const std::string& species_name)
{
    ParmParse pp(species_name);
    std::string which_profile_s;
    pp.getarr("predefined_profile_params", params);
    pp.query("predefined_profile_name", which_profile_s);
    if (which_profile_s == "parabolic_channel"){
        which_profile = predefined_profile_flag::parabolic_channel;
    }
}

ParseDensityProfile::ParseDensityProfile(std::string parse_density_function)
    : _parse_density_function(parse_density_function)
{
    parser_density.define(parse_density_function);
    parser_density.registerVariables({"x","y","z"});

    ParmParse pp("my_constants");
    std::set<std::string> symbols = parser_density.symbols();
    symbols.erase("x");
    symbols.erase("y");
    symbols.erase("z"); // after removing variables, we are left with constants
    for (auto it = symbols.begin(); it != symbols.end(); ) {
        Real v;
        if (pp.query(it->c_str(), v)) {
            parser_density.setConstant(*it, v);
            it = symbols.erase(it);
        } else {
            ++it;
        }
    }
    for (auto const& s : symbols) { // make sure there no unknown symbols
        amrex::Abort("ParseDensityProfile: Unknown symbol "+s);
    }
}

Real ParseDensityProfile::getDensity(Real x, Real y, Real z) const
{
    return parser_density.eval(x,y,z);
}

