#include <PlasmaInjector.H>

using namespace amrex;

InjectorMomentum::~InjectorMomentum ()
{
    switch (type)
    {
    case InjectorMomentumType::parser:
    {
        object.parser.m_ux_parser.clear();
        object.parser.m_uy_parser.clear();
        object.parser.m_uz_parser.clear();
    }
    }
}

std::size_t
InjectorMomentum::sharedMemoryNeeded () const noexcept
{
    switch (type)
    {
    case InjectorMomentumType::parser:
    {
        return amrex::Gpu::numThreadsPerBlockParallelFor() * sizeof(double);
    }
    default:
        return 0;
    }
}

bool
InjectorMomentum::useRandom () const noexcept
{
    switch (type)
    {
    case InjectorMomentumType::gaussian:
    {
        return true;
    }
    default:
        return false;
    }
}

ConstantMomentumDistribution::ConstantMomentumDistribution(Real ux,
                                                           Real uy,
                                                           Real uz)
    : _ux(ux), _uy(uy), _uz(uz)
{}

void ConstantMomentumDistribution::getMomentum(vec3& u, Real x, Real y, Real z) {
    u[0] = _ux;
    u[1] = _uy;
    u[2] = _uz;
}

CustomMomentumDistribution::CustomMomentumDistribution(const std::string& species_name)
{
  ParmParse pp(species_name);
  pp.getarr("custom_momentum_params", params);
}

GaussianRandomMomentumDistribution::GaussianRandomMomentumDistribution(Real ux_m,
                                                                       Real uy_m,
                                                                       Real uz_m,
                                                                       Real ux_th,
                                                                       Real uy_th,
                                                                       Real uz_th)
    : _ux_m(ux_m), _uy_m(uy_m), _uz_m(uz_m), _ux_th(ux_th), _uy_th(uy_th), _uz_th(uz_th)
{
}

void GaussianRandomMomentumDistribution::getMomentum(vec3& u, Real x, Real y, Real z) {
    Real ux_th = amrex::RandomNormal(0.0, _ux_th);
    Real uy_th = amrex::RandomNormal(0.0, _uy_th);
    Real uz_th = amrex::RandomNormal(0.0, _uz_th);

    u[0] = _ux_m + ux_th;
    u[1] = _uy_m + uy_th;
    u[2] = _uz_m + uz_th;
}
RadialExpansionMomentumDistribution::RadialExpansionMomentumDistribution(Real u_over_r) : _u_over_r( u_over_r )
{
}

void RadialExpansionMomentumDistribution::getMomentum(vec3& u, Real x, Real y, Real z) {
  u[0] = _u_over_r * x;
  u[1] = _u_over_r * y;
  u[2] = _u_over_r * z;
}

ParseMomentumFunction::ParseMomentumFunction(std::string parse_momentum_function_ux,
                                             std::string parse_momentum_function_uy,
                                             std::string parse_momentum_function_uz)
    : _parse_momentum_function_ux(parse_momentum_function_ux),
      _parse_momentum_function_uy(parse_momentum_function_uy),
      _parse_momentum_function_uz(parse_momentum_function_uz)
{
    parser_ux.define(parse_momentum_function_ux);
    parser_uy.define(parse_momentum_function_uy);
    parser_uz.define(parse_momentum_function_uz);

    amrex::Array<std::reference_wrapper<WarpXParser>,3> parsers{parser_ux, parser_uy, parser_uz};
    ParmParse pp("my_constants");
    for (auto& p : parsers) {
        auto& parser = p.get();
        parser.registerVariables({"x","y","z"});
        std::set<std::string> symbols = parser.symbols();
        symbols.erase("x");
        symbols.erase("y");
        symbols.erase("z"); // after removing variables, we are left with constants
        for (auto it = symbols.begin(); it != symbols.end(); ) {
            Real v;
            if (pp.query(it->c_str(), v)) {
                parser.setConstant(*it, v);
                it = symbols.erase(it);
            } else {
                ++it;
            }
        }
        for (auto const& s : symbols) { // make sure there no unknown symbols
            amrex::Abort("ParseMomentumFunction: Unknown symbol "+s);
        }
    }
}

void ParseMomentumFunction::getMomentum(vec3& u, Real x, Real y, Real z)
{
    u[0] = parser_ux.eval(x,y,z);
    u[1] = parser_uy.eval(x,y,z);
    u[2] = parser_uz.eval(x,y,z);
}

