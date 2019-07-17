#include "PlasmaInjector.H"

#include <sstream>
#include <functional>

#include <WarpXConst.H>
#include <WarpX_f.H>
#include <AMReX.H>
#include <WarpX.H>

using namespace amrex;

namespace {
    void StringParseAbortMessage(const std::string& var,
                                 const std::string& name) {
        std::stringstream stringstream;
        std::string string;
        stringstream << var << " string '" << name << "' not recognized.";
        string = stringstream.str();
        amrex::Abort(string.c_str());
    }

    Real parseChargeName(const ParmParse pp, const std::string& name) {
        Real result;
        if (name == "q_e") {
            return PhysConst::q_e;
        } else if (pp.query("charge", result)) {
            return result;
        } else {
            StringParseAbortMessage("Charge", name);
            return 0.0;
        }
    }

    Real parseChargeString(const ParmParse pp, const std::string& name) {
        if(name.substr(0, 1) == "-")
            return -1.0 * parseChargeName(pp, name.substr(1, name.size() - 1));
        return parseChargeName(pp, name);
    }

    Real parseMassString(const ParmParse pp, const std::string& name) {
        Real result;
        if (name == "m_e") {
            return PhysConst::m_e;
        } else if (name == "m_p"){
            return PhysConst::m_p;
        } else if (name == "inf"){
	    return std::numeric_limits<double>::infinity();
        } else if (pp.query("mass", result)) {
            return result;
        } else {
            StringParseAbortMessage("Mass", name);
            return 0.0;
        }
    }
}

PlasmaInjector::PlasmaInjector(){
    part_pos = NULL;
}

PlasmaInjector::PlasmaInjector(int ispecies, const std::string& name)
    : species_id(ispecies), species_name(name)
{
    ParmParse pp(species_name);

    // parse charge and mass
    std::string charge_s;
    pp.get("charge", charge_s);
    std::transform(charge_s.begin(),
                   charge_s.end(),
                   charge_s.begin(),
                   ::tolower);
    charge = parseChargeString(pp, charge_s);

    std::string mass_s;
    pp.get("mass", mass_s);
    std::transform(mass_s.begin(),
                   mass_s.end(),
                   mass_s.begin(),
                   ::tolower);
    mass = parseMassString(pp, mass_s);

    // parse injection style
    std::string part_pos_s;
    pp.get("injection_style", part_pos_s);
    std::transform(part_pos_s.begin(),
                   part_pos_s.end(),
                   part_pos_s.begin(),
                   ::tolower);
    if (part_pos_s == "python") {
        return;
    } else if (part_pos_s == "singleparticle") {
        pp.getarr("single_particle_pos", single_particle_pos, 0, 3);
        pp.getarr("single_particle_vel", single_particle_vel, 0, 3);
        for (auto& x : single_particle_vel) {
            x *= PhysConst::c;
        }
        pp.get("single_particle_weight", single_particle_weight);
        add_single_particle = true;
        return;
    } else if (part_pos_s == "gaussian_beam") {
        pp.get("x_m", x_m);
        pp.get("y_m", y_m);
        pp.get("z_m", z_m);
        pp.get("x_rms", x_rms);
        pp.get("y_rms", y_rms);
        pp.get("z_rms", z_rms);
        pp.get("q_tot", q_tot);
        pp.get("npart", npart);
        pp.query("do_symmetrize", do_symmetrize);
        gaussian_beam = true;
        parseMomentum(pp);
    }
    else if (part_pos_s == "nrandompercell") {
        pp.query("num_particles_per_cell", num_particles_per_cell);
        part_pos.reset(new RandomPosition(num_particles_per_cell));
        parseDensity(pp);
        parseMomentum(pp);
    } else if (part_pos_s == "nuniformpercell") {
        num_particles_per_cell_each_dim.resize(3);
        pp.getarr("num_particles_per_cell_each_dim", num_particles_per_cell_each_dim);
#if ( AMREX_SPACEDIM == 2 )
        num_particles_per_cell_each_dim[2] = 1;
#endif
        part_pos.reset(new RegularPosition(num_particles_per_cell_each_dim));
        num_particles_per_cell = num_particles_per_cell_each_dim[0] *
                                 num_particles_per_cell_each_dim[1] *
                                 num_particles_per_cell_each_dim[2];
        parseDensity(pp);
        parseMomentum(pp);
    } else {
        StringParseAbortMessage("Injection style", part_pos_s);
    }

    pp.query("radially_weighted", radially_weighted);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(radially_weighted, "ERROR: Only radially_weighted=true is supported");

    // parse plasma boundaries
    xmin = std::numeric_limits<amrex::Real>::lowest();
    ymin = std::numeric_limits<amrex::Real>::lowest();
    zmin = std::numeric_limits<amrex::Real>::lowest();

    xmax = std::numeric_limits<amrex::Real>::max();
    ymax = std::numeric_limits<amrex::Real>::max();
    zmax = std::numeric_limits<amrex::Real>::max();

    pp.query("xmin", xmin);
    pp.query("ymin", ymin);
    pp.query("zmin", zmin);
    pp.query("xmax", xmax);
    pp.query("ymax", ymax);
    pp.query("zmax", zmax);

}

void PlasmaInjector::parseDensity(ParmParse pp){
    // parse density information
    std::string rho_prof_s;
    pp.get("profile", rho_prof_s);
    std::transform(rho_prof_s.begin(),
                   rho_prof_s.end(),
                   rho_prof_s.begin(),
                   ::tolower);
    if (rho_prof_s == "constant") {
        pp.get("density", density);
        rho_prof.reset(new ConstantDensityProfile(density));
    } else if (rho_prof_s == "custom") {
        rho_prof.reset(new CustomDensityProfile(species_name));
    } else if (rho_prof_s == "predefined") {
        rho_prof.reset(new PredefinedDensityProfile(species_name));
    } else if (rho_prof_s == "parse_density_function") {
        pp.get("density_function(x,y,z)", str_density_function);
        rho_prof.reset(new ParseDensityProfile(str_density_function));
    } else {
        StringParseAbortMessage("Density profile type", rho_prof_s);
    }
}

void PlasmaInjector::parseMomentum(ParmParse pp){
    // parse momentum information
    std::string mom_dist_s;
    pp.get("momentum_distribution_type", mom_dist_s);
    std::transform(mom_dist_s.begin(),
                   mom_dist_s.end(),
                   mom_dist_s.begin(),
                   ::tolower);
    if (mom_dist_s == "constant") {
        Real ux = 0.;
        Real uy = 0.;
        Real uz = 0.;
        pp.query("ux", ux);
        pp.query("uy", uy);
        pp.query("uz", uz);
        mom_dist.reset(new ConstantMomentumDistribution(ux, uy, uz));
    } else if (mom_dist_s == "custom") {
        mom_dist.reset(new CustomMomentumDistribution(species_name));
    } else if (mom_dist_s == "gaussian") {
        Real ux_m = 0.;
        Real uy_m = 0.;
        Real uz_m = 0.;
        Real ux_th = 0.;
        Real uy_th = 0.;
        Real uz_th = 0.;
        pp.query("ux_m", ux_m);
        pp.query("uy_m", uy_m);
        pp.query("uz_m", uz_m);
        pp.query("ux_th", ux_th);
        pp.query("uy_th", uy_th);
        pp.query("uz_th", uz_th);
        mom_dist.reset(new GaussianRandomMomentumDistribution(ux_m, uy_m, uz_m, 
                                                              ux_th, uy_th, uz_th));
    } else if (mom_dist_s == "radial_expansion") {
        Real u_over_r = 0.;
        pp.query("u_over_r", u_over_r);
        mom_dist.reset(new RadialExpansionMomentumDistribution(u_over_r));
    } else if (mom_dist_s == "parse_momentum_function") {
        pp.get("momentum_function_ux(x,y,z)", str_momentum_function_ux);
        pp.get("momentum_function_uy(x,y,z)", str_momentum_function_uy);
        pp.get("momentum_function_uz(x,y,z)", str_momentum_function_uz);
        mom_dist.reset(new ParseMomentumFunction(str_momentum_function_ux, 
                                                 str_momentum_function_uy, 
                                                 str_momentum_function_uz));
    } else {
        StringParseAbortMessage("Momentum distribution type", mom_dist_s);
    }
}

void PlasmaInjector::getPositionUnitBox(vec3& r, int i_part, int ref_fac) {
    return part_pos->getPositionUnitBox(r, i_part, ref_fac);
}

void PlasmaInjector::getMomentum(vec3& u, Real x, Real y, Real z) {
    mom_dist->getMomentum(u, x, y, z);
    u[0] *= PhysConst::c;
    u[1] *= PhysConst::c;
    u[2] *= PhysConst::c;
}

bool PlasmaInjector::insideBounds(Real x, Real y, Real z) {
  if (x >= xmax || x < xmin ||
      y >= ymax || y < ymin ||
      z >= zmax || z < zmin ) return false;
  return true;
}

Real PlasmaInjector::getDensity(Real x, Real y, Real z) {
    return rho_prof->getDensity(x, y, z);
}

InjectorPosition*
PlasmaInjector::getInjectorPosition ()
{
    if (inj_pos == nullptr)
    {
        if (RandomPosition const* ppp =
            dynamic_cast<RandomPosition const*>(part_pos.get()))
        {
            inj_pos.reset(new InjectorPosition(InjectorPositionType::random,
                                               xmin, xmax, ymin, ymax, zmin, zmax));
        }
        else if (RegularPosition const* ppp =
                 dynamic_cast<RegularPosition const*>(part_pos.get()))
        {
            inj_pos.reset(new InjectorPosition(InjectorPositionType::regular,
                                               xmin, xmax, ymin, ymax, zmin, zmax,
                                               Dim3{ppp->_num_particles_per_cell_each_dim[0],
                                                    ppp->_num_particles_per_cell_each_dim[1],
#if (AMREX_SPACEDIM == 3)
                                                    ppp->_num_particles_per_cell_each_dim[2]
#else
                                                    1.0;
#endif
                                               }));
        }
        else
        {
            amrex::Abort("PlasmaInjector::getInjectorPosition: how did this happend");
        }
    }
    return inj_pos.get();
}

InjectorDensity*
PlasmaInjector::getInjectorDensity ()
{
    if(inj_rho == nullptr)
    {
        if (ConstantDensityProfile const* p =
            dynamic_cast<ConstantDensityProfile const*>(rho_prof.get()))
        {
            inj_rho.reset(new InjectorDensity(InjectorDensityType::constant,
                                              p->_density));
        }
        else if (ParseDensityProfile const* p =
                 dynamic_cast<ParseDensityProfile const*>(rho_prof.get()))
        {
            inj_rho.reset(new InjectorDensity(InjectorDensityType::parser,
                                              p->parser_density));
        }
        else
        {
            amrex::Abort("PlasmaInjector::getInjectorDensity: how did this happen?");
        }
    }
    return inj_rho.get();
}

InjectorMomentum*
PlasmaInjector::getInjectorMomentum ()
{
    if (inj_mom == nullptr)
    {
        if (ConstantMomentumDistribution const* p =
            dynamic_cast<ConstantMomentumDistribution const*>(mom_dist.get()))
        {
            inj_mom.reset(new InjectorMomentum(InjectorMomentumType::constant,
                                              p->_ux, p->_uy, p->_uz));
        }
        else if (ParseMomentumFunction const* p =
                 dynamic_cast<ParseMomentumFunction const*>(mom_dist.get()))
        {
            inj_mom.reset(new InjectorMomentum(InjectorMomentumType::parser,
                                               p->parser_ux, p->parser_uy,
                                               p->parser_uz));
        }
        else
        {
            amrex::Abort("PlasmaInjector::getInjectorMomentum: how did this happen?");
        }
    }
    return inj_mom.get();
}

