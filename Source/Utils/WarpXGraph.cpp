#include <WarpX.H>
#include "Particles/MultiParticleContainer.H"
#include "ablastr/utils/Communication.H"
#include <AMReX_Graph.H>

std::string
WarpX::GraphFabName ()
{
    return g_fabname;
}

// Costs is based on fine, so default to fine.
std::string
WarpX::GraphFabName (int lev, PatchType pt)
{
    return (g_fabname + ((pt == PatchType::fine) ? "_f_" : "_c_" ) + std::to_string(lev));
}

void
WarpX::GraphSetup ()
{
    BL_PROFILE("GraphSetup()");
    graph.clear();

    const int nlevels = Efield_fp.size();
    g_temp.resize(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab* Ex = Efield_fp[lev][0].get();

        g_temp[lev] = std::make_unique<amrex::LayoutData<double>>(boxArray(lev), Ex->DistributionMap());
    }
}

void
WarpX::GraphSetup (const int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    BL_PROFILE("GraphSetup(lev, ba, dm)");
    graph.clear();

    if (lev >= g_temp.size()) {
        g_temp.resize(lev+1);
    }

    g_temp[lev] = std::make_unique<amrex::LayoutData<double>>(ba, dm);
}


void
WarpX::GraphClearTemps ()
{
    BL_PROFILE("GraphClearTemps()");
    const int nlevels = g_temp.size();

    for (int lev = 0; lev < nlevels; ++lev)
    {
        const auto iarr = g_temp[lev]->IndexArray();

        for (int i : iarr) {
            (*g_temp[lev])[i] = 0.0;
        }
    }
}

void
WarpX::GraphAddTemps (int lev, const std::string& fab_name, const std::string& wgt_name, const double scaling)
{
    // sizeof(amrex::Real) is expected size of Fab data. Change/make an input if needed.

    graph.addLayoutData(*g_temp[lev], fab_name, sizeof(amrex::Real),
                        wgt_name, scaling);
}

void
WarpX::GraphAddTemps (const std::string& fab_name, const std::string& wgt_name, const std::vector<double>& scaling)
{
    const int nlevels = g_temp.size();
    for (int lev = 0; lev < nlevels; ++lev)
    {
        GraphAddTemps(lev, fab_name + std::to_string(lev), wgt_name, scaling[lev]);
    }
}

void
WarpX::GraphAddTemps (const std::string& fab_name, const std::string& wgt_name, const double scaling)
{
    const int nlevels = g_temp.size();
    for (int lev = 0; lev < nlevels; ++lev)
    {
        GraphAddTemps(lev, fab_name, wgt_name, scaling);
    }
}


void
WarpX::GraphAddCellsandParticles ()
{
    BL_PROFILE("GraphAddCellsandParticles()");
    const int nlevels = g_temp.size();

    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab* Ex = Efield_fp[lev][0].get();
        std::string fab_name = GraphFabName(lev);

        GraphClearTemps();
        // =========================
        // Add cells
        for (amrex::MFIter mfi(*Ex, false); mfi.isValid(); ++mfi) {
            const amrex::Box& gbx = mfi.growntilebox();
            (*g_temp[lev])[mfi.index()] += gbx.numPts();
        }

        std::string weight_name = "ncells";
        graph.addLayoutData(*g_temp[lev], fab_name, sizeof(amrex::Real),
                            weight_name, costs_heuristic_cells_wt);

        // =========================
        // Add particles by species
        const auto& mypc_ref = GetInstance().GetPartContainer();
        const auto nSpecies = mypc_ref.nSpecies();

        for (int i_s = 0; i_s < nSpecies; ++i_s) {
            GraphClearTemps();

            auto & myspc = mypc_ref.GetParticleContainer(i_s);
            for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti) {
                (*g_temp[lev])[pti.index()] += pti.numParticles();
            }

            weight_name = "n" + mypc_ref.GetSpeciesNames()[i_s];
            graph.addLayoutData(*g_temp[lev], fab_name, sizeof(amrex::Real),
                            weight_name, costs_heuristic_particles_wt);
        }
    }
}

void
WarpX::GraphAddLoadBalance (const int lev, const bool do_load_balance, const amrex::Vector<int>& new_dm,
                            const amrex::Real currentEfficiency, const amrex::Real proposedEfficiency)
{
    BL_PROFILE("GraphAddLoadBalance()");
    std::string fab_name = GraphFabName(lev);
    std::string wgt_name = "costs";

    amrex::LayoutData<double> temp((g_temp[lev])->boxArray(), (g_temp[lev])->DistributionMap());

    for (int i=0; i<temp.local_size(); ++i) {
        temp.data()[i] = double(g_temp[lev]->data()[i]);
    }

    // Add costs, with scale equaling currentEfficiency
    graph.addLayoutData(temp, fab_name, sizeof(amrex::Real), wgt_name, double(currentEfficiency));

    // Add generated mapping, with scale equaling proposedEfficiency
    // Use the name to show whether load balance occurred.
    wgt_name = (do_load_balance) ? "new_dm" : "unused_dm";
    graph.addNodeWeight(fab_name, wgt_name, new_dm, proposedEfficiency, true);
}

void
WarpX::GraphAddFillBoundaryE(const amrex::IntVect& ng, const bool nodal_sync, const double scaling)
{
    BL_PROFILE("GraphAddFillBoundary()");
    // Long term alternative: optional name string. Named calls get added to Graph.
    // Generalize with an enum to select the MF to use?
    // i.e. WarpXData::E, WarpXData::B ...

    const int nlevels = WarpX::finest_level;

    for (int lev=0; lev <= nlevels; ++lev)
    {
        std::size_t buffer_size = (do_single_precision_comms) ?
                             sizeof(ablastr::utils::communication::comm_float_type) : 0;

        int nComp = Efield_fp[lev][0]->nComp();


        if (nodal_sync) {
            graph.addFillBoundaryAndSync("FB&S_Ex_f_" + std::to_string(lev),
                                          GraphFabName(lev, PatchType::fine), scaling,
                                          *Efield_fp[lev][0], 0, nComp,
                                          ng, Geom(lev).periodicity(),
                                          buffer_size);
            if (lev>0) {
                graph.addFillBoundaryAndSync("FB&S_Ex_c_" + std::to_string(lev),
                                              GraphFabName(lev, PatchType::coarse), scaling,
                                              *Efield_cp[lev][0], 0, nComp,
                                              ng, Geom(lev-1).periodicity(),
                                              buffer_size);
            }
        } else {
            graph.addFillBoundary("FB_Ex_f_" + std::to_string(lev),
                                  GraphFabName(lev, PatchType::fine), scaling,
                                  *Efield_fp[lev][0], 0, nComp, ng,
                                  Geom(lev).periodicity(), false, buffer_size);
            if (lev>0) {
                graph.addFillBoundary("FB_Ex_c_" + std::to_string(lev),
                                      GraphFabName(lev, PatchType::coarse), scaling,
                                      *Efield_cp[lev][0], 0, nComp, ng,
                                      Geom(lev-1).periodicity(), false, buffer_size);
            }
        }
    }
}

void
WarpX::GraphPrintandClear (const std::string& graph_name, const int precision = 5)
{
    BL_PROFILE("GraphPrintandClear()");
    graph.print_table(graph_name, precision);
    graph.clear();

    amrex::ParallelDescriptor::Barrier();
}
