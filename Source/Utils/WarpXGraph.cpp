#include "AMReX_Graph.H"
#include "Utils/WarpXGraph.H"

std::string
WarpX::GraphFabName()
{
    return ("Ex_lev");
}

std::string
WarpX::GraphFabName(int lev)
{
    return ("Ex_lev" + std::to_string(lev));
}

void
WarpX::GraphResetTemps()
{
    for (int lev = 0; lev <= finest_level(); ++lev)
    {
        const auto iarr = g_temp[lev]->IndexArray();
        for (int i : iarr) {
            (*g_temp[lev])(i) = 0.0;
        }
    }
}

void
WarpX::GraphAddTemps(int lev, const std::string& fab_name, const std::string& wgt_name, const double scaling)
{
    graph.addLayoutData(*g_temp[lev], fab_name + std::to_string(lev), sizeof(amrex::Real),
                        wgt_name, scaling);
}

void
WarpX::GraphAddTemps(const std::string& fab_name, const std::string& wgt_name, const std::vector<double>& scaling)
{
    for (int lev = 0; lev <= finest_level(); ++lev)
    {
        GraphAddTemps(lev, fab_name, wgt_name, scaling[lev]);
    }
}

void
WarpX::GraphAddTemps(const std::string& fab_name, const std::string& wgt_name, const double& scaling)
{
    for (int lev = 0; lev <= finest_level(); ++lev)
    {
        GraphAddTemps(lev, fab_name, wgt_name, scaling);
    }
}


void
WarpX::GraphAddCellsandParticles()
{
    for (int lev = 0; lev <= finest_level(); ++lev) {

        amrex::MultiFab* Ex = Efield_fp[lev][0].get();
        std::string fab_name = GraphFabName(lev);

        g_temp[lev] = std::make_unique<LayoutData<Real>>(boxArray(lev), Ex.DistributionMapping());

        GraphResetTemps();
        // =========================
        // Add cells
        for (MFIter mfi(*Ex, false); mfi.isValid(); ++mfi) {
            const Box& gbx = mfi.growntilebox();
            (*g_temp[lev])[mfi.index()] += gbx.numPts();
        }

        std::string weight_name = "ncells";
        graph.addLayoutData(*g_temp[lev], fab_name, sizeof(amrex::Real),
                            weight_name, costs_heuristic_cells_wt);

        // =========================
        // Add particles by species
        const auto& mypc_ref = GetInstance().GetPartContainer();
        const auto nSpecies = mypc_ref.nSpecies();

        for (int i_s = 0; i_s < nSpecies, ++i_s) {
            GraphResetTemps();

            auto & myspc = mypc_ref.GetParticleContainer(i_s);
            for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti) {
                (*g_temp[lev])[mfi.index()] += pti.numParticles();
            }

            weight_name = "n" + mypc_ref.GetSpeciesNames()[i_s];
            graph.addLayoutData(*g_temp[lev], fab_name, sizeof(amrex::Real),
                            weight_name, costs_heuristic_particles_wt);
        }
    }
}

void
WarpX::GraphAddLoadBalance(const int lev, const bool do_load_balance, const Vector<int>& new_dm,
                           const Real currentEfficiency, const Real proposedEfficiency)
{
    std::string fab_name = GraphFabName(lev);
    std::string wgt_name = "costs";

    // Add costs, with scale equaling currentEfficiency
    graph.addLayoutData(*costs[lev], fab_name, sizeof(amrex::Real),
                        wgt_name, double(currentEfficiency));

    // Add generated mapping, with scale equaling proposedEfficiency
    // Use the name to show whether load balance occurred.
    wgt_name = (do_load_balance) ? "new_dm" : "unused_dm";
    graph.addNodeWeight(fab_name, wgt_name, new_dm, double(propsedEfficiency), true);
}

void
WarpX::GraphPrintandClear(const std::string& graph_name)
{
    graph.print_table(graph_name);
    graph.clear();
}
