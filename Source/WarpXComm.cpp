
#include <WarpX.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FillPatchUtil_F.H>
#include <WarpX_f.H>

#include <algorithm>
#include <cstdlib>

using namespace amrex;

void
WarpX::UpdateAuxiliaryData ()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {// xxxxx pml stuff???
        const auto& period = Geom(lev).periodicity();
        Bfield_fp[lev][0]->FillBoundary(period);
        Bfield_fp[lev][1]->FillBoundary(period);
        Bfield_fp[lev][2]->FillBoundary(period);
        Efield_fp[lev][0]->FillBoundary(period);
        Efield_fp[lev][1]->FillBoundary(period);
        Efield_fp[lev][2]->FillBoundary(period);
    }

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const auto& crse_period = Geom(lev-1).periodicity();
        const int ng = Bfield_cp[lev][0]->nGrow();
        const DistributionMapping& dm = Bfield_cp[lev][0]->DistributionMap();

        // B field
        {
            MultiFab dBx(Bfield_cp[lev][0]->boxArray(), dm, 1, ng);
            MultiFab dBy(Bfield_cp[lev][1]->boxArray(), dm, 1, ng);
            MultiFab dBz(Bfield_cp[lev][2]->boxArray(), dm, 1, ng);
            dBx.copy(*Bfield_aux[lev-1][0], 0, 0, 1, 0, ng, crse_period);
            dBy.copy(*Bfield_aux[lev-1][1], 0, 0, 1, 0, ng, crse_period);
            dBz.copy(*Bfield_aux[lev-1][2], 0, 0, 1, 0, ng, crse_period);
            // xxxxx PML stuff?
            Bfield_cp[lev][0]->FillBoundary(crse_period);
            Bfield_cp[lev][1]->FillBoundary(crse_period);
            Bfield_cp[lev][2]->FillBoundary(crse_period);
            MultiFab::Subtract(dBx, *Bfield_cp[lev][0], 0, 0, 1, ng);
            MultiFab::Subtract(dBy, *Bfield_cp[lev][1], 0, 0, 1, ng);
            MultiFab::Subtract(dBz, *Bfield_cp[lev][2], 0, 0, 1, ng);
            
            const Real* dx = Geom(lev-1).CellSize();
            const int ref_ratio = refRatio(lev-1)[0];
            {
                std::array<FArrayBox,BL_SPACEDIM> bfab;
                for (MFIter mfi(*Bfield_aux[lev][0]); mfi.isValid(); ++mfi)
                {
                    Box ccbx = mfi.fabbox();
                    ccbx.enclosedCells();
                    ccbx.coarsen(ref_ratio).refine(ref_ratio); // so that ccbx is coarsenable
                    
                    const FArrayBox& cxfab = dBx[mfi];
                    bfab[0].resize(amrex::convert(ccbx,Bx_nodal_flag));
#if (BL_SPACEDIM > 1)
                    const FArrayBox& cyfab = dBy[mfi];
                    bfab[1].resize(amrex::convert(ccbx,By_nodal_flag));
#endif
#if (BL_SPACEDIM > 2)
                    const FArrayBox& czfab = dBz[mfi];
                    bfab[2].resize(amrex::convert(ccbx,Bz_nodal_flag));
#endif
                    amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                                 D_DECL(BL_TO_FORTRAN_ANYD(bfab[0]),
                                                        BL_TO_FORTRAN_ANYD(bfab[1]),
                                                        BL_TO_FORTRAN_ANYD(bfab[2])),
                                                 D_DECL(BL_TO_FORTRAN_ANYD(cxfab),
                                                        BL_TO_FORTRAN_ANYD(cyfab),
                                                        BL_TO_FORTRAN_ANYD(czfab)),
                                                 dx, &ref_ratio);
                    
                    for (int idim = 0; idim < BL_SPACEDIM; ++idim)
                    {
                        FArrayBox& aux = (*Bfield_aux[lev][idim])[mfi];
                        FArrayBox& fp  =  (*Bfield_fp[lev][idim])[mfi];
                        const Box& bx = aux.box();
                        aux.copy(fp, bx, 0, bx, 0, 1);
                        aux.plus(bfab[idim], bx, bx, 0, 0, 1);
                    }
                }
            }
        }

        // E field
        {
            MultiFab dEx(Efield_cp[lev][0]->boxArray(), dm, 1, ng);
            MultiFab dEy(Efield_cp[lev][1]->boxArray(), dm, 1, ng);
            MultiFab dEz(Efield_cp[lev][2]->boxArray(), dm, 1, ng);
            dEx.copy(*Efield_aux[lev-1][0], 0, 0, 1, 0, ng, crse_period);
            dEy.copy(*Efield_aux[lev-1][1], 0, 0, 1, 0, ng, crse_period);
            dEz.copy(*Efield_aux[lev-1][2], 0, 0, 1, 0, ng, crse_period);
            // xxxxx PML stuff?
            Efield_cp[lev][0]->FillBoundary(crse_period);
            Efield_cp[lev][1]->FillBoundary(crse_period);
            Efield_cp[lev][2]->FillBoundary(crse_period);
            MultiFab::Subtract(dEx, *Efield_cp[lev][0], 0, 0, 1, ng);
            MultiFab::Subtract(dEy, *Efield_cp[lev][1], 0, 0, 1, ng);
            MultiFab::Subtract(dEz, *Efield_cp[lev][2], 0, 0, 1, ng);
            
            const Real* dx = Geom(lev-1).CellSize();
            const int ref_ratio = refRatio(lev-1)[0];
            {
                std::array<FArrayBox,BL_SPACEDIM> efab;
                for (MFIter mfi(*Efield_aux[lev][0]); mfi.isValid(); ++mfi)
                {
                    Box ccbx = mfi.fabbox();
                    ccbx.enclosedCells();
                    ccbx.coarsen(ref_ratio).refine(ref_ratio); // so that ccbx is coarsenable
                    
                    const FArrayBox& cxfab = dEx[mfi];
                    efab[0].resize(amrex::convert(ccbx,Ex_nodal_flag));
#if (BL_SPACEDIM > 1)
                    const FArrayBox& cyfab = dEy[mfi];
                    efab[1].resize(amrex::convert(ccbx,Ey_nodal_flag));
#endif
#if (BL_SPACEDIM > 2)
                    const FArrayBox& czfab = dEz[mfi];
                    efab[2].resize(amrex::convert(ccbx,Ez_nodal_flag));
#endif
                    amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                        D_DECL(BL_TO_FORTRAN_ANYD(efab[0]),
                                               BL_TO_FORTRAN_ANYD(efab[1]),
                                               BL_TO_FORTRAN_ANYD(efab[2])),
                                        D_DECL(BL_TO_FORTRAN_ANYD(cxfab),
                                               BL_TO_FORTRAN_ANYD(cyfab),
                                               BL_TO_FORTRAN_ANYD(czfab)),
                                        dx, &ref_ratio);
                    
                    for (int idim = 0; idim < BL_SPACEDIM; ++idim)
                    {
                        FArrayBox& aux = (*Efield_aux[lev][idim])[mfi];
                        FArrayBox& fp  =  (*Efield_fp[lev][idim])[mfi];
                        const Box& bx = aux.box();
                        aux.copy(fp, bx, 0, bx, 0, 1);
                        aux.plus(efab[idim], bx, bx, 0, 0, 1);
                    }
                }
            }
        }
    }
}

void
WarpX::FillBoundaryB ()
{
// xxxxx
#if 0
    FillBoundaryB(0, true); // level 0

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const IntVect& ref_ratio = refRatio(lev-1);

        const std::array<MultiFab,BL_SPACEDIM> crse {
                    MultiFab(*Bfield[lev-1][0], amrex::make_alias, 0, 1),
#if (BL_SPACEDIM == 3)
                    MultiFab(*Bfield[lev-1][1], amrex::make_alias, 0, 1),
#endif
                    MultiFab(*Bfield[lev-1][2], amrex::make_alias, 0, 1) };

        std::array<MultiFab,BL_SPACEDIM> fine {
                    MultiFab(*Bfield[lev][0], amrex::make_alias, 0, 1),
#if (BL_SPACEDIM == 3)
                    MultiFab(*Bfield[lev][1], amrex::make_alias, 0, 1),
#endif
                    MultiFab(*Bfield[lev][2], amrex::make_alias, 0, 1) };

        amrex::InterpCrseFineBndryEMfield(amrex::InterpB, crse, fine, Geom(lev-1), Geom(lev), ref_ratio[0]);

        FillBoundaryB(lev, true);
    }
#endif
}

void
WarpX::FillBoundaryE ()
{
// xxxxx
#if 0
    FillBoundaryE(0, true); // level 0

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const IntVect& ref_ratio = refRatio(lev-1);

        const std::array<MultiFab,BL_SPACEDIM> crse {
                    MultiFab(*Efield[lev-1][0], amrex::make_alias, 0, 1),
#if (BL_SPACEDIM == 3)
                    MultiFab(*Efield[lev-1][1], amrex::make_alias, 0, 1),
#endif
                    MultiFab(*Efield[lev-1][2], amrex::make_alias, 0, 1) };

        std::array<MultiFab,BL_SPACEDIM> fine {
                    MultiFab(*Efield[lev][0], amrex::make_alias, 0, 1),
#if (BL_SPACEDIM == 3)
                    MultiFab(*Efield[lev][1], amrex::make_alias, 0, 1),
#endif
                    MultiFab(*Efield[lev][2], amrex::make_alias, 0, 1) };

        amrex::InterpCrseFineBndryEMfield(amrex::InterpE, crse, fine, Geom(lev-1), Geom(lev), ref_ratio[0]);

        FillBoundaryE(lev, true);
    }
#endif
}

void
WarpX::FillBoundaryE(int lev, bool force)
{
// xxxxx
#if 0
    if (force || WarpX::nox > 1 || WarpX::noy > 1 || WarpX::noz > 1)
    {
// xxxxx
#if 0
        if (do_pml && lev == 0) {
            WarpX::ExchangeWithPML(*Efield[lev][0], *pml_E[0], geom[lev]);
            WarpX::ExchangeWithPML(*Efield[lev][1], *pml_E[1], geom[lev]);
            WarpX::ExchangeWithPML(*Efield[lev][2], *pml_E[2], geom[lev]);

            (*pml_E[0]).FillBoundary( geom[lev].periodicity() );
            (*pml_E[1]).FillBoundary( geom[lev].periodicity() );
            (*pml_E[2]).FillBoundary( geom[lev].periodicity() );
        }
#endif

        (*Efield[lev][0]).FillBoundary( geom[lev].periodicity() );
        (*Efield[lev][1]).FillBoundary( geom[lev].periodicity() );
        (*Efield[lev][2]).FillBoundary( geom[lev].periodicity() );
    }
#endif
}

void
WarpX::FillBoundaryB(int lev, bool force)
{
// xxxxx
#if 0
    if (force || WarpX::nox > 1 || WarpX::noy > 1 || WarpX::noz > 1)
    {
// xxxxx
#if 0
        if (do_pml && lev == 0) {
            WarpX::ExchangeWithPML(*Bfield[lev][0], *pml_B[0], geom[lev]);
            WarpX::ExchangeWithPML(*Bfield[lev][1], *pml_B[1], geom[lev]);
            WarpX::ExchangeWithPML(*Bfield[lev][2], *pml_B[2], geom[lev]);

            (*pml_B[0]).FillBoundary( geom[lev].periodicity() );
            (*pml_B[1]).FillBoundary( geom[lev].periodicity() );
            (*pml_B[2]).FillBoundary( geom[lev].periodicity() );
        }
#endif

        (*Bfield[lev][0]).FillBoundary( geom[lev].periodicity() );
        (*Bfield[lev][1]).FillBoundary( geom[lev].periodicity() );
        (*Bfield[lev][2]).FillBoundary( geom[lev].periodicity() );
    }
#endif
}

void
WarpX::AverageDownB ()
{
// xxxxx
#if 0
    for (int lev = finest_level; lev > 0; --lev)
    {
        const IntVect& ref_ratio = refRatio(lev-1);
        const Geometry& crse_geom = Geom(lev-1);
#if (BL_SPACEDIM == 3)
        Array<const MultiFab*> fine {Bfield[lev][0].get(), Bfield[lev][1].get(), Bfield[lev][2].get()};
#else
        Array<const MultiFab*> fine {Bfield[lev][0].get(), Bfield[lev][2].get()};
#endif
        Array<std::unique_ptr<MultiFab> > crse(BL_SPACEDIM);
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            BoxArray cba = fine[idim]->boxArray();
            const DistributionMapping& dm = fine[idim]->DistributionMap();
            cba.coarsen(ref_ratio);
            crse[idim].reset(new MultiFab(cba, dm, 1, 0));
        }

        amrex::average_down_faces(fine, amrex::GetArrOfPtrs(crse), ref_ratio);

#if (BL_SPACEDIM == 3)
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            Bfield[lev-1][idim]->copy(*crse[idim], crse_geom.periodicity());
        }
#else
        Bfield[lev-1][0]->copy(*crse[0], crse_geom.periodicity());
        Bfield[lev-1][2]->copy(*crse[1], crse_geom.periodicity());
        amrex::average_down(*Bfield[lev][1], *Bfield[lev-1][1], 0, 1, ref_ratio);
#endif        
    }
#endif
}

void
WarpX::AverageDownE ()
{
//xxxxx    AverageDownEorCurrent(Efield);
}

void
WarpX::AverageDownCurrent ()
{
//xxxxx    AverageDownEorCurrent(current);
}

void
WarpX::AverageDownEorCurrent (Array<std::array<std::unique_ptr<amrex::MultiFab>,3> >& data)
{
    for (int lev = finest_level; lev > 0; --lev)
    {
        const IntVect& ref_ratio = refRatio(lev-1);
        const auto& crse_period = Geom(lev-1).periodicity();
#if (BL_SPACEDIM == 3)
        Array<const MultiFab*> fine {data[lev][0].get(), data[lev][1].get(), data[lev][2].get()};
#else
        Array<const MultiFab*> fine {data[lev][0].get(), data[lev][2].get()};
#endif
        
        Array<std::unique_ptr<MultiFab> > crse(BL_SPACEDIM);
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            BoxArray cba = fine[idim]->boxArray();
            const DistributionMapping& dm = fine[idim]->DistributionMap();
            cba.coarsen(ref_ratio);
            crse[idim].reset(new MultiFab(cba, dm, 1, 0));
        }

        amrex::average_down_edges(fine, amrex::GetArrOfPtrs(crse), ref_ratio);

#if (BL_SPACEDIM == 3)
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            data[lev-1][idim]->copy(*crse[idim], crse_period);
        }
#else
        data[lev-1][0]->copy(*crse[0], crse_period);
        data[lev-1][2]->copy(*crse[1], crse_period);
//xxxxx amrex::average_down_nodes(*data[lev][1], *data[lev-1][1], 0, 1, ref_ratio);
#endif        
    }
}

void
WarpX::SyncCurrent ()
{

            // Restrict fine patch current onto the coarse patch

            // sumBoundary????

            // Add fine level's coarse patch current onto this level's fine patch


// xxxxx
#if 0
    for (int lev = finest_level; lev >= 0; --lev)
    {
        const auto& period = Geom(lev).periodicity();
        current[lev][0]->SumBoundary(period);
        current[lev][1]->SumBoundary(period);
        current[lev][2]->SumBoundary(period);

        // xxxxx
#if 0
        if (lev > 0) {
            bndry4fine[lev-1]->addTo(*current[lev][0], *current[lev][1], *current[lev][2]);
        }
        if (lev < finest_level) {
            bndry4crse[lev+1]->addTo(*current[lev][0], *current[lev][1], *current[lev][2]);
        }
#endif
    }

    AverageDownCurrent();
#endif
}


