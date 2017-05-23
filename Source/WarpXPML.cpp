
#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpXPML.H>

using namespace amrex;

#if 0
void
WarpX::InitPML ()
{
    if (!do_pml) return;

    const Geometry& gm0 = Geom(0);
    const Box& domainbox = gm0.Domain();
    Box grownbox = domainbox;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        if (!Geometry::isPeriodic(i)) {
            grownbox.grow(i,pml_ncell);
        }
    }

    int block_size = maxGridSize(0);
    while (block_size < pml_ncell) {
        block_size += blockingFactor(0);
    }

    Array<IntVect> shift;
    {
        int len[3] = {0,0,0};
        int jmp[3] = {1,1,1};

        for (int i = 0; i < BL_SPACEDIM; ++i) {
            if (!Geometry::isPeriodic(i)) {
                len[i] = jmp[i] = domainbox.length(i);
            }
        }

        for (int i = -len[0]; i <= len[0]; i += jmp[0]) {
        for (int j = -len[1]; j <= len[1]; j += jmp[1]) {
        for (int k = -len[2]; k <= len[2]; k += jmp[2]) {
            if (i != 0 || j != 0 || k!= 0) {
                shift.push_back(IntVect(D_DECL(i,j,k)));
            }
        }
        }
        }
    }

    BoxList bl;
    for (const IntVect& iv : shift) {
        Box bbx = domainbox;
        bbx.shift(iv);
        bbx &= grownbox;
        BoxList bltmp(bbx);
        bltmp.maxSize(block_size);
        bl.catenate(bltmp);
    }
    pml_ba.define(bl);
    
    // xxxxx for now let's just use a simple distributionmapping
    pml_dm.RoundRobinProcessorMap(pml_ba.size(), ParallelDescriptor::NProcs());

    const int ng = Efield[0][0]->nGrow();
    pml_E[0].reset(new MultiFab(amrex::convert(pml_ba,Ex_nodal_flag),pml_dm,2,ng));
    pml_E[1].reset(new MultiFab(amrex::convert(pml_ba,Ey_nodal_flag),pml_dm,2,ng));
    pml_E[2].reset(new MultiFab(amrex::convert(pml_ba,Ez_nodal_flag),pml_dm,2,ng));
    pml_B[0].reset(new MultiFab(amrex::convert(pml_ba,Bx_nodal_flag),pml_dm,2,ng));
    pml_B[1].reset(new MultiFab(amrex::convert(pml_ba,By_nodal_flag),pml_dm,2,ng));
    pml_B[2].reset(new MultiFab(amrex::convert(pml_ba,Bz_nodal_flag),pml_dm,2,ng));

    pml_E[0]->setVal(0.0);
    pml_E[1]->setVal(0.0);
    pml_E[2]->setVal(0.0);
    pml_B[0]->setVal(0.0);
    pml_B[1]->setVal(0.0);
    pml_B[2]->setVal(0.0);

    {
        const Real* dx = gm0.CellSize();
        const int* dlo = domainbox.loVect();
        const int* dhi = domainbox.hiVect();
        const int* glo = grownbox.loVect();
        const int* ghi = grownbox.hiVect();
        const IntVect& sz = grownbox.size();
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            pml_sigma     [idim].resize(sz[idim]+1, 0.0);
            pml_sigma_star[idim].resize(sz[idim]  , 0.0);
            pml_sigma     [idim].m_lo = glo[idim];
            pml_sigma     [idim].m_hi = ghi[idim]+1;
            pml_sigma_star[idim].m_lo = glo[idim];
            pml_sigma_star[idim].m_hi = ghi[idim];

            const Real fac = 4.0*PhysConst::c/(dx[idim]*static_cast<Real>(pml_ncell*pml_ncell));

            for (int ind = glo[idim]; ind < dlo[idim]; ++ind) {
                Real offset = static_cast<Real>(dlo[idim] - ind);
                pml_sigma[idim][ind-glo[idim]] = fac*(offset*offset);
            }
            for (int ind = dhi[idim]+2; ind < ghi[idim]+2; ++ind) {
                Real offset = static_cast<Real>(ind - (dhi[idim]+1));
                pml_sigma[idim][ind-glo[idim]] = fac*(offset*offset);
            }

            for (int icc = glo[idim]; icc < dlo[idim]; ++icc) {
                Real offset = static_cast<Real>(dlo[idim] - icc) - 0.5;
                pml_sigma_star[idim][icc-glo[idim]] = fac*(offset*offset);
            }
            for (int icc = dhi[idim]+1; icc < ghi[idim]+1; ++icc) {
                Real offset = static_cast<Real>(icc - dhi[idim]) - 0.5;
                pml_sigma_star[idim][icc-glo[idim]] = fac*(offset*offset);                
            }
        }
    }
}


void
WarpX::ComputePMLFactors (int lev, Real dt)
{
    ComputePMLFactorsE(lev, dt);
    ComputePMLFactorsB(lev, dt);
}

void
WarpX::ComputePMLFactorsE (int lev, Real dt)
{
    const Real* dx = Geom(lev).CellSize();

    const std::array<Real,BL_SPACEDIM> dtsdx {D_DECL(dt/dx[0], dt/dx[1], dt/dx[2])};

    const Real c2 = PhysConst::c*PhysConst::c;
    const std::array<Real,BL_SPACEDIM> dtsdx_c2 {D_DECL(dtsdx[0]*c2, dtsdx[1]*c2, dtsdx[2]*c2)};

    for (int idim = 0; idim < BL_SPACEDIM; ++idim)
    {
        pml_sigma_fac1[idim].assign(pml_sigma[idim].size(), 1.0);
        pml_sigma_fac2[idim].assign(pml_sigma[idim].size(), dtsdx_c2[idim]);
        pml_sigma_fac1[idim].m_lo = pml_sigma[idim].m_lo;
        pml_sigma_fac1[idim].m_hi = pml_sigma[idim].m_hi;
        pml_sigma_fac2[idim].m_lo = pml_sigma[idim].m_lo;
        pml_sigma_fac2[idim].m_hi = pml_sigma[idim].m_hi;

        if (!Geometry::isPeriodic(idim))
        {
            for (int i = 0; i < pml_ncell; ++i)
            {
                pml_sigma_fac1[idim][i] = std::exp(-pml_sigma[idim][i]*dt);
                pml_sigma_fac2[idim][i] = (1.0-pml_sigma_fac1[idim][i])
                    / (pml_sigma[idim][i]*dt) * dtsdx_c2[idim];
            }

            for (int iiii = 0; iiii < pml_ncell; ++iiii)
            {
                int i = (iiii - pml_ncell) + static_cast<int>(pml_sigma_fac1[idim].size());
                pml_sigma_fac1[idim][i] = std::exp(-pml_sigma[idim][i]*dt);
                pml_sigma_fac2[idim][i] = (1.0-pml_sigma_fac1[idim][i])
                    / (pml_sigma[idim][i]*dt) * dtsdx_c2[idim];
            }
        }
    }
}

void
WarpX::ComputePMLFactorsB (int lev, Real dt)
{
    const Real* dx = Geom(lev).CellSize();

    const std::array<Real,BL_SPACEDIM> dtsdx {D_DECL(dt/dx[0], dt/dx[1], dt/dx[2])};

    for (int idim = 0; idim < BL_SPACEDIM; ++idim)
    {
        pml_sigma_star_fac1[idim].assign(pml_sigma_star[idim].size(), 1.0);
        pml_sigma_star_fac2[idim].assign(pml_sigma_star[idim].size(), dtsdx[idim]);
        pml_sigma_star_fac1[idim].m_lo = pml_sigma_star[idim].m_lo;
        pml_sigma_star_fac1[idim].m_hi = pml_sigma_star[idim].m_hi;
        pml_sigma_star_fac2[idim].m_lo = pml_sigma_star[idim].m_lo;
        pml_sigma_star_fac2[idim].m_hi = pml_sigma_star[idim].m_hi;
        
        if (!Geometry::isPeriodic(idim))
        {
            for (int i = 0; i < pml_ncell; ++i)
            {
                pml_sigma_star_fac1[idim][i] = std::exp(-pml_sigma_star[idim][i]*dt);
                pml_sigma_star_fac2[idim][i] = (1.0-pml_sigma_star_fac1[idim][i])
                    / (pml_sigma_star[idim][i]*dt) * dtsdx[idim];
            }

            for (int iiii = 0; iiii < pml_ncell; ++iiii)
            {
                int i = (iiii - pml_ncell) + static_cast<int>(pml_sigma_star_fac1[idim].size());

                pml_sigma_star_fac1[idim][i] = std::exp(-pml_sigma_star[idim][i]*dt);
                pml_sigma_star_fac2[idim][i] = (1.0-pml_sigma_star_fac1[idim][i])
                    / (pml_sigma_star[idim][i]*dt) * dtsdx[idim];
            }
        }
    }
}

void
WarpX::ExchangeWithPML (MultiFab& regmf, MultiFab& pmlmf, const Geometry& gm)
{
    if (!do_pml) return;

    // Copy from PML to regular data

    MultiFab totpmlmf(pmlmf.boxArray(), pmlmf.DistributionMap(), 1, 0);
    MultiFab::LinComb(totpmlmf, 1.0, pmlmf, 0, 1.0, pmlmf, 1, 0, 1, 0);

    MultiFab tmpregmf(regmf.boxArray(), regmf.DistributionMap(), 1, regmf.nGrow());
    tmpregmf.copy(totpmlmf, 0, 0, 1, 0, regmf.nGrow(), gm.periodicity());

    Box dom = gm.Domain();
    dom.convert(regmf.ixType());
    for (int idim=0; idim < BL_SPACEDIM; ++idim) {
        if (Geometry::isPeriodic(idim)) {
            dom.grow(idim, regmf.nGrow());
        }
    }
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(regmf); mfi.isValid(); ++mfi)
    {
        const BoxList& bl = amrex::boxDiff(mfi.fabbox(), dom);
        for (const Box& bx : bl)
        {
            regmf[mfi].copy(tmpregmf[mfi], bx, 0, bx, 0, 1);
        }
    }

    // Copy from regular data to PML's first component
    pmlmf.copy (regmf, 0, 0, 1, 0, pmlmf.nGrow(), gm.periodicity());
    // Zero out the second component
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(pmlmf); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.fabbox();
        bx &= dom;
        if (bx.ok()) {
            pmlmf[mfi].setVal(0.0, bx, 1, 1);
        }
    }
}

#endif
