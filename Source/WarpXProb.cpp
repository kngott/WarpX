#include <WarpX.H>

using namespace amrex;

void
WarpX::InitLevelData (int lev, Real time)
{
    for (int i = 0; i < 3; ++i) {
	current_fp[lev][i]->setVal(0.0);
	Efield_fp[lev][i]->setVal(0.0);
	Bfield_fp[lev][i]->setVal(0.0);
    }

    if (lev > 0) {
        for (int i = 0; i < 3; ++i) {
            Efield_aux[lev][i]->setVal(0.0);
            Bfield_aux[lev][i]->setVal(0.0);

            current_cp[lev][i]->setVal(0.0);
            Efield_cp[lev][i]->setVal(0.0);
            Bfield_cp[lev][i]->setVal(0.0);
        }
    }
}
