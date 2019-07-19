#ifndef WARPX_GPU_PARSER_H_
#define WARPX_GPU_PARSER_H_

#include <WarpXParser.H>
#include <AMReX_Gpu.H>

class GpuParser
{
public:
    GpuParser (WarpXParser const& wp);
    void clear ();

    AMREX_GPU_HOST_DEVICE
    double
    operator() (double x, double y, double z) const noexcept
    {
#ifdef AMREX_USE_GPU
        amrex::Gpu::SharedMemory<double> gsm;
        double* p = gsm.dataPtr();
        int tid = threadIdx.x + threadIdx.y*blockDim.x + threadIdx.z*(blockDim.x*blockDim.y);
        p[tid*3] = x;
        p[tid*3+1] = y;
        p[tid*3+2] = z;
#else
        var[0] = x; var[1] = y; var[2] = z;
#endif
        return wp_ast_eval(m_parser.ast);
    }
private:
    struct wp_parser m_parser;
#ifndef AMREX_USE_GPU
    mutable double var[3];
#endif
};

#endif
