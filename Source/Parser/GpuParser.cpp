#include <GpuParser.H>

GpuParser::GpuParser (WarpXParser const& wp)
{
#ifdef _OPENMP
    struct wp_parser* a_wp = wp.m_parser[0];
#else
    struct wp_parser* a_wp = wp.m_parser;
#endif
    m_parser.sz_mempool = wp_ast_size((struct wp_node*)a_wp);
    m_parser.p_root = (struct wp_node*)
        amrex::The_Managed_Arena()->alloc(m_parser.sz_mempool);
    m_parser.p_free = m_parser.p_root;
    m_parser.ast = wp_parser_ast_dup(&m_parser, a_wp->ast, 0); // 0: don't free the source
#ifdef AMREX_USE_GPU
    wp_parser_regvar_gpu(&m_parser, "x", 0);
    wp_parser_regvar_gpu(&m_parser, "y", 1);
    wp_parser_regvar_gpu(&m_parser, "z", 2);
#else
    wp_parser_regvar(&m_parser, "x", var);
    wp_parser_regvar(&m_parser, "y", var+1);
    wp_parser_regvar(&m_parser, "z", var+2);
#endif
}

void
GpuParser::clear ()
{
    amrex::The_Managed_Arena()->free(m_parser.ast);
}


