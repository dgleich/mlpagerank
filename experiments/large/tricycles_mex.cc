#include <mex.h>

#include <vector>
#include <set>
#include <map>

// everything in this script has a horrible name.

// adjust "n+1" in the constructor to handle Matlab's "1" based indices.

struct tris_calculator {
    std::vector< std::vector< int > > out_edges;
    std::vector< std::vector< int > > in_edges;
    size_t n;
    size_t nedges;
    int* esrc;
    int* edst;

    std::vector< int > is_recip;
    
    std::vector< int > c1;
    std::vector< int > c2;
    std::vector< int > c3;

    triweights_calculator(size_t n_, int* src, int* dst, size_t nedges_) 
    : out_edges(n_+1), in_edges(n_+1), n(n_+1), nedges(nedges_), esrc(src), 
      edst(dst), is_recip(nedges)
    {
        std::vector< std::map< int, size_t > > recip(n);

        for (size_t i=0; i<nedges; ++i) {
            int s = src[i];
            int d = dst[i];

            // the 
            out_edges[s].push_back(i);
            in_edges[d].push_back(i);

            // check if the reciprocal edge (d, s) exists
            if (recip[d].count(s) > 0) {
                is_recip[i] = 1;
                // look up the reciprocal edge
                is_recip[recip[d][s]] = 1;
            } else {
                // insert this edge
                recip[s][d] = i;
            }
        }
    }

    double _assign(size_t eid1, size_t eid2, size_t eid3, double *array) {
        array[eid1] += 1;
        array[eid2] += 1;
        array[eid3] += 1;
    }

    void build() {
        // for each edge...
        for (size_t ei=0; ei<nedges; ++ei) {
            int es = esrc[ei];
            int ed = edst[ei];

            std::map<int, size_t> neighs;

            // for this edge, look at the list of out-neighbors from dst
            // and look at the set of in-neighbors for src.
            // any vertex in common completes the cycle.

            for (int i=0; i<out_edges[ed].size(); ++i) {
                int outdst_id = out_edges[ed][i];
                int dstdst = edst[outdst_id];
                if (dstdst != es) {
                    neighs[dstdst] = outdst_id;
                }
            }

            for (int i=0; i<in_edges[es].size() && neighs.size() > 0; ++i) {
                int insrc_id = in_edges[es][i]; // the id of the edge to src
                int srcsrc = esrc[insrc_id]; // the vertex id of the edge to src.

                if (srcsrc == ed) { continue; }
                
                if (neighs.count(srcsrc) > 0) {
                    //mexPrintf("%i %i %i    %i %i %i\n", es, ed, srcsrc, insrc_id, ei, neighs[srcsrc]);
                    // we have a triangle (insrc_id, ei, neighs[srcsrc])
                    int nrecip = 0;

                    c1.push_back(es);
                    c2.push_back(ed);
                    c3.push_back(srcsrc);
                }
            }
        }
    }
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t n = (size_t) mxGetScalar(prhs[0]);
    int *ai = (int*) mxGetPr(prhs[1]);
    int *aj = (int*) mxGetPr(prhs[2]);
    size_t nedges = (size_t) mxGetNumberOfElements(prhs[1]);

    triweights_calculator calc(n, ai, aj, nedges);

    calc.build();
    
    size_t ntris = calc.c1.size();

    plhs[0] = mxCreateDoubleMatrix(ntris, 3, mxREAL);    
    double *T = mxGetPr(plhs[0]);
    for (size_t i=0; i<ntris; ++i) {
        T[i] = calc.c1[i];
        T[i+ntris] = calc.c2[i];
        T[i+2*ntris] = calc.c3[i];
    }
}

