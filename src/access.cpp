#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector do_get_access_indices(IntegerVector dims,
                                    List to_access)
{
    int n_dim = dims.size();
    int n_before[n_dim];
    int n_before_access[n_dim];
    int access_dims[n_dim];

    n_before[0] = 1;
    for (int i=1; i<n_dim; i++)
    {
        n_before[i] = n_before[i-1] * dims[i-1];
    }

    n_before_access[0] = 1;
    for (int i=0; i<n_dim; i++)
    {
        access_dims[i] = ((IntegerVector) to_access[i]).size();
        if (i>0)
            n_before_access[i] = n_before_access[i-1] * access_dims[i-1];
    }
    int n = n_before_access[n_dim-1] * access_dims[n_dim-1];

    int dim_index;
    IntegerVector rv(n);

    for (int i=0; i<n; i++)
    {
        rv[i] = 0;
        for (int d=0; d<n_dim; d++)
        {
            dim_index = (i / n_before_access[d]) % access_dims[d];
            rv[i] += (((IntegerVector) to_access[d])[dim_index]-1) * n_before[d];
        }
    }

    return rv;
}
