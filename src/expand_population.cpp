#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector do_expand_population(NumericVector src,
                                   IntegerVector target_dims,
                                   IntegerVector src_to_target_dim_map)
{
    //Set up the mapping to source indexing
    int n_src_dims = src_to_target_dim_map.length();
    int src_total_index;
    int src_dim_index;

    int num_before_in_target[target_dims.length()];
    num_before_in_target[0] = 1;
    for (int i=1; i<target_dims.length(); i++)
        num_before_in_target[i] = num_before_in_target[i-1] * target_dims[i-1];

    int num_in_src_dim[n_src_dims];
    int num_before_in_src[n_src_dims];
    num_before_in_src[0] = 1;
    int num_before_in_target_by_src_index[n_src_dims];

    for (int j=0; j<n_src_dims; j++)
    {
        num_in_src_dim[j] = target_dims[src_to_target_dim_map[j]];
        num_before_in_target_by_src_index[j] = num_before_in_target[src_to_target_dim_map[j]];
        if (j>0)
            num_before_in_src[j] = num_before_in_src[j-1] * num_in_src_dim[j-1];
    }

    //Calculate the length of the rv
    int len = 1;
    for (int i=0; i<target_dims.length(); i++)
        len *= target_dims[i];
    NumericVector rv(len);

    //Loop through and pull from src into rv
    for (int i=0; i<len; i++)
    {
        src_total_index = 0;
        for (int j=0; j<n_src_dims; j++)
        {
            src_dim_index =  (i / num_before_in_target_by_src_index[j]) % num_in_src_dim[j];
            src_total_index += src_dim_index * num_before_in_src[j];
        }

        rv[i] = src[src_total_index];
    }

    //Return
    return (rv);
}
