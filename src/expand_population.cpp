#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector ATTEMPTED_NEW_INCORRECT_do_expand_population(NumericVector src,
                                   IntegerVector target_dims,
                                   IntegerVector src_to_target_dim_map)
{
    //Set up the mapping to source indexing
    int n_src_dims = src_to_target_dim_map.length();

    // n before in target
    int n_target_dims = target_dims.length();
    int num_before_in_target[target_dims.length()];
    num_before_in_target[0] = 1;
    for (int i=1; i<n_target_dims; i++)
        num_before_in_target[i] = num_before_in_target[i-1] * target_dims[i-1];

    // target is in source
    bool target_is_in_src[n_target_dims];
    for (int i=0; i<n_target_dims; i++)
    {
        target_is_in_src[i] = false;
        for (int j=0; j<n_src_dims; j++)
        {
            if (src_to_target_dim_map[j]==i)
            {
                target_is_in_src[i] = true;
                break;
            }
        }
    }

    // n before in src
    int num_in_src_dim[n_src_dims];
    int num_before_in_src[n_src_dims];
    num_before_in_src[0] = 1;

    for (int j=0; j<n_src_dims; j++)
    {
        num_in_src_dim[j] = target_dims[src_to_target_dim_map[j]];
        if (j>0)
            num_before_in_src[j] = num_before_in_src[j-1] * num_in_src_dim[j-1];
    }

    // invert the mapping to go from target to src
    int target_to_src_dim_map[n_target_dims];
    for (int j=0; j<n_target_dims; j++)
    {
        for (int i=0; i<n_src_dims; i++)
        {
            if (src_to_target_dim_map[i]==j)
            {
                target_to_src_dim_map[j] = i;
                break;
            }
        }
    }

    //Calculate the length of the rv
    int len = 1;
    for (int i=0; i<target_dims.length(); i++)
        len *= target_dims[i];

    // set up the rv
    NumericVector rv(len);

    // set up to crunch the indices
    int indices_into_src[len];
    int i_in_src_dim;
    int d_src;

    if (target_is_in_src[0])
    {
        d_src = target_to_src_dim_map[0];
        for (int i_in_target_dim=0; i_in_target_dim<target_dims[0]; i_in_target_dim++)
        {
            //for now;
            i_in_src_dim = i_in_target_dim;
            indices_into_src[i_in_target_dim] = i_in_src_dim * num_before_in_src[d_src];
        }
    }
    else
    {
        for (int i_in_target_dim=0; i_in_target_dim<target_dims[0]; i_in_target_dim++)
            indices_into_src[i_in_target_dim] = 0;
    }

    for (int d_target=1; d_target<n_target_dims; d_target++)
    {
        if (target_is_in_src[d_target])
        {
            d_src = target_to_src_dim_map[d_target];
            for (int i_in_target_dim=target_dims[d_target]-1; i_in_target_dim>0; i_in_target_dim--)
            {
                //for now;
                i_in_src_dim = i_in_target_dim;

                for (int j_before=0; j_before<num_before_in_target[d_target]; j_before++)
                {
                    indices_into_src[num_before_in_target[d_target] * i_in_target_dim + j_before] = i_in_src_dim * num_before_in_src[d_src] +
                        indices_into_src[j_before]; //the value for the previous indices
                }
            }
        }
        else
        {
            for (int i_in_d=target_dims[d_target]-1; d_target>0; d_target--)
            {
                for (int j_before=0; j_before<num_before_in_target[d_target]; j_before++)
                {
                    indices_into_src[num_before_in_target[d_target] * i_in_d + j_before] = indices_into_src[j_before]; //the value for the previous indices
                }
            }
        }
    }

    //Loop through and pull from src into rv
    for (int i=0; i<len; i++)
    {
        rv[i] = src[indices_into_src[i]];
    }

    //Return
    return (rv);
}

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
