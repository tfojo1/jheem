#include <Rcpp.h>
using namespace Rcpp;

//---------------//
//-- CONSTANTS --//
//---------------//

int N_HIV_DIMS = 8;
int N_NONHIV_DIMS = 6;
int N_GENERAL_DIMS = 5;

int AGE=0;
int RACE=1;
int SUBPOPULATION=2;
int SEX=3;
int RISK=4;
int NONHIV_SUBSET=5;
int CONTINUUM=5;
int CD4=6;
int HIV_SUBSET=7;

int FIXED_SIZE_BY_STRATA = 1;

//-------------------//
//-- ARRAY HELPERS --//
//-------------------//

//-- Multiply two arrays --//
void multiply_arr(double *m1, double *m2,
                            double *rv, int size)
{
    for (int i=0; i<size; i++)
        rv[i] = m1[i] * m2[i];
}

void multiply_arr(const NumericVector &m1, const NumericVector &m2,
                  double *rv)
{
    int size = m1.size();
    for (int i=0; i<size; i++)
        rv[i] = m1[i] * m2[i];
}

void multiply_arr(double *m1, const NumericVector &m2,
                  double *rv)
{
    int size = m2.size();
    for (int i=0; i<size; i++)
        rv[i] = m1[i] * m2[i];
}

//-- Multiply an array by another array and by a scalar --//
void double_multiply_arr(double *m1, double *m2, double m3,
                         double *rv, int size)
{
    for (int i=0; i<size; i++)
        rv[i] = m1[i] * m2[i] * m3;
}

//-- Divide arrays - if the denom is zero, put a zero in --//
void divide_arr_zero_if_zero_denom(double *numerator, double *denominator,
                double *rv, int size)
{
    for (int i=0; i<size; i++)
    {
        if (denominator[i]==0)
            rv[i] = 0;
        else
            rv[i] = numerator[i] / denominator[i];
    }
}




void copy_arr(double *to_copy, double *rv, int size)
{
    for (int i=0; i<size; i++)
        rv[i] = to_copy[i];
}

void clear_arr(double *arr, int size)
{
    for (int i=0; i<size; i++)
        arr[i] = 0;
}


//-- Add and subtract arrays--//
void increment_arr(double *rv, double *to_add, int size)
{
    for (int i=0; i<size; i++)
        rv[i] += to_add[i];
}

void decrement_arr(double *rv, double *to_add, int size)
{
    for (int i=0; i<size; i++)
        rv[i] -= to_add[i];
}

void add_arr(double *v1, double *v2, double *rv, int size)
{
    for (int i=0;i<size; i++)
        rv[i] = v1[i] + v2[i];
}

void add_arr(const NumericVector &v1, double *v2, double *rv)
{
    for (int i=0;i<v1.size(); i++)
        rv[i] = v1[i] + v2[i];
}

//-----------------------//
//-- LOW-LEVEL HELPERS --//
//-----------------------//

//-- Take the product of an array --//
int prod(const IntegerVector &v)
{
    int rv = 1;
    for (int i=0; i<v.size(); i++)
        rv = rv * v[i];
    return rv;
}

int prod_if(const IntegerVector &v, const LogicalVector &mask)
{
    int rv = 1;
    for (int i=0; i<v.size(); i++)
    {
        if (mask[i])
            rv = rv * v[i];
    }
    return rv;
}

//-- Cumulative product --//
IntegerVector cumprod_rightshifted(const IntegerVector &v)
{
    IntegerVector rv(v.size());
    rv[0] = 1;
    for (int i=1; i<v.size(); i++)
        rv[i] = rv[i-1] * v[i-1];

    return rv;
}

IntegerVector cumprod_righshifted_if(const IntegerVector &v, const IntegerVector &mask)
{
    IntegerVector rv(v.size());
    int current = 1;

    for (int i=0; i<v.size(); i++)
    {
        if (mask[i])
        {
            rv[i] = current;
            current = current * v[i];
        }
        else
            rv[i] = 0;
    }

    return rv;
}
IntegerVector cumprod_righshifted_if(const IntegerVector &v, const LogicalVector &mask)
{
    IntegerVector rv(v.size());
    int current = 1;

    for (int i=0; i<v.size(); i++)
    {
        if (mask[i])
        {
            rv[i] = current;
            current = current * v[i];
        }
        else
            rv[i] = 0;
    }

    return rv;
}

int max(const IntegerVector &v)
{
    int rv = v[0];
    for (int i=1; i<v.size(); i++)
    {
        if (v[i] > rv)
            rv = v[i];
    }

    return rv;
}

//-----------------------//
//-- MID-LEVEL HELPERS --//
//-----------------------//

void collapse_or_expand_8d(double *arr, double *rv,
                           const IntegerVector &full_dims, const LogicalVector &mask,
                           bool expand=false, bool overwrite=true)
{
    IntegerVector from_mult,to_mult;
    int N,from_index,to_index;
    if (expand)
    {
        from_mult = cumprod_righshifted_if(full_dims, mask);
        to_mult = cumprod_rightshifted(full_dims);
        N = prod(full_dims);
    }
    else
    {
        from_mult = cumprod_rightshifted(full_dims);
        to_mult = cumprod_righshifted_if(full_dims, mask);
        N = prod_if(full_dims, mask);
    }

    if (overwrite)
        clear_arr(rv, N);

    for (int i7=0; i7<full_dims[7]; i7++)
    {
        for (int i6=0; i6<full_dims[6]; i6++)
        {
            for (int i5=0; i5<full_dims[5]; i5++)
            {
                for (int i4=0; i4<full_dims[4]; i4++)
                {
                    for (int i3=0; i3<full_dims[3]; i3++)
                    {
                        for (int i2=0; i2<full_dims[2]; i2++)
                        {
                            for (int i1=0; i1<full_dims[1]; i1++)
                            {
                                for (int i0=0; i0<full_dims[0]; i0++)
                                {
                                    from_index = i0*from_mult[0] + i1*from_mult[1] +
                                        i2*from_mult[2] + i3*from_mult[3] + i4*from_mult[4] +
                                        i5*from_mult[5] + i6*from_mult[6] + i7*from_mult[7];

                                    to_index = i0*to_mult[0] + i1*to_mult[1] +
                                        i2*to_mult[2] + i3*to_mult[3] + i4*to_mult[4] +
                                        i5*to_mult[5] + i6*to_mult[6] + i7*to_mult[7];

                                    rv[to_index] += arr[from_index];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void collapse_or_expand_6d(double *arr, double *rv,
                           const IntegerVector &full_dims, const LogicalVector &mask,
                           bool expand=false, bool overwrite=true)
{
    IntegerVector from_mult,to_mult;
    int N,from_index,to_index;
    if (expand)
    {
        from_mult = cumprod_righshifted_if(full_dims, mask);
        to_mult = cumprod_rightshifted(full_dims);
        N = prod(full_dims);
    }
    else
    {
        from_mult = cumprod_rightshifted(full_dims);
        to_mult = cumprod_righshifted_if(full_dims, mask);
        N = prod_if(full_dims, mask);
    }

    if (overwrite)
        clear_arr(rv, N);

    for (int i5=0; i5<full_dims[5]; i5++)
    {
        for (int i4=0; i4<full_dims[4]; i4++)
        {
            for (int i3=0; i3<full_dims[3]; i3++)
            {
                for (int i2=0; i2<full_dims[2]; i2++)
                {
                    for (int i1=0; i1<full_dims[1]; i1++)
                    {
                        for (int i0=0; i0<full_dims[0]; i0++)
                        {
                            from_index = i0*from_mult[0] + i1*from_mult[1] +
                                i2*from_mult[2] + i3*from_mult[3] + i4*from_mult[4] +
                                i5*from_mult[5];

                            to_index = i0*to_mult[0] + i1*to_mult[1] +
                                i2*to_mult[2] + i3*to_mult[3] + i4*to_mult[4] +
                                i5*to_mult[5];

                            rv[to_index] += arr[from_index];
                        }
                    }
                }
            }
        }
    }
}

void collapse_or_expand_5d(double *arr, double *rv,
                           const IntegerVector &full_dims, const LogicalVector &mask,
                           bool expand=false, bool overwrite=true)
{
    IntegerVector from_mult,to_mult;
    int N,from_index,to_index;
    if (expand)
    {
        from_mult = cumprod_righshifted_if(full_dims, mask);
        to_mult = cumprod_rightshifted(full_dims);
        N = prod(full_dims);
    }
    else
    {
        from_mult = cumprod_rightshifted(full_dims);
        to_mult = cumprod_righshifted_if(full_dims, mask);
        N = prod_if(full_dims, mask);
    }

    if (overwrite)
        clear_arr(rv, N);

    for (int i4=0; i4<full_dims[4]; i4++)
    {
        for (int i3=0; i3<full_dims[3]; i3++)
        {
            for (int i2=0; i2<full_dims[2]; i2++)
            {
                for (int i1=0; i1<full_dims[1]; i1++)
                {
                    for (int i0=0; i0<full_dims[0]; i0++)
                    {
                        from_index = i0*from_mult[0] + i1*from_mult[1] +
                            i2*from_mult[2] + i3*from_mult[3] + i4*from_mult[4];

                        to_index = i0*to_mult[0] + i1*to_mult[1] +
                            i2*to_mult[2] + i3*to_mult[3] + i4*to_mult[4];

                        rv[to_index] += arr[from_index];
                    }
                }
            }
        }
    }
}

void rowsum(double *arr, double *rv, int nrow, int ncol, bool overwrite=true)
{
    if (overwrite)
        clear_arr(rv, nrow);

    for (int col=0; col<ncol; col++)
    {
        for (int row=0; row<nrow; row++)
        {
            rv[row] += arr[row + nrow*col];
        }
    }
}


// v is a vector of N_DROP * N_INVARIANT * N_FROM (col major)
// rv is a vector of length N_INVARIANT * N_TO (col major)
// mat represents a (col major) arrangement of N_DROP * N_INVARIANT matrices of size N_FROM x N_TO
//
// multiply such that rv[invariant, to] = v[[drop]][[invariant]][from] %*% mat[[drop]][[invariant]][from,to]
void do_piecewise_matrix_multiply_with_drop(double *v, double *mat, double *rv,
                                  int N_DROP, int N_INVARIANT, int N_FROM, int N_TO, bool overwrite=true)
{
    int V_DROP_MULT = 1;
    int V_INV_MULT = N_DROP;
    int V_FROM_MULT = V_INV_MULT * N_INVARIANT;

    int RV_INV_MULT = 1;
    int RV_TO_MULT = N_INVARIANT;

    int MAT_DROP_MULT = 1;
    int MAT_INV_MULT = N_DROP;
    int MAT_FROM_MULT = MAT_INV_MULT * N_INVARIANT;
    int MAT_TO_MULT = MAT_FROM_MULT * N_FROM;

    if (overwrite)
        clear_arr(rv, N_INVARIANT * N_TO);

    for (int to=0; to<N_TO; to++)
    {
        for (int from=0; from<N_FROM; from++)
        {
            for (int inv=0; inv<N_INVARIANT; inv++)
            {
                for (int drop=0; drop<N_DROP; drop++)
                {
                    rv[inv*RV_INV_MULT + to*RV_TO_MULT] += v[drop*V_DROP_MULT + inv*V_INV_MULT + from*V_FROM_MULT] *
                        mat[drop*MAT_DROP_MULT + inv*MAT_INV_MULT + from*MAT_FROM_MULT + to*MAT_TO_MULT];
                }
            }
        }
    }
}

void do_piecewise_matrix_multiply(double *v, double *mat, double *rv,
                                            int N_INVARIANT, int N_FROM, int N_TO, bool overwrite=true)
{
    int V_INV_MULT = 1;
    int V_FROM_MULT = V_INV_MULT * N_INVARIANT;

    int RV_INV_MULT = 1;
    int RV_TO_MULT = N_INVARIANT;

    int MAT_INV_MULT = 1;
    int MAT_FROM_MULT = MAT_INV_MULT * N_INVARIANT;
    int MAT_TO_MULT = MAT_FROM_MULT * N_FROM;

    if (overwrite)
        clear_arr(rv, N_INVARIANT * N_TO);

    for (int to=0; to<N_TO; to++)
    {
        for (int from=0; from<N_FROM; from++)
        {
            for (int inv=0; inv<N_INVARIANT; inv++)
            {
                rv[inv*RV_INV_MULT + to*RV_TO_MULT] += v[inv*V_INV_MULT + from*V_FROM_MULT] *
                    mat[inv*MAT_INV_MULT + from*MAT_FROM_MULT + to*MAT_TO_MULT];
            }
        }
    }
}

void do_matrix_multiply(double *v, double *mat, double *rv,
                        int N_FROM, int N_TO, bool overwrite=true)
{
    if (overwrite)
        clear_arr(rv, N_TO);

    for (int to=0; to<N_TO; to++)
    {
        for (int from=0; from<N_FROM; from++)
        {
            rv[to] += v[from] * mat[from + to*N_FROM];
        }
    }
}

void do_transition_multiply(double *transitions, double *population, double *rv,
                            const IntegerVector &pop_dims, int transition_dim, bool overwrite=false)
{
    int N_BEFORE = 1;
    int N_AFTER = 1;
    int N_TRANSITION = pop_dims[transition_dim];

    for (int i=0; i<transition_dim; i++)
        N_BEFORE *= pop_dims[i];

    for (int i=transition_dim+1; i<pop_dims.length(); i++)
        N_AFTER *= pop_dims[i];

//    int N_AFTER_TRANS = N_AFTER*N_TRANSITION;
    int N_BEFORE_TRANS = N_BEFORE * N_TRANSITION;
    int N_POP = N_AFTER*N_TRANSITION*N_BEFORE;

    if (overwrite)
        clear_arr(rv, N_POP);

    for (int after=0; after < N_AFTER; after++)
    {
        for (int transition_to=0; transition_to<N_TRANSITION; transition_to++)
        {
            for (int transition_from=0; transition_from<N_TRANSITION; transition_from++)
            {
                for (int before=0; before<N_BEFORE; before++)
                {
                    rv[before + transition_to*N_BEFORE + after*N_BEFORE_TRANS] +=
                    population[before + transition_from*N_BEFORE + after*N_BEFORE_TRANS] *
                    transitions[before + transition_from*N_BEFORE + after*N_BEFORE_TRANS + transition_to*N_POP];
                }
            }
        }
    }

}

//--------------------------//
//-- THE MAIN DX FUNCTION --//
//--------------------------//

// [[Rcpp::export]]
NumericVector dx_function(NumericVector hiv_positive,
                          NumericVector hiv_negative,
                          IntegerVector DIMS_HIV,
                          IntegerVector DIMS_NONHIV,
                          NumericVector AGING_RATES_HIV,
                          NumericVector AGING_RATES_NONHIV,
                          NumericVector FERTILITY_RATES_HIV,
                          NumericVector FERTILITY_RATES_NONHIV,
                          NumericVector FRACTION_BIRTHS_INFECTED,
                          LogicalVector BIRTHS_DIMMASK_NONHIV,
                          NumericVector BIRTH_PROPORTIONS_NONHIV,
                          LogicalVector BIRTHS_DIMMASK_HIV_TO_NONHIV,
                          NumericVector BIRTH_PROPORTIONS_HIV_TO_NONHIV,
                          LogicalVector BIRTHS_DIMMASK_HIV_TO_HIV,
                          NumericVector BIRTH_PROPORTIONS_HIV_TO_HIV,
                          int NUM_TRANSMISSION_ROUTES,
                          List CONTACT_MATRICES,
                          List TRANSMISSIBILITIES,
                          List SUSCEPTIBILITIES,
                          NumericVector GLOBAL_TRANSMISSION_RATES,
                          List CONTACT_DIMMASKS,
                          NumericVector NEW_INFECTION_PROPORTIONS,
                          IntegerVector COLLAPSE_INCIDENCE_INDICES,
                          List HIV_TRANSITION_ARRAYS,
                          IntegerVector HIV_TRANSITION_DIMENSIONS,
                          List NONHIV_TRANSITION_ARRAYS,
                          IntegerVector NONHIV_TRANSITION_DIMENSIONS,
                          IntegerVector TRACKED_HIV_TRANSITION_INDICES,
                          IntegerVector TRACKED_HIV_TRANSITION_FROM,
                          IntegerVector TRACKED_HIV_TRANSITION_ARRAY_INDICES,
                          IntegerVector TRACKED_HIV_TRANSITION_COLLAPSE_INDICES,
                          IntegerVector TRACKED_NONHIV_TRANSITION_INDICES,
                          IntegerVector TRACKED_NONHIV_TRANSITION_FROM,
                          IntegerVector TRACKED_NONHIV_TRANSITION_ARRAY_INDICES,
                          IntegerVector TRACKED_NONHIV_TRANSITION_COLLAPSE_INDICES,
                          NumericVector HIV_SPECIFIC_MORTALITY_RATES,
                          NumericVector GENERAL_MORTALITY_RATES_HIV,
                          NumericVector GENERAL_MORTALITY_RATES_NONHIV,
                          int FIXED_SIZE,
                          LogicalVector FIXED_SIZE_STRATA_DIMMASK,
                          bool MODEL_BIRTHS,
                          bool MODEL_MATERNAL_TRANSMISSION,
                          bool TRACK_HIV_SPECIFIC_MORTALITY,
                          LogicalVector COLLAPSE_DIMMASK_HIV_SPECIFIC_MORTALITY,
                          bool TRACK_HIV_OVERALL_MORTALITY,
                          LogicalVector COLLAPSE_DIMMASK_HIV_OVERALL_MORTALITY,
                          bool TRACK_NONHIV_MORTALITY,
                          LogicalVector COLLAPSE_DIMMASK_NONHIV_MORTALITY)
{

//----------------------------------//
//-- SET UP VARIABLES and STORAGE --//
//----------------------------------//

    //-- Declare Dimension and Size Variables --//
    IntegerVector DIMS_INCIDENT = IntegerVector::create(DIMS_HIV[AGE], DIMS_HIV[RACE], DIMS_HIV[SUBPOPULATION], DIMS_HIV[SEX], DIMS_HIV[RISK],
                                DIMS_NONHIV[NONHIV_SUBSET], DIMS_HIV[CONTINUUM], DIMS_HIV[CD4], DIMS_HIV[HIV_SUBSET]);
    IntegerVector DIMS_GENERAL = IntegerVector::create(DIMS_HIV[AGE], DIMS_HIV[RACE], DIMS_HIV[SUBPOPULATION], DIMS_HIV[SEX], DIMS_HIV[RISK]);
    int N_HIV = prod(DIMS_HIV);
    int N_NONHIV = prod(DIMS_NONHIV);
    int N_GENERAL = DIMS_HIV[AGE] * DIMS_HIV[RACE] * DIMS_HIV[SUBPOPULATION] * DIMS_HIV[SEX] * DIMS_HIV[RISK];
    int N_INCIDENT = max(COLLAPSE_INCIDENCE_INDICES)+1;//prod(DIMS_INCIDENT);

    int N_TRACKED_HIV_TRANSITIONS = 0;
    if (TRACKED_HIV_TRANSITION_COLLAPSE_INDICES.size() > 0)
        N_TRACKED_HIV_TRANSITIONS = max(TRACKED_HIV_TRANSITION_COLLAPSE_INDICES) + 1;

    int N_TRACKED_NONHIV_TRANSITIONS = 0;
    if (TRACKED_NONHIV_TRANSITION_COLLAPSE_INDICES.size() > 0)
        N_TRACKED_NONHIV_TRANSITIONS = max(TRACKED_NONHIV_TRANSITION_COLLAPSE_INDICES) + 1;

    int N_MAX = N_HIV;
    if (N_NONHIV > N_HIV)
        N_MAX = N_NONHIV;

    int N_HIV_SPECIFIC_MORTALITY = prod_if(DIMS_HIV, COLLAPSE_DIMMASK_HIV_SPECIFIC_MORTALITY);
    int N_HIV_OVERALL_MORTALITY = prod_if(DIMS_HIV, COLLAPSE_DIMMASK_HIV_OVERALL_MORTALITY);
    int N_NONHIV_MORTALITY = prod_if(DIMS_NONHIV, COLLAPSE_DIMMASK_NONHIV_MORTALITY);

    //-- Allocate scratch arrays --//
    int N_BIG_SCRATCH = 4*N_MAX;
    double big_scratch[N_BIG_SCRATCH];

    double *scratch1 = big_scratch;
    double *scratch2 = scratch1 + N_MAX;
    double *scratch3 = scratch2 + N_MAX;
    double *scratch4 = scratch3 + N_MAX;

    //-- Offsets/Indices for Storage into the Delta array --//
    int HIV_OFFSET = 0;
    int NONHIV_OFFSET = N_HIV;
    int INCIDENCE_OFFSET = NONHIV_OFFSET + N_NONHIV;
    int n_track = INCIDENCE_OFFSET + N_INCIDENT;

    int HIV_SPECIFIC_MORTALITY_OFFSET = 0;
    if (TRACK_HIV_SPECIFIC_MORTALITY)
    {
        HIV_SPECIFIC_MORTALITY_OFFSET = n_track;
        n_track += N_HIV_SPECIFIC_MORTALITY;
    }

    int HIV_OVERALL_MORTALITY_OFFSET = 0;
    if (TRACK_HIV_OVERALL_MORTALITY)
    {
        HIV_OVERALL_MORTALITY_OFFSET = n_track;
        n_track += N_HIV_OVERALL_MORTALITY;
    }

    int NONHIV_MORTALITY_OFFSET = 0;
    if (TRACK_NONHIV_MORTALITY)
    {
        NONHIV_MORTALITY_OFFSET = n_track;
        n_track += N_NONHIV_MORTALITY;
    }

    int HIV_TRACKED_TRANSITION_OFFSET = n_track;
    n_track += N_TRACKED_HIV_TRANSITIONS; //TRACKED_HIV_TRANSITION_INDICES.size();

    int NONHIV_TRACKED_TRANSITION_OFFSET = n_track;
    n_track += N_TRACKED_NONHIV_TRANSITIONS; //TRACKED_NONHIV_TRANSITION_INDICES.size();

    int N_AGE = DIMS_HIV[AGE];
    int N_NONAGE_HIV = N_HIV / N_AGE;
    int N_NONAGE_NONHIV = N_NONHIV / N_AGE;
    int N_RACE = DIMS_HIV[RACE];

    //-- Allocate the return vector (delta) and index its components --//
    NumericVector delta(n_track);
    clear_arr(delta.begin(), n_track);

    double *delta_hiv = delta.begin() + HIV_OFFSET;
    double *delta_nonhiv = delta.begin() + NONHIV_OFFSET;
    double *incident_cases = delta.begin() + INCIDENCE_OFFSET;

    double *rv_hiv_specific_mortality;
    if (TRACK_HIV_SPECIFIC_MORTALITY)
        rv_hiv_specific_mortality = delta.begin() + HIV_SPECIFIC_MORTALITY_OFFSET;

    double *rv_hiv_overall_mortality;
    if (TRACK_HIV_OVERALL_MORTALITY)
        rv_hiv_overall_mortality = delta.begin() + HIV_OVERALL_MORTALITY_OFFSET;

    double *rv_nonhiv_mortality;
    if (TRACK_NONHIV_MORTALITY)
        rv_nonhiv_mortality = delta.begin() + NONHIV_MORTALITY_OFFSET;

    double *tracked_transitions_hiv = delta.begin() + HIV_TRACKED_TRANSITION_OFFSET;
    double *tracked_transitions_nonhiv = delta.begin() + NONHIV_TRACKED_TRANSITION_OFFSET;

    //-- Sum the population --//
    double general_population[N_GENERAL];
    LogicalVector GENERAL_FROM_NONHIV_MASK = LogicalVector::create(1,1,1,1,1,0);
    LogicalVector GENERAL_FROM_HIV_MASK = LogicalVector::create(1,1,1,1,1,0,0,0);
    collapse_or_expand_8d(hiv_positive.begin(), general_population, DIMS_HIV, GENERAL_FROM_HIV_MASK, false, true);
    collapse_or_expand_6d(hiv_negative.begin(), general_population, DIMS_NONHIV, GENERAL_FROM_NONHIV_MASK, false, false);


    //------------//
    //-- BIRTHS --//
    //------------//

    //-- Births --//
    if (MODEL_BIRTHS)
    {
        //Calculate non-hiv births
        double *births_to_nonhiv = scratch1;
        double *births_from_nonhiv = scratch2;
        double *births_from_nonhiv_summed = scratch3;
        int n_age_for_births = 1;
        if (BIRTHS_DIMMASK_NONHIV[AGE])
            n_age_for_births = N_AGE;
        int n_race_for_births = 1;
        if (BIRTHS_DIMMASK_NONHIV[RACE])
            n_race_for_births = N_RACE;

        multiply_arr(hiv_negative, FERTILITY_RATES_NONHIV, births_from_nonhiv);
        collapse_or_expand_6d(births_from_nonhiv, births_from_nonhiv_summed, DIMS_NONHIV, BIRTHS_DIMMASK_NONHIV, false, true);
        do_piecewise_matrix_multiply_with_drop(births_from_nonhiv_summed, BIRTH_PROPORTIONS_NONHIV.begin(), births_to_nonhiv,
                                               n_age_for_births, n_race_for_births,
                                               prod_if(DIMS_NONHIV, BIRTHS_DIMMASK_NONHIV)/n_age_for_births/n_race_for_births,
                                               N_NONAGE_NONHIV/n_race_for_births,
                                               true); //the last argument tells it to overwrite
        //NB: births_to_nonhiv is in scratch1

        //Calculate hiv births to non-hiv
        double *births_from_hiv = scratch2;
        double *births_hiv_to_nonhiv = scratch3;
        double *births_hiv_to_nonhiv_summed = scratch4;
        n_age_for_births = 1;
        if (BIRTHS_DIMMASK_HIV_TO_NONHIV[AGE])
            n_age_for_births = N_AGE;
        n_race_for_births = 1;
        if (BIRTHS_DIMMASK_HIV_TO_NONHIV[RACE])
            n_race_for_births = N_RACE;

        multiply_arr(hiv_positive, FERTILITY_RATES_HIV, births_from_hiv);
        if (MODEL_MATERNAL_TRANSMISSION)
            multiply_arr(births_from_hiv, FRACTION_BIRTHS_INFECTED, births_hiv_to_nonhiv);
        else
            births_hiv_to_nonhiv = births_from_hiv;
        collapse_or_expand_8d(births_hiv_to_nonhiv, births_hiv_to_nonhiv_summed, DIMS_HIV, BIRTHS_DIMMASK_HIV_TO_NONHIV, false, true);
        do_piecewise_matrix_multiply_with_drop(births_hiv_to_nonhiv_summed, BIRTH_PROPORTIONS_HIV_TO_NONHIV.begin(), births_to_nonhiv,
                                               n_age_for_births, n_race_for_births,
                                               prod_if(DIMS_HIV, BIRTHS_DIMMASK_HIV_TO_NONHIV)/n_age_for_births/n_race_for_births,
                                               N_NONAGE_NONHIV/n_race_for_births,
                                               false); //the last argument tells it NOT to overwrite
        //NB: births_to_nonhiv is in scratch1
        //    births_from_hiv is in scratch2 ONLY IF modeling maternal transmission

        //Fold in all non-hiv births
        for (int i=0; i<N_NONAGE_NONHIV; i++)
            delta_nonhiv[i*N_AGE] += births_to_nonhiv[i];

        // Calculate HIV-infected births
        if (MODEL_MATERNAL_TRANSMISSION)
        {
            double *births_hiv_to_hiv_summed = scratch3;
            double *births_to_hiv = scratch4;
            n_age_for_births = 1;
            if (BIRTHS_DIMMASK_HIV_TO_HIV[AGE])
                n_age_for_births = N_AGE;
            n_race_for_births = 1;
            if (BIRTHS_DIMMASK_HIV_TO_HIV[RACE])
                n_race_for_births = N_RACE;

            decrement_arr(births_from_hiv, births_hiv_to_nonhiv, N_HIV);
            collapse_or_expand_8d(births_from_hiv, births_hiv_to_hiv_summed, DIMS_HIV, BIRTHS_DIMMASK_HIV_TO_HIV, false, true);
            do_piecewise_matrix_multiply_with_drop(births_hiv_to_hiv_summed, BIRTH_PROPORTIONS_HIV_TO_HIV.begin(), births_to_hiv,
                                                   n_age_for_births, n_race_for_births,
                                                   prod_if(DIMS_HIV, BIRTHS_DIMMASK_HIV_TO_HIV)/n_age_for_births/n_race_for_births,
                                                   N_HIV/n_race_for_births,
                                                   true); //the last argument tells it to overwrite

            //Fold in all hiv births
            for (int i=0; i<N_NONAGE_HIV; i++)
                delta_hiv[i*N_AGE] += births_to_nonhiv[i];
        }
    }

    //---------------//
    //-- MORTALITY --//
    //---------------//

    double *hiv_overall_mortality = scratch1;
    double *hiv_specific_mortality = scratch2;
    double *nonhiv_mortality = scratch3;

    //-- Calculate Mortality --//
    multiply_arr(hiv_positive, HIV_SPECIFIC_MORTALITY_RATES, hiv_specific_mortality);

    multiply_arr(hiv_positive, GENERAL_MORTALITY_RATES_HIV, hiv_overall_mortality);
    increment_arr(hiv_overall_mortality, hiv_specific_mortality, N_HIV);

    multiply_arr(hiv_negative, GENERAL_MORTALITY_RATES_NONHIV, nonhiv_mortality);

    //-- Fold into delta_hiv and delta_nonhiv --//
    decrement_arr(delta_hiv, hiv_overall_mortality, N_HIV);
    decrement_arr(delta_nonhiv, nonhiv_mortality, N_NONHIV);


    //-----------//
    //-- AGING --//
    //-----------//

    //-- Aging HIV-positive --//
    double *aging_up = scratch4; //need to keep this distinct from hiv_overall_mortality (scratch1) and nonhiv_mortality (scratch3) above
    multiply_arr(hiv_positive, AGING_RATES_HIV, aging_up);
    decrement_arr(delta_hiv, aging_up, N_HIV);
    for (int nonage=0; nonage<N_NONAGE_HIV; nonage++)
    {
        for (int age=1; age<N_AGE; age++)
            delta_hiv[nonage*N_AGE + age] += aging_up[nonage*N_AGE + (age-1)];
    }

    if (TRACK_HIV_OVERALL_MORTALITY)
    {
        for (int nonage=0; nonage<N_NONAGE_HIV; nonage++)
        {
            hiv_overall_mortality[nonage*N_AGE + N_AGE-1] += aging_up[nonage*N_AGE + N_AGE-1];
        }
    }

    //-- Aging HIV-negative --//
    multiply_arr(hiv_negative, AGING_RATES_NONHIV, aging_up);


    decrement_arr(delta_nonhiv, aging_up, N_NONHIV);
    for (int nonage=0; nonage<N_NONAGE_NONHIV; nonage++)
    {
        for (int age=1; age<N_AGE; age++)
            delta_nonhiv[nonage*N_AGE + age] += aging_up[nonage*N_AGE + (age-1)];
    }

    if (TRACK_NONHIV_MORTALITY)
    {
        for (int nonage=0; nonage<N_NONAGE_NONHIV; nonage++)
        {
            nonhiv_mortality[nonage*N_AGE + N_AGE-1] += aging_up[nonage*N_AGE + N_AGE-1];
        }
    }

    //-- Collapse and save mortality if warranted (now that we have folded in aging out)--//
    // NB: these were cleared (as components of delta) up at the top
    if (TRACK_HIV_SPECIFIC_MORTALITY)
        collapse_or_expand_8d(hiv_specific_mortality, rv_hiv_specific_mortality,
                              DIMS_HIV, COLLAPSE_DIMMASK_HIV_SPECIFIC_MORTALITY,
                              false, false);

    if (TRACK_HIV_OVERALL_MORTALITY)
        collapse_or_expand_8d(hiv_overall_mortality, rv_hiv_overall_mortality,
                              DIMS_HIV, COLLAPSE_DIMMASK_HIV_OVERALL_MORTALITY,
                              false, false);

    if (TRACK_NONHIV_MORTALITY)
        collapse_or_expand_6d(nonhiv_mortality, rv_nonhiv_mortality,
                              DIMS_NONHIV, COLLAPSE_DIMMASK_NONHIV_MORTALITY,
                              false, false);


    //----------------------//
    //-- HIV TRANSMISSION --//
    //----------------------//

    double *infection_hazards = scratch1;
    clear_arr(infection_hazards, N_NONHIV);

    double *transmissibility = scratch3;
    double *force_of_infection = scratch4;
    double *denominator = scratch3;
    double *one_hazard = scratch3;

    for (int route=0; route<NUM_TRANSMISSION_ROUTES; route++)
    {
        LogicalVector route_dimmask_general = (LogicalVector) CONTACT_DIMMASKS[route];
        LogicalVector route_dimmask_hiv = LogicalVector::create(route_dimmask_general[0],route_dimmask_general[1],
                                                                route_dimmask_general[2],route_dimmask_general[3],
                                                                route_dimmask_general[4],0,0,0);
        int ndim_route = sum(route_dimmask_general);
        int N_route = prod_if(DIMS_GENERAL, route_dimmask_general);

        multiply_arr(((NumericVector)TRANSMISSIBILITIES[route]).begin(), hiv_positive, transmissibility);
        collapse_or_expand_8d(transmissibility, force_of_infection, DIMS_HIV, route_dimmask_hiv, false, true);

        collapse_or_expand_5d(general_population, denominator, DIMS_GENERAL, route_dimmask_general, false, true);

        divide_arr_zero_if_zero_denom(force_of_infection, denominator, force_of_infection, N_GENERAL); //in scratch4

        do_matrix_multiply(force_of_infection, ((NumericVector)CONTACT_MATRICES[route]).begin(), one_hazard,
                           N_route, N_NONHIV, true);

        NumericVector susceptibility = ((NumericVector)SUSCEPTIBILITIES[route]);

        for (int i=0; i<N_NONHIV; i++)
            infection_hazards[i] += one_hazard[i] * susceptibility[i] * GLOBAL_TRANSMISSION_RATES[route];

    }
    double *incident_cases_from = infection_hazards; //in scratch1
    multiply_arr(infection_hazards, hiv_negative, incident_cases_from);

    //Split the incident cases by proportions
    double *incident_cases_to = scratch2;
    clear_arr(incident_cases_to, N_HIV);

    IntegerVector from_mult = cumprod_rightshifted(DIMS_NONHIV);
    IntegerVector to_mult = cumprod_rightshifted(DIMS_HIV);
    IntegerVector all_mult = cumprod_rightshifted(DIMS_INCIDENT);

    int from_index,to_index,all_index;
    for (int hsub=0; hsub<DIMS_HIV[HIV_SUBSET]; hsub++)
    {
        for (int cd4=0; cd4<DIMS_HIV[CD4]; cd4++)
        {
            for (int cont=0; cont<DIMS_HIV[CONTINUUM]; cont++)
            {
                for (int nhsub=0; nhsub<DIMS_NONHIV[NONHIV_SUBSET]; nhsub++)
                {
                    for (int risk=0; risk<DIMS_HIV[RISK]; risk++)
                    {
                        for (int sex=0; sex<DIMS_HIV[SEX]; sex++)
                        {
                            for (int sub=0; sub<DIMS_HIV[SUBPOPULATION]; sub++)
                            {
                                for (int race=0; race<DIMS_HIV[RACE]; race++)
                                {
                                    for (int age=0; age<DIMS_HIV[AGE]; age++)
                                    {
                                        from_index = age*from_mult[AGE] + race*from_mult[RACE] +
                                            sub*from_mult[SUBPOPULATION] + sex*from_mult[SEX] +
                                            risk*from_mult[RISK] + nhsub*from_mult[NONHIV_SUBSET];

                                        to_index = age*to_mult[AGE] + race*to_mult[RACE] +
                                            sub*to_mult[SUBPOPULATION] + sex*to_mult[SEX] +
                                            risk*to_mult[RISK] + cont*to_mult[CONTINUUM] +
                                            cd4*to_mult[CD4] + hsub*to_mult[HIV_SUBSET];

                                        all_index = age*all_mult[AGE] + race*all_mult[RACE] +
                                            sub*all_mult[SUBPOPULATION] + sex*all_mult[SEX] +
                                            risk*all_mult[RISK] + nhsub*all_mult[NONHIV_SUBSET] +
                                            cont*all_mult[CONTINUUM+1] +
                                            cd4*all_mult[CD4+1] + hsub*all_mult[HIV_SUBSET+1];

                                        incident_cases[COLLAPSE_INCIDENCE_INDICES[all_index]] += incident_cases_from[from_index] * NEW_INFECTION_PROPORTIONS[all_index];
                                        incident_cases_to[to_index] += incident_cases_from[from_index] * NEW_INFECTION_PROPORTIONS[all_index];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    decrement_arr(delta_nonhiv, incident_cases_from, N_NONHIV);
    increment_arr(delta_hiv, incident_cases_to, N_HIV);

    //-----------------//
    //-- TRANSITIONS --//
    //-----------------//

    double *transition_array;

    //-- HIV-positive Transitions --//
    double *hiv_transitions_to = scratch1;
    double *rowsummed_hiv_transition_array = scratch2;
    double *hiv_transitions_from = scratch3;

    clear_arr(hiv_transitions_to, N_HIV);
    clear_arr(rowsummed_hiv_transition_array, N_HIV);
    for (int i=0; i<HIV_TRANSITION_ARRAYS.length(); i++)
    {
        transition_array = ((NumericVector) HIV_TRANSITION_ARRAYS[i]).begin();
        do_transition_multiply(transition_array, hiv_positive.begin(), hiv_transitions_to,
                               DIMS_HIV, HIV_TRANSITION_DIMENSIONS[i], false);

        rowsum(transition_array, rowsummed_hiv_transition_array,
               N_HIV, DIMS_HIV[HIV_TRANSITION_DIMENSIONS[i]], false);
    }

    multiply_arr(rowsummed_hiv_transition_array, hiv_positive, hiv_transitions_from);

    increment_arr(delta_hiv, hiv_transitions_to, N_HIV);
    decrement_arr(delta_hiv, hiv_transitions_from, N_HIV);

    //-- HIV-negative Transitions --//
    double *nonhiv_transitions_to = scratch1;
    double *rowsummed_nonhiv_transition_array = scratch2;
    double *nonhiv_transitions_from = scratch3;


    clear_arr(hiv_transitions_to, N_NONHIV);
    clear_arr(rowsummed_nonhiv_transition_array, N_NONHIV);
    for (int i=0; i<NONHIV_TRANSITION_ARRAYS.length(); i++)
    {
        transition_array = ((NumericVector) NONHIV_TRANSITION_ARRAYS[i]).begin();
        do_transition_multiply(transition_array, hiv_negative.begin(), nonhiv_transitions_to,
                               DIMS_NONHIV, NONHIV_TRANSITION_DIMENSIONS[i], false);

        rowsum(transition_array, rowsummed_nonhiv_transition_array,
               N_NONHIV, DIMS_NONHIV[NONHIV_TRANSITION_DIMENSIONS[i]], false);
    }

    multiply_arr(rowsummed_nonhiv_transition_array, hiv_negative, nonhiv_transitions_from);

    increment_arr(delta_nonhiv, nonhiv_transitions_to, N_NONHIV);
    decrement_arr(delta_nonhiv, nonhiv_transitions_from, N_NONHIV);

    //-- Tracked Transitions --//
    // NB: the tracked transitions were cleared up top
    for (int i=0; i<TRACKED_HIV_TRANSITION_INDICES.size(); i++)
        tracked_transitions_hiv[TRACKED_HIV_TRANSITION_COLLAPSE_INDICES[i]] += hiv_positive[TRACKED_HIV_TRANSITION_FROM[i]] *
            ((NumericVector) HIV_TRANSITION_ARRAYS[TRACKED_HIV_TRANSITION_ARRAY_INDICES[i]])[TRACKED_HIV_TRANSITION_INDICES[i]];

    for (int i=0; i<TRACKED_NONHIV_TRANSITION_INDICES.size(); i++)
        tracked_transitions_nonhiv[TRACKED_NONHIV_TRANSITION_COLLAPSE_INDICES[i]] += hiv_negative[TRACKED_NONHIV_TRANSITION_FROM[i]] *
            ((NumericVector) NONHIV_TRANSITION_ARRAYS[TRACKED_NONHIV_TRANSITION_ARRAY_INDICES[i]])[TRACKED_NONHIV_TRANSITION_INDICES[i]];



    //---------------------------//
    //-- CONSTANT STRATA SIZES --//
    //---------------------------//

    if (FIXED_SIZE==FIXED_SIZE_BY_STRATA)
    {
        // To keep strata size constant, we need to calculate, for each fixed stratum j, the offsets such that
        //   delta_j - offset_j = 0
        // The offset for stratum i, which belongs to superstratum j with fixed size is
        // offset_i = delta_j * (pre-delta_i + delta_i) / (pre-delta_j + delta_j)

        int N_FIXED_STRATA = prod_if(DIMS_GENERAL, FIXED_SIZE_STRATA_DIMMASK);

        LogicalVector FIXED_STRATA_FROM_HIV_DIMMASK = LogicalVector::create(FIXED_SIZE_STRATA_DIMMASK[0],FIXED_SIZE_STRATA_DIMMASK[1],
                                                                            FIXED_SIZE_STRATA_DIMMASK[2],FIXED_SIZE_STRATA_DIMMASK[3],
                                                                            FIXED_SIZE_STRATA_DIMMASK[4],0,0,0);
        LogicalVector FIXED_STRATA_FROM_NONHIV_DIMMASK = LogicalVector::create(FIXED_SIZE_STRATA_DIMMASK[0],FIXED_SIZE_STRATA_DIMMASK[1],
                                                                               FIXED_SIZE_STRATA_DIMMASK[2],FIXED_SIZE_STRATA_DIMMASK[3],
                                                                               FIXED_SIZE_STRATA_DIMMASK[4],0);
        double *net_change_by_stratum = scratch1;
        double *strata_after_change = scratch2;
        double *offset_multiplier_by_stratum = net_change_by_stratum;

        double *offset_multiplier_hiv = scratch3;
        double *hiv_after_delta = scratch4;
        double *offset_hiv = hiv_after_delta;

        double *offset_multiplier_nonhiv = scratch3;
        double *nonhiv_after_delta = scratch4;
        double *offset_nonhiv = nonhiv_after_delta;

        // Calculate summed delta for each fixed stratum
        collapse_or_expand_8d(delta_hiv, net_change_by_stratum, DIMS_HIV, FIXED_STRATA_FROM_HIV_DIMMASK, false, true);
        collapse_or_expand_6d(delta_nonhiv, net_change_by_stratum, DIMS_NONHIV, FIXED_STRATA_FROM_NONHIV_DIMMASK, false, false);
        // now net_change_by_stratum = delta_j (ie, sum of all deltas within stratum j)

        // Calculate the offset multiplier for each fixed stratum
        collapse_or_expand_5d(general_population, strata_after_change, DIMS_GENERAL, FIXED_SIZE_STRATA_DIMMASK, false, true);
        increment_arr(strata_after_change, net_change_by_stratum, N_FIXED_STRATA);
        divide_arr_zero_if_zero_denom(net_change_by_stratum, strata_after_change, offset_multiplier_by_stratum, N_FIXED_STRATA);
        // now offset_by_stratum = delta_j / (pre-delta_j + delta_j)

        // Fold the offset into HIV+
        collapse_or_expand_8d(offset_multiplier_by_stratum, offset_multiplier_hiv, DIMS_HIV, FIXED_STRATA_FROM_HIV_DIMMASK, true, true);
        add_arr(hiv_positive, delta_hiv, hiv_after_delta);
        multiply_arr(hiv_after_delta, offset_multiplier_hiv, offset_hiv, N_HIV);
        decrement_arr(delta_hiv, offset_hiv, N_HIV);

        // Fold the offset into HIV-
        collapse_or_expand_6d(offset_multiplier_by_stratum, offset_multiplier_nonhiv, DIMS_NONHIV, FIXED_STRATA_FROM_NONHIV_DIMMASK, true, true);
        add_arr(hiv_negative, delta_nonhiv, nonhiv_after_delta);
        multiply_arr(nonhiv_after_delta, offset_multiplier_nonhiv, offset_nonhiv, N_NONHIV);
        decrement_arr(delta_nonhiv, offset_nonhiv, N_NONHIV);
    }

    //------------//
    //-- RETURN --//
    //------------//

    return delta;
}

// [[Rcpp::export]]
NumericVector test(List l)
{
    Rcout << l.size() << "\n";

    double *test = ((NumericVector)l[0]).begin();
    Rcout << test[35] << "\n";

    NumericVector v1 = l[0];
//    Rcout << v1 << "\n";

    NumericVector v2 = l[1];
 //   Rcout << v2 << "\n";

    return v2;
}

//----------------------------//
//-- OTHER EXPORTED HELPERS --//
//----------------------------//

// [[Rcpp::export]]
NumericVector do_marginal_sums_hiv_positive(NumericVector arr, IntegerVector dims, LogicalVector keep_mask)
{
    int N = prod_if(dims, keep_mask);
    NumericVector rv(N);

    collapse_or_expand_8d(arr.begin(), rv.begin(), dims, keep_mask);

    return rv;
}

// [[Rcpp::export]]
NumericVector do_marginal_sums_hiv_negative(NumericVector arr, IntegerVector dims, LogicalVector keep_mask)
{
    int N = prod_if(dims, keep_mask);
    NumericVector rv(N);

    collapse_or_expand_6d(arr.begin(), rv.begin(), dims, keep_mask);

    return rv;
}

// [[Rcpp::export]]
NumericVector do_marginal_sums_general(NumericVector arr, IntegerVector dims, LogicalVector keep_mask)
{
    int N = prod_if(dims, keep_mask);
    NumericVector rv(N);

    collapse_or_expand_6d(arr.begin(), rv.begin(), dims, keep_mask);

    return rv;
}
