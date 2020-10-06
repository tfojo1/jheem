
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List pull_time_varying_parameters(List param_sets, double time)
{
    List rv(param_sets.length());

    List param_set, values;
    NumericVector elem, times, elem_before, elem_after;
    int index_before, index_after;
    bool stepwise;
    double weight_before, weight_after;

    for (int i=0; i<param_sets.length(); i++)
    {
        param_set = (List) param_sets[i];
        times = as<NumericVector>(param_set["times"]);
        values = as<List>(param_set["values"]);
        stepwise = as<bool>(param_set["stepwise"]);
        int n_times = times.length();

        if (n_times==1)
            elem = (NumericVector) values[0];
        else
        {
            index_before = -1; //the last index such that times[t] <= time
            for (int t=0; t<n_times; t++)
            {
                if (times[t] <= time)
                    index_before = t;
                else
                    break;
            }

            if (index_before==n_times-1)
                index_after = -1;
            else
                index_after = index_before + 1; //the first index such that times[t] > time

            if (index_before==-1) //there were no times before
                elem = (NumericVector) values[index_after];
            else if (times[index_before]==time || stepwise || index_after==-1 ||
                     times[index_before]==R_PosInf || times[index_before]==R_NegInf)
                elem = (NumericVector) values[index_before];
            else if (times[index_after]==time || times[index_after]==R_PosInf || times[index_after]==R_NegInf)
                elem = (NumericVector) values[index_after];
            else
            {
                elem_before = (NumericVector) values[index_before];
                elem_after = (NumericVector) values[index_after];
                elem = NumericVector(elem_before.length());

                weight_before = (times[index_after] - time) / (times[index_after] - times[index_before]);
                weight_after = (time - times[index_before]) / (times[index_after] - times[index_before]);

                for (int j=0; j<elem_before.length(); j++)
                    elem[j] = elem_before[j] * weight_before + elem_after[j] * weight_after;
            }
        }

        rv[i] = elem;
    }

    return (rv);
}

// [[Rcpp::export]]
int test_pull_elem_a(List l)
{
    int rv = (int) l["a"];

    return (rv);
}

// [[Rcpp::export]]
void print_vector(NumericVector v)
{
    for (int i=0; i<v.length(); i++)
    {
        Rcout << v[i];
        if (v[i]==R_PosInf)
            Rcout << " POS INF ";
        if (v[i]==R_NegInf)
            Rcout << " NEG INF ";
        Rcout << "\n";
    }
}
