#include <Rcpp.h>
using namespace Rcpp;





//' Helper function for from_ranks_to_integrable_values
//'
//' @param rank_interval_mult the value of the normalized density of Y_A.
//' @param j_max maximum number of sample points.
//'
// [[Rcpp::export]]
NumericVector cpp_helper_from_ranks_to_integrable_values(NumericVector rank_interval_mult, int j_max)
{
  NumericVector j_vec(j_max +1);


  double p_corresponding_to_j;
  int r_last = 0;
  int r_current = 0;
  double value_in_j_vec = rank_interval_mult[0];
  int r_max = rank_interval_mult.length()-1;

  // j goes from 0 to j_max
  for (int j = 0; j <= j_max; j++)
  {
    p_corresponding_to_j = (double) j / (double) j_max;

    // Find which rank corresponds to position j in j_vec
    for (int k = r_last; k <= r_max; k++)
    {
      if ((double) k / (double) (r_max+1) <= p_corresponding_to_j) {
        r_current = k;
      }else{
        break;
      }
    }



    // measure the amount increased in cumulative, or the actual density in non cumulative
    // only if the r_current was updated.
    if(r_current != r_last)
    {
      double sumOfRanks = 0;
      for (int k = r_last+1; k <= r_current; k++)
      {
        sumOfRanks += rank_interval_mult[k];
      }


      if(sumOfRanks == 0)
      {
        value_in_j_vec = 0;
      }
      else
      {
        value_in_j_vec = sumOfRanks / (r_current - r_last);
      }

    }

    j_vec[j] = value_in_j_vec;
    r_last = r_current;
  }

  return(j_vec);
}
