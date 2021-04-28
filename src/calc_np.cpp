#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]

double fcalc_np(double c_i, NumericVector &c_j) {
  //initialize c_j_sum
  double c_j_sum = 0;
  int n = c_j.size();

  // calculate the sum of c_j
  for(int i = 0; i < n; ++i) {
    c_j_sum += c_j[i];
  }
  double g_i = c_i * log(c_i / (c_j_sum));
  return g_i;
}


// [[Rcpp::export]]

NumericVector fcalc_np_all(List neighbors, StringVector v,
                           NumericVector exp) {

  int n = v.size(); // define size of node list
  NumericVector np_vec(n); np_vec.names() = v; //initialize np_vec

  // initialize the variables we will use in the loop- does this save time?
  String v_i;
  NumericVector neighbors_i;
  StringVector neighbors_named_i;
  NumericVector c_j;
  double c_i;
  int neighbors_n;

  for(int i = 0; i < n; ++i) {

    // do all the indexing to create a named list of neighbors for a given node v
    // then created named NumericVector with the expression for those neighbors
    // guarantee this could be done in fewer lines but oh well

    v_i = v[i];
    neighbors_i = neighbors[v_i];
    neighbors_n = neighbors_i.size();

    // need to correct for zero indexing in cpp vs 1 indexing in R
    for(int j = 0; j< neighbors_n; ++j) {
      neighbors_i[j] = neighbors_i[j] - 1;
    }
    neighbors_named_i = v[neighbors_i];

    //define c_j
    c_j = exp[neighbors_named_i];

    // define c_i
    c_i = exp[v_i];

    //calculate network potential for each node and assign to new vector
    np_vec[v_i] = fcalc_np(c_i = c_i, c_j = c_j);

  }

  return np_vec;

}

