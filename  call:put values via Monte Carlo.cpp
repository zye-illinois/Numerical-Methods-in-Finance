//
//  main.cpp
//  IE525_hw1
//
//  Created by Zhiyi Ye on 3/27/21.
//
#include <algorithm>    // Needed for the "max" function
#include <cmath>
#include <iostream>
#include<time.h>
using namespace std;
using namespace std::chrono;

// A simple implementation of the Box-Muller algorithm, used to generate gaussian random numbers
double gaussian_box_muller() {
  double x = 0.0;
  double y = 0.0;
  double euclid_sq = 0.0;

  // Continue generating two uniform random variables
  // until the square of their "euclidean distance"
  // is less than unity
  do {
    x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    euclid_sq = x*x + y*y;
  } while (euclid_sq >= 1.0);

  return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}
int num_sims = 100000;
// Pricing a European vanilla call option with a Monte Carlo method
double monte_carlo_call_price(const int& num_sims, const double& S, const double& K, const double& r, double& q, const double& sigma,const double& T) {
  double S_adjust = S* exp(T*(r-q-0.5*sigma*sigma));
// current asset price
  double S_cur = 0.0;
  double payoff_sum = 0.0;
    double payoff_avg = 0.0;

  for (int i=0; i<num_sims; i++) {
    double gauss_bm = gaussian_box_muller();
    S_cur = S_adjust * exp(sqrt(sigma*sigma*T)*gauss_bm);
    payoff_sum += std::max(S_cur - K, 0.0);
    payoff_avg = payoff_sum/num_sims*exp(-r*T);
  }

  return (payoff_avg);
}

float standard_error(const int& num_sims, const double& S, const double& K, const double& r, double& q, const double& sigma,const double& T) {
    double S_adjust = S* exp(T*(r-q-0.5*sigma*sigma));
  // current asset price
    double S_cur = 0.0;
    double payoff_sum = 0.0;
    double payoff_avg = 0.0;
    double se=0;
    double pay_off =0;
    for (int i=0; i<num_sims; i++) {
      double gauss_bm = gaussian_box_muller();
      S_cur = S_adjust * exp(sqrt(sigma*sigma*T)*gauss_bm);
      pay_off = std::max(S_cur-K,0.0);
      payoff_sum+= pay_off;
      payoff_avg = payoff_sum/num_sims;
     se = sqrt((pay_off-payoff_avg)*(pay_off-payoff_avg)/(i-1));
        }
    return(se);
}

int main(int argc, char **argv) {
  // First we create the parameter list
    clock_t t1,t2;
    t1=clock();
  double S = 1868.99;  // Current value of price
  double K = 1870.0;  // Strike price
  double r = 0.003866;   // Risk-free rate
  double sigma = 0.2979;    // Volatility of the underlying
  double T = 0.01923;
  double q = 0.0232;
  double z = 1.96;
  double ci = z*sigma/num_sims;
  double call = monte_carlo_call_price(num_sims,S, K, r, q,sigma, T);
  double se =standard_error(num_sims,S, K, r, q,sigma, T);
  double low = call-ci;
  double high = call+ci;
  // Then we calculate the call/put values via Monte Carlo
  // Finally we output the parameters and prices
  std::cout << "Number of Paths: " << num_sims << std::endl;
  std::cout << "The estimated Call Price is: $" << call << std::endl;
  std::cout<<"The estimated standard error is: "<< se << std::endl;
  std::cout<<"The 95% confidence interval is: [" <<low<< " "<<high<<"]."<<std::endl;
  t2=clock();
  float seconds = ((float)t2-float(t1))/CLOCKS_PER_SEC;
  cout<<"This program took me: "<<seconds<<" seconds to run"<<endl;
  return 0;
    }

