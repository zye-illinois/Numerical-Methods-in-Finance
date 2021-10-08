// QMC Project 1_2
//  Created by zhiyi Ye on 4/30/21.

#include "sobol.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;
using namespace std::chrono;

double se;
double expc; //expectation
double K=100.0;
double T=1.00;
double r=0.1;
double q=0;
double sigma=0.2;
double S0=100;
int m; // Monitoring
double inv_m ;
double dT ;
double UNIRAN()
{
    double a;
    a = ((double)(rand()) + 1.) / ((double)(RAND_MAX)+1.);
    if (abs(a - 1.0) < 0.0000000001)
    {
        a = UNIRAN();
    }
    return a;
}
double NORMRAN()
{
    double u1 = UNIRAN();
    return cos(8.*atan(1.)*u1)*sqrt(-2.*log(u1));
}
double gaussian_box_muller() {
  double x = 0.0;
  double y = 0.0;
  double euclid_sq = 0.0;
  do {
    x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    euclid_sq = x*x + y*y;
  } while (euclid_sq >= 1.0);

  return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}
// QMC with sobol sequences
//Max function
double Max(double c, double d) {
    return (d < c )? c:d;
}
//Minimum function
double Min(double c, double d)
{
    return (d < c )? d:c;
}
double get_cdf(double p)
{
    if (p>15) {
        return 1;
    }
    if (p<-15) {
        return 0;
    }
    double A[15]=
    {1.253314137315500,0.6556795424187985,0.4213692292880545,0.3045902987101033,0.2366523829135607,0.1928081047153158,0.1623776608968675,0.1401041834530502,0.1231319632579329,0.1097872825783083,0.09902859647173193,0.09017567550106468,0.08276628650136917,0.0764757610162485,0.07106958053885211};
    double c=0.918938533204672;
    double y,a,b,q,s,h;
    int j,z;
    j=floor(Min((abs(p)+0.5), 14));
    z=j;
    h=abs(p)-z;
    a=A[j];
    b=z*a-1;
    q=1;
    s=a+h*b;
    for (int i=2; i<24-j; i=i+2)
    {
        a=(a+z*b)/i;
        b=(b+z*a)/(i+1);
        q=q*pow(h,2);
        s=s+q*(a+h*b);
    }
    y=s*exp(-0.5*pow(p,2)-c);
    if (p>0) {
        y=1-y;
    }
    return y;
}
double a[]={2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637};
double b[]={-8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833};
double d[]={0.3374754822726147, 0.9761690190917186, 0.1607979714918209, 0.0276438810333863, 0.0038405729373609, 0.0003951896511919, 0.0000321767881768, 0.0000002888167364, 0.0000003960315187};

double Beasley(double l)
{
    double y = l-0.5;
    double r;
    double x;
    if (abs(y)<0.42)
    {
        r=y*y;
        x=y*(((a[3]*r+a[2])*r+a[1])*r+a[0])/((((b[3]*r+b[2])*r+b[1])*r+b[0])*r+1);
    }
    else
    {
        r=l;
        if (y>0) {
            r=1-l;
        }
        r=log(-log(r));
        x=d[0]+r*(d[1]+r*(d[2]+r*(d[3]+r*(d[4]+r*(d[5]+r*(d[6]+r*(d[7]+r*d[8])))))));
        
        if (y<0) {
            x=-x;
        }
    }
    return x;
}

double get_gaussian_det(double s)
{
    double l=s;
    double b0=Beasley(l);
    double cdf = get_cdf(b0);
    double b1 = b0 - (cdf-l)/(exp(-pow(b0, 2)/2)/sqrt(2*3.141592654));
    return b1;
}
double QMC_Sobol(int num, int batch)
{
    double z, z2;
    z=0;
    z2=0;
    srand((int)time(0));
    for (int j = 0 ; j<batch; j++)
    {
        double x;
        x =0;
        double U[m];
        for (int i =0 ; i<m; i++)
        {
            srand(((int)time(0)*10000*(j+1))*batch);
            U[i]=UNIRAN();
        }
        for (unsigned long long i = 0; i < num; ++i)
        {
    
            double P[m]; //Price
            const double s = sobol::sample(i+1, 0);
            double inter = s+U[0] - floor(s+U[0]);
            const double RandV = (inter);
            P[0]=S0*exp((r-q-pow(sigma,2)*0.5)*dT*(1)+sigma*sqrt(dT*(1))*RandV);
            for (unsigned d = 1; d < m; ++d)
            {
                const double s = sobol::sample(i+1, d);
                double inter = s+U[d] - floor(s+U[d]);
                const double RandV = get_gaussian_det(inter);
                P[d]=P[d-1]*exp((r-q-pow(sigma,2)*0.5)*dT*(1)+sigma*sqrt(dT*(1))*RandV);
            }
            double tem = 0;
            for (int j = 0; j<m; j++) {
                tem += P[j];
            }
            tem = tem/m;
            tem = Max(tem-K, 0);
            x += tem;
        }
        x = x/num;
        double x2= pow(x, 2);
        z += x;
        z2 +=x2;
    }
    expc = z/batch;
    se = sqrt(((z2/batch)-pow(expc, 2))/(batch-1));
    return expc;
}

int main(int argc, const char * argv[])
{
    int trials =100000;
    sscanf (argv[1], "%d", &m);
    inv_m = 1/double(m);
    dT = T/double(m);
    sscanf (argv[1], "%d", &m);
    int batchs =10;
    clock_t t1,t2 = 0;
    t1=clock();
    double Asian_Call;
    Asian_Call= QMC_Sobol(trials, batchs);
    double discount=exp(-r*T);
    cout << "Total trials are: " <<trials<<" * "<<batchs<< " = "<<batchs*trials<<endl;
    cout << "-------------------------------" <<endl;
    cout << "The price is: " << discount*expc<<endl;
    cout << "-------------------------------" <<endl;
    cout << "The Standard Error is: "<< discount*se<<endl;
    cout << "-------------------------------" <<endl;
    t2=clock();
    float seconds = ((float)t2-float(t1))/CLOCKS_PER_SEC;
    cout<<"This program took me: "<<seconds<<" seconds to run"<<endl;

}
