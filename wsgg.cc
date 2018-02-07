#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////

double cc[]={ 7.412956e-001, -5.244441e-001,  5.822860e-001, -2.096994e-001,  2.420312e-002,
             -9.412652e-001,  2.799577e-001, -7.672319e-001,  3.204027e-001, -3.910174e-002,
              8.531866e-001,  8.230754e-002,  5.289430e-001, -2.468463e-001,  3.109396e-002,
             -3.342806e-001,  1.474987e-001, -4.160689e-001,  1.697627e-001, -2.040660e-002,
              4.314362e-002, -6.886217e-002,  1.109773e-001, -4.208608e-002,  4.918817e-003,
              1.552073e-001, -4.862117e-001,  3.668088e-001, -1.055508e-001,  1.058568e-002,
              6.755648e-001,  1.409271e+000, -1.383449e+000,  4.575210e-001, -5.019760e-002,
             -1.125394e+000, -5.913199e-001,  9.085441e-001, -3.334201e-001,  3.842361e-002,
              6.040543e-001, -5.533854e-002, -1.733014e-001,  7.916083e-002, -9.893357e-003,
             -1.105453e-001,  4.646634e-002, -1.612982e-003, -3.539835e-003,  6.121277e-004,
              2.550242e-001,  3.805403e-001, -4.249709e-001,  1.429446e-001, -1.574075e-002,
             -6.065428e-001,  3.494024e-001,  1.853509e-001, -1.013694e-001,  1.302441e-002,
              8.123855e-001, -1.102009e+000,  4.046178e-001, -8.118223e-002,  6.298101e-003,
             -4.532290e-001,  6.784475e-001, -3.432603e-001,  8.830883e-002, -8.415221e-003,
              8.693093e-002, -1.306996e-001,  7.414464e-002, -2.029294e-002,  2.010969e-003,
             -3.451994e-002,  2.656726e-001, -1.225365e-001,  3.001508e-002, -2.820525e-003,
              4.112046e-001, -5.728350e-001,  2.924490e-001, -7.980766e-002,  7.996603e-003,
             -5.055995e-001,  4.579559e-001, -2.616436e-001,  7.648413e-002, -7.908356e-003,
              2.317509e-001, -1.656759e-001,  1.052608e-001, -3.219347e-002,  3.386965e-003,
             -3.754908e-002,  2.295193e-002, -1.600472e-002,  5.046318e-003, -5.364326e-004};  

double dd[]={ 3.404288e-002,  6.523048e-002, -4.636852e-002,  1.386835e-002, -1.444993e-003,
              3.509457e-001,  7.465138e-001, -5.293090e-001,  1.594423e-001, -1.663261e-002,
              4.570740e+000,  2.168067e+000, -1.498901e+000,  4.917165e-001, -5.429990e-002,
              1.098169e+002, -5.092359e+001,  2.343236e+001, -5.163892e+000,  4.393889e-001};

vector<vector<vector<double> > > c(4, vector<vector<double> >(5, vector<double>(5, 0.0)));
vector<vector<double> >          d(5, vector<double>(5, 0.0));

void set_c_d(){

    int ni = 4;
    int nj = 5;
    int nk = 5;
    
    for(int i=0; i<ni; i++)
        for(int j=0; j<nj; j++)
            for(int k=0; k<nk; k++)
                c[i][j][k] = cc[ i*nj*nk + j*(nk) + k ];

    for(int i=0; i<ni; i++)
        for(int k=0; k<nk; k++)
            d[i][k] = dd[ i*nk + k ];

    return;
}


//////////////////////////////////////////////////////////////////////////////////////////////
/** compute gray gas absorption coefficients and weighting factors
  * number of gray gases = 4 plus 1 for the clear gas.
  * @param T    \input temperature (K)
  * @param P    \input pressure (Pa)
  * @param XCO2 \input CO2 mole fraction
  * @param XH2O \input H2O mole fraction
  * @param K    \output array of absorption coefficients K[0] = 0 (clear gas), units (1/m).
  * @param a    \output array of weights: sum = 1.
**/

void get_k_a(const double T, 
             const double P, 
                   double XCO2, 
             const double XH2O, 
             vector<double> &K, 
             vector<double> &a){

    //------------------------
    static bool c_d_are_set = false;
    if(!c_d_are_set){
        set_c_d();
        c_d_are_set = true;
    }
    //------------------------

    if(abs(XCO2) < 1E-6) XCO2 = 1E-6;
    double Mr = XH2O/XCO2;
    //if(Mr<0.01) Mr = 0.01;
    //if(Mr>4)    Mr = 4;

    //if(T<500)  T = 500;
    //if(T>2400) T = 2400;
    double Tr = T/1200;

    int nGG = 5;           // 4 gray gases and 1 clear gas
    int ni  = 4;
    int nj  = 5;
    int nk  = 5;

    K = vector<double>(nGG,0.0);
    a = vector<double>(nGG,0.0);
    vector<vector<double> > b(ni, vector<double>(nj,0.0));

    for(int i=0; i<ni; i++)
        for(int j=0; j<nj; j++)
            for(int k=0; k<nk; k++)
                //b[i][j] += cc[i*nj*nk+j*nk+k]*pow(Mr,k);
                b[i][j] += c[i][j][k]*pow(Mr,k);

    a[0] = 1.0;
    for(int i=1; i<nGG; i++){
        for(int j=0; j<nj; j++)
            a[i] += b[i-1][j]*pow(Tr,j);
        a[0] -= a[i];
    }


    for(int i=1; i<nGG; i++){
        for(int k=0; k<nk; k++)
            //K[i] += dd[(i-1)*nk+k]*pow(Mr,k);
            K[i] += d[i-1][k]*pow(Mr,k);
        K[i] *= (P/101325)*(XH2O+XCO2);
    }

    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////

//int main() {
//
//    vector<double> K;
//    vector<double> a;
//    get_k_a(1000,101325,0.1,0.1, K, a);
//    
//    cout << endl << "K =";
//    for(int i=0; i<K.size(); i++)
//        cout << endl << K[i];
//    
//    cout << endl << "a =";
//    for(int i=0; i<a.size(); i++)
//        cout << endl << a[i];
//
//    return 0;
//
//}
