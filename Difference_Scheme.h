#include<math.h>
using namespace std;
void first_windward(double *f,double &fm,int pm) {fm = f[0];}
void VanLeerLim(double *ft,double &ftm,int pm){
    //pm = 1左值，pm=-1右值
    int j=1;
    double dL = ft[j]-ft[j-1*pm],dR=ft[j+1*pm]-ft[j];
    double eps = 1e-12;
    double phi = (abs(dL*dR) + dL*dR+eps)/(abs(dL*dR)+dR*dR+eps);
    ftm = ft[j]+0.5*phi*dR;
}
double minmod(double a,double b){
    if(!((a>0.0)^(b>0.0))) return (abs(a)>abs(b)?b:a);
    else return 0.0;
}
void NND(double *ft,double &ftm,int pm){
    //pm = 1左值，pm=-1右值
    int j=1;
    ftm = ft[j] + 0.5*minmod(ft[j]-ft[j-1*pm],ft[j+1*pm]-ft[j]);
}