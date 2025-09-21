
#include <iostream>
#include <stdio.h>
#include "math.h"

float A=0.0, B=1.0, h=0.5;
float f(float x){
    float f0=-1.0/12.0, f1=-1.0/6.0, f2=1.0/6.0, f3=0.5, f4=0, f5=0; //koef. f(x)
    return f0*pow(x,5)+f1*pow(x,4)+f2*pow(x,3)+f3*pow(x,2)+f4*x+f5;
}
//Runge-Kutta
int k;
const int n2=2; //   =(b-a)/h
float x2[n2], y2[n2]; //razbienie otr.
float K1, K2, K3, K4;
float PravCh(float x){
    return sin(x)*cos(x);
}
//interpol. polinom
float P2(float x, float* x2, float* y2, int n2){
    float _P2 = 0; int k,i;
    float Proizv = 1;
    for(i=0;i<n2;i++){
        for(k=0;k<n2;k++) {
            if(k!=i)Proizv *= (x-x2[k])/(x2[i]-x2[k]);
        }
        _P2 += y2[i]*Proizv;
    }
    return _P2;
}
//maximum
float eps=0.0001; //pogreshn
float x, max, aa, bb;
float fi, X1, X2, Y1, Y2;
float F(float x, float* x2, float* y2, int n2){
    return fabs(P2(x,x2,y2,n2)-f(x));
}
//
float xRez[20], FRez[20];//для записи результата

int main(int argc, char *argv[])
{
    
    FILE* file1;  int i=0;
    //Runge-Kutta
    for(k=0;k<n2;k++) x2[k]=A+k*h;
    y2[0]=A;
    for(k=0;k<n2;k++){
        K1=h*PravCh(x2[k]); K2=h*PravCh(x2[k]+0.5*h);
        K3=h*PravCh(x2[k]+0.5*h); K4=h*PravCh(x2[k]+h);
        y2[k+1]=y2[k]+(1.0/6.0)*(K1+2*K2+2*K3+K4);
    }
    //
    max=0.0, x=0.0, aa=A, bb=B, fi = 1.618;
    X1=aa, X2=bb;
    while(X2-X1 > eps) {
        X1 = bb-(bb-aa)/fi; X2 = aa+(bb-aa)/fi;
        Y1 = F(X1,x2,y2,n2); Y2 = F(X2,x2,y2,n2);
        if(Y1<=Y2) aa = X1;
        else bb = X2;
    }
    x = (aa+bb)/2.0;
    max = F(x,x2,y2,n2);
    std::cout << "\t x="<<x<<" max="<<max<<"\n";
    //
    for(i=0;i<=20;i++){
        xRez[i] = A+i*(B-A)/20.0; FRez[i] = F(xRez[i],x2,y2,n2);
    }
    fopen_s(&file1, "result_x.txt", "w");//Открывается файл для записи результатов
    for(i=0;i<=20;i++)fprintf(file1,"%f\n",xRez[i]); fclose(file1);
    fopen_s(&file1, "result_F.txt", "w");//Открывается файл для записи результатов
    for(i=0;i<=20;i++)fprintf(file1,"%f\n",FRez[i]); fclose(file1);
    fopen_s(&file1, "result_P2.txt", "w");//Открывается файл для записи результатов
    for(i=0;i<=20;i++)fprintf(file1,"%f\n",P2(xRez[i],x2,y2,n2));  fclose(file1);
    fopen_s(&file1, "result_ff.txt", "w");//Открывается файл для записи результатов
    for(i=0;i<=20;i++)fprintf(file1,"%f\n",f(xRez[i]));  fclose(file1);
    std::cout << "files ready"<<"\n";

    
}
