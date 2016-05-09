/* Recurrence relation, Legendre polynomial

 Numerical recipes in C, 1992
 Press, Flannery, Teukolsky, Vetterling

 Chapter 6., Special Functions, 
 Section 6.8 Spherical Harmonics, p.252-254

 The Legendre polynomial recurrence relation is written as

 (n+1)*P_{n+1}(x)=(2n+1)*x*P_n(x)-n*P_{n-1}(x)

 float plgndr(int l, int m, float x)
 Computes the associated Legendre polynomial P^m_l(x). Here m and l
 are integers satisfying 0<=m<=l, while x lies in the range
 -1<=x<=1.

 P(0,x) = 1
 P(1,x) = x
 P(n,x) = (2*n-1)/n * x * P(n-1,x) - (n-1)/n * P(n-2,x)
 
*/
/*
 Recall: Floating point Types
 
 float 
 Single precision floating point number typical size: 32 bits
 double
 Double precision floating point number typical size: 64 bits
 long double
 Possibly even bigger floating point number (somewhat obscure)
 
*/
#include <stdio.h>    // standard file input and output header
#include <stdlib.h>   // utility functons such as malloc() and rand()
#include <string.h>
#include <math.h>     // mathematical functions such as sin() and cos()}
#include <stddef.h>   // for error function

/* function declarations */
double plgndr(int l, int m, double x);


/* program begins */
int main()
{
    /* local variable definition */
    int l = 2; /* Asign values to variables here */
    int m = 2;
    double x = 0.5;   // x must be between -1 and +1
    double ret;       // return function
    /* call function */
    ret = plgndr(l,m,x); //declare function
    printf( "The value is : %f\n", ret); //prints value
    return 0;
}

/* functions */
void nrerror(char error_text[]) //Numerical Recipes standard error handler
{
    fprintf(stderr, "Run time error....\n");
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, "Now exiting the system\n");
    exit(1);
}
double plgndr(int l, int m, double x) //originally written as float
/* Computes the associated Legendre polynomial P^m_l(x). Here m and l
   are integers satisfying 0<=m<=l, while x lies in
   the range -1<=x<=1. */
{
    void nrerror(char error_text[]);
    double fact,pll,pmm,pmmp1,somx2;
    int i,ll;
    
    if (m < 0 || m > l || fabs(x) > 1.0)
        nrerror("Bad arguments in routine plgndr");
    pmm=1.0;          //Compute P^m_m
    if (m > 0){
        somx2=sqrt((1.0-x)*(1.0+x));
        fact=1.0;
        for (i=1; i<=m; i++){
            pmm *= -fact*somx2;
            fact += 2.0;
        }
    }
    if (l == m)
        return pmm;
    else{             //Compute P^m_{m+1}
        pmmp1=x*(2*m+1)*pmm;
        if (l == (m+1))
            return pmmp1;
        else {         //Compute P^m_l, l > m+1
            for (ll=m+2; ll<=l; ll++){
                pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
                pmm=pmmp1;
                pmmp1=pll;
            }
            return pll;
        }
    }
}

    
    

