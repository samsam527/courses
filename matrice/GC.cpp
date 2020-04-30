
#include "matrice.h"
#include <iostream>


void GC(matrice& x, matrice& A, matrice& b, matrice& x0, int iterations_max)
{
    int dim = A.getnbl();
    matrice x_k(dim,1), x_km1(dim,1), r_k(dim,1), r_km1(dim,1), r_km2(dim,1),
    p_k(dim,1), p_km1(dim,1);

    //etape :  k = 0;
    x_k = x0;
    r_k = b - A * x_k;

    double seuil = 1e-15;
    for(int k=1; k<=iterations_max && r_k.norm()>seuil; k++)
    {
        //std::cout<<"residu "<< r_k.norm()<<std::endl;
        r_km2 = r_km1;
        r_km1 = r_k;
        x_km1 = x_k;
        p_km1 = p_k;
        if(k==1)
        {
            p_k = r_km1;
        }
        else
        {
            double a = r_km1.norm();
            a*=a;
            double b = r_km2.norm();
            b*=b;
            double c_k = a/b;
            p_k = r_km1 + p_km1*c_k;
        }
        double a_k = r_km1.norm();
        a_k *= a_k;
        matrice p_k_t = p_k.transpose();
        matrice B = p_k_t * A * p_k;
        a_k /= B[0][0];
        x_k = x_km1 + p_k * a_k;
        r_k = r_km1 - A * p_k * a_k;
    }

    x = x_k;
}
