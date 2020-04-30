#include "matrice.h"
#include <iostream>
#include <math.h>
using namespace std;

matrice::matrice()
{
    nbl = nbc = 1;
    T = new double * [nbl];
    T[0] = new double[nbc];
}

matrice::matrice(int nbl, int nbc)
{
    this->nbl = nbl;
    this->nbc = nbc;

    T = new double * [nbl];
    for (int i = 0; i< nbl; i++)
        T[i] = new double [nbc];
}

matrice::matrice(const matrice &A)
{
    nbl = A.nbl;
    nbc = A.nbc;

    T = new double * [nbl];
    for (int i = 0; i< nbl; i++)
        T[i] = new double [nbc];

     for (int i = 0; i< nbl; i++)
        for (int j = 0; j< nbc; j++)
            T[i][j]=A.T[i][j];

}

matrice::matrice(int nbl, int nbc, double *tab)
{
    this->nbl = nbl;
    this->nbc = nbc;

    T = new double * [nbl];
    for (int i = 0; i< nbl; i++)
        T[i] = new double [nbc];

//    for (int i = 0; i< nbl; i++)
//        for (int j = 0; j< nbc; j++)
//            T[i][j] = tab[i*nbc+j];

    int compteur=0;
    for (int i = 0; i< nbl; i++)
        for (int j = 0; j< nbc; j++)
            T[i][j] = tab[compteur++];
}

matrice::~matrice()
{
    for (int i = 0; i< nbl; i++)
        delete[] T[i];
    delete[] T;
}

matrice &matrice::operator=(const matrice &A)
{
    for(int i=0;i<nbl;i++)
        for(int j=0;j<nbc;j++)
            T[i][j]=A.T[i][j];

    return *this;
}

matrice matrice::operator+(const matrice &A) const
{
    matrice resultat(nbl, nbc);
    for (int i(0); i< nbl;i++)
    {
        for (int j(0); j<nbc;j++)
        {
            resultat[i][j] = T[i][j] + A[i][j];
        }
    }
    return resultat;
}

matrice matrice::operator-(const matrice &A) const
{
    matrice resultat(nbl, nbc);
    for (int i(0); i< nbl;i++)
    {
        for (int j(0); j<nbc;j++)
        {
            resultat[i][j] = T[i][j] - A[i][j];
        }
    }
    return resultat;
}

double* matrice::operator[](int i) const
{
    return T[i];
}

void matrice ::  afficher() const
{
    for (int i=0; i<nbl;i++)
    {
        for (int j(0); j<nbc;j++)
            cout << T[i][j] << " " ;
        cout << endl;
    }
}

matrice& matrice::operator+=(const matrice &A)
{
    for (int i(0); i< nbl;i++)
        for (int j(0); j<nbc;j++)
            T[i][j] += A[i][j];

    return *this ;
}

matrice& matrice::operator-=(const matrice &A)
{
    for (int i(0); i< nbl;i++)
        for (int j(0); j<nbc;j++)
            T[i][j] -= A[i][j];

    return *this ;
}

matrice matrice::operator*(double x) const
{
    matrice A(*this);
     for (int i(0); i< nbl;i++)
        for (int j(0); j<nbc;j++)
            A[i][j]*=x;
    return A;
}

matrice matrice::operator/(double x) const
{
    matrice A = (*this)*(1.0/x);
    return A;
}

matrice matrice::operator*(const matrice &A) const
{
    matrice B(nbl,A.nbc);
     for (int i(0); i<nbl;i++)
        for (int j(0); j<A.nbc;j++)
        {
            B[i][j]=0;
            for(int k=0;k<nbc;k++)
                B[i][j]+=T[i][k]*A[k][j];
        }
    return B;
}

matrice& matrice::operator*=(const matrice &A)
{
    *this = (*this)*A;
    return *this;
}

matrice& matrice::operator*=(double x)
{
    for (int i(0); i< nbl;i++)
        for (int j(0); j<nbc;j++)
            T[i][j]*=x;
    return *this;
}

matrice& matrice::operator/=(double x)
{
    *this *= 1.0/x;
    return *this;
}

matrice matrice::sousmatrice(int init_ligne, int init_colonne, int dim_nbl, int dim_nbc) const
{
    matrice mat(dim_nbl, dim_nbc);

    for(int i=0;i<dim_nbl;i++)
        for(int j=0;j<dim_nbc;j++)
            mat[i][j] = T[i+init_ligne][j+init_colonne];

    return mat;
}

matrice matrice::diag()
{
   matrice A(*this);

   for(int i=0;i<nbl;i++)
        for(int j=0;j<nbc;j++)
            if (i!=j)
                A[i][j]=0;

   return A;

}

void matrice::multipli(double *v_out, double *v_in)
{
    for(int i(0); i<nbl;i++)
    {
       double somme(0);
       for(int j(0);j<nbc;j++)
            somme += T[i][j]*v_in[j];
       v_out[i] = somme;
    }
}

double* matrice::multipli(double *v_in)
{
    double *v = new double[nbl];
    multipli(v, v_in);
    return v;
}

matrice matrice::transpose() const
{
    matrice P(nbc,nbl);
    for (int i(0); i<nbc;i++)
    {
        for (int j(0); j< nbl; j++)
            P[i][j] = (*this)[j][i];
    }
    return P;
}

matrice matrice::exp(int n) const
{
    matrice T(nbl,nbc);
    for(int i=0;i<nbl;i++)
        for(int j=0;j<nbc;j++)
            T[i][j] = (i==j ? 1 : 0);

    matrice S = T;
    for(int i=1;i<n;i++)
    {
        T *= (*this)/i;
        S += T;
    }
    return S;
}

double matrice::norm() const
{
    double s=0;
    for(int i=0;i<nbl;i++)
        for(int j=0;j<nbc;j++)
            s += T[i][j]*T[i][j];
    return sqrt(s);
}
