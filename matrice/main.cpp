#include <iostream>
#include "matrice.h"
#include <stdlib.h>
#include <math.h>

using namespace std;
void GC(matrice& x, matrice& A, matrice& b, matrice& x0, int iterations_max);

int main()
{
    /*
    cout << "Hello world!" << endl;
    int nbl=5, nbc=5;
    matrice A(nbl,nbc);
    matrice B = A; // ou matrice B(A);
    matrice C(nbl,nbc);
    C = A; //operator=
    int compteur = 0;
    for(int i=0;i<nbl;i++)
        for(int j=0;j<nbc;j++)
            A[i][j] = ++compteur; //<->{compteur+=1; A[i][j] = compteur;}
            //A[i][j] = compteur++ <-> {A[i][j] = compteur; compteur+=1;}


    A.afficher();
    cout << endl;
    A.transpose().afficher();


    double *v1 = new double[nbc];
    for(int i=0;i<nbc;i++)
        v1[i] = rand()/((double) RAND_MAX);
    double *v2 = A.multipli(v1);

    A.multipli(v2, v1);
     // v2 = A*v1
    matrice V1(nbc,1,v1);
    cout << endl;
    cout << "Vecteur v1" << endl;
    V1.afficher();

    cout << endl;
    matrice V2(nbl, 1, v2);
    cout << "Vecteur v2" << endl;
    V2.afficher();



    cout<<endl;

    matrice sousmat = A.sousmatrice(1, 2, 3, 3);
    //sousmat.afficher();

    cout<<endl;

    //A.diag().afficher();
    matrice E(nbl,nbc);
    for(int i=0;i<nbl;i++)
        for(int j=0;j<nbc;j++)
            E[i][j]= rand()/((double) RAND_MAX);

    E.afficher();
    cout<<endl;

    E.exp().afficher();
    cout<<endl;
    */

    //A, --> exp(A) mathématiquement
    // comparer exp(A) avec A.exp()

    /*
    matrice P(2,2);
    P[0][0] = 1;
    P[0][1] = 2;
    P[1][0] = 3;
    P[1][1] = 4;
    matrice Pm(2,2);
    Pm[0][0] = P[1][1];
    Pm[1][1] = P[0][0];
    Pm[1][0] = -P[1][0];
    Pm[0][1] = -P[0][1];
    P.afficher();
    cout<<endl;
    Pm.afficher();
    cout<<endl;
    Pm /= (P[0][0]*P[1][1] - P[0][1]*P[1][0]);
    matrice H = Pm*P;
    H.afficher();
    cout<<endl;

    matrice D(2,2);
    D = D.diag();
    D[0][0]=3;
    D[1][1]=1;
    D.afficher();
    cout<<endl;

    matrice expD = D;
    expD[0][0] = exp(D[0][0]);
    expD[1][1] = exp(D[1][1]);

    matrice A = Pm*D*P;
    A.afficher();
    cout<<endl;

    matrice expA = Pm*expD*P;
    matrice expA2 = A.exp();
    expA.afficher();
    cout<<endl;

    expA2.afficher();
    cout<<endl;
    */

// Tester GC
    int dim = 10;
    matrice B(dim,dim);
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            B[i][j] = 2*rand()/((double) RAND_MAX);

    matrice A = B.transpose()*B;
    matrice x(dim,1);
    for(int i=0;i<dim;i++)
        x[i][0] = rand()/((double) RAND_MAX);
    matrice b = A*x;

    matrice x0(dim,1);
    for(int i=0;i<dim;i++)
        x0[i][0] = 1;

    int iterations_max = 12;

    matrice x1(dim,1);

    GC(x1, A, b, x0, iterations_max);

    cout<<"Solotion par GC :"<<endl;
    x1.afficher();

    cout<<"Vraie solotion :"<<endl;
    x.afficher();

    return 0;
}
