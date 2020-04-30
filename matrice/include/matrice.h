#ifndef MATRICE_H
#define MATRICE_H


class matrice
{
    public:
        matrice();
        matrice(int nbl, int nbc);
        matrice(const matrice &A);
        matrice(int nbl, int nbc, double *tab);
        virtual ~matrice();
        matrice &operator=(const matrice &A);
        matrice operator+(const matrice &A) const;
        matrice operator*(double x) const;
        matrice operator/(double x) const;
        matrice operator*(const matrice &A) const;
        matrice operator-(const matrice &A) const;
        double* operator[](int i) const;
        matrice& operator+=(const matrice &A);
        matrice& operator-=(const matrice &A);
        matrice& operator*=(const matrice &A);
        matrice& operator*=(double x);
        matrice& operator/=(double x);
        void afficher() const;
        matrice sousmatrice(int init_ligne, int init_colonne, int dim_nbl, int dim_nbc) const;
        matrice diag();
        double* multipli(double *v_in);
        void multipli(double *v_out, double *v_in);
        matrice transpose() const;
        matrice exp(int n=100) const;
        int getnbl() {return nbl;}
        int getnbc() {return nbc;}
        double norm() const;


    protected:

    private:
        double **T;
        int nbl, nbc;
};

#endif // MATRICE_H
