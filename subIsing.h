#ifndef SUBISING_H_INCLUDED
#define SUBISING_H_INCLUDED

#include "Ising.h"
using namespace std;

typedef vector< vector< vector<int> > > Tensor;
typedef vector< vector<int> > Matrix;
typedef vector<int> Row;

typedef vector<Ising> Rowlattice;
typedef vector< vector<Ising> > Matlattice;

class subIsing
{

    public:

    subIsing(Ising &lattice);

    int energy;

    void update_energy();

    double sweeps(int nsweep, int seed1);
    double sweeps_parallel(int nsweep, int seed1);

    void print_spins(int i, int j);
    void print_fields(int i, int j);
    void print_bonds(int i, int j);
    void print_inout(int i, int j);

    void combine();
    void print_combine_spins();
    void print_combine_fields();

    private:

    int n;
    int m;
    Matrix nmat;
    Matrix mmat;

    int k1;
    int k2;
    bool corner;

    Matrix spins;
    Matrix fields;
    Tensor bonds;
    Matlattice matlattice;

    void separate_bonds();
    void contact(Ising &mid, Ising &down, Ising &up);
    void contact(Ising &mid, Ising &left, Ising &right, Ising &down, Ising &up);
    void auto_contact();

    void exchange_fields(Ising &mid, Ising &down, Ising &up);
    void exchange_fields(Ising &mid, Ising &left, Ising &right, Ising &down, Ising &up);
    void auto_exchange();

};

#endif // SUBISING_H_INCLUDED
