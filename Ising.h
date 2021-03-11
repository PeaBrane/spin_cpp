#ifndef ISING_H_INCLUDED
#define ISING_H_INCLUDED

#include "Cycle.h"

using namespace std;

typedef vector< vector< vector<int> > > Tensor;
typedef vector< vector<int> > Matrix;
typedef vector<int> Row;

class Ising
{

    public:

    Ising(int n0, int m0, int nthreads0, int seed0);
    Ising(Matrix &spins0, Tensor &bonds0, Matrix &fields0);

    Matrix spins;
    Matrix fields;
    Tensor bonds;
    Matrix bonds_right;
    Matrix bonds_up;

    Row b_left;
    Row b_right;
    Row b_down;
    Row b_up;

    Row in_left;
    Row in_right;
    Row in_down;
    Row in_up;
    Row out_left;
    Row out_right;
    Row out_down;
    Row out_up;

    void eye_spins();
    void rand_spins();

    void ferro_bonds();
    double add_loops(int nloops, int seed0);
    double add_loops_parallel(int nloops, int seed0);

    void separate_bonds();
    void update_fields();
    void update_energy();

    void sweep_bound(double beta, int seed0, bool corner);
    double sweeps(int nsweep, int seed0, bool r);

    void reset(int seed0);

    void print_spins();
    void print_fields();
    void print_bonds();

    int n;
    int m;
    int nthreads;
    int energy;
    int energy_sol = 0;

    private:

    int seed;

    void add_bond(Row point1, Row point2, int bond);
    int add_loop(Matrix cycle);

    void flip(int i,int j,double beta, mt19937& rng);
    void flip_bound(int i,int j,double beta, mt19937& rng, bool corner);
    void sweep_row(int i,double beta, mt19937& rng);
    void sweep(double beta, mt19937& rng);
    void sweep_rand(double beta, mt19937& rng);

};

#endif // 2D_ISING_H_INCLUDED
