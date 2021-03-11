#include "subIsing.h"

typedef vector< vector<int> > Matrix;
typedef vector<int> Row;

int main(int argc, char** argv)
{

    if (argc != 3)
    {
        cout << "Input parameter: <n> <loop_density>" << endl;
        cout << "See README for details." << endl;
        return 0;
    }

    int n = atoi(argv[1]);
    int m = n;

    double loop_density = atof(argv[2]);

    int nloops = (int)floor(n*loop_density);

    vector<double> time;

    cout << endl;
    cout << "This example tests the parallel efficiency of the loop algorithm." << endl;
    cout << "The test is done on a " << n << "-by-" << m << " lattice with " << nloops << " loops." << endl << endl;

    for (int i=1; i<=24; i++)
    {

        int nthreads = i;
        omp_set_num_threads(nthreads);

        Ising lattice(n,m,nthreads,0);
        lattice.rand_spins();

        cout << "Using " << nthreads << " threads..." << endl;
        double t = 0;
        if (nthreads == 1) {t = lattice.add_loops(nloops,0);}
        else {t = lattice.add_loops_parallel(nloops,0);}
        time.push_back(t);
        cout << "Time used: " << t << " seconds." << endl << endl;

    }

    ofstream myfile;
    myfile.open("test_gen.csv");

    for (int i=1; i<=24; i++)
    {
        myfile << time[i-1] << ",";
    }

    return 0;

}

