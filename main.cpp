#include "subIsing.h"

typedef vector< vector<int> > Matrix;
typedef vector<int> Row;

int main(int argc, char** argv)
{

    if (argc != 7)
    {
        cout << "Input parameter: <n> <m> <loop_density> <SA sweeps> <# of threads> <seed>" << endl;
        cout << "See README for details." << endl;
        return 0;
    }

    int n = atoi(argv[1]);
    int m = atoi(argv[2]);

    double loop_density = atof(argv[3]);
    int nsweep = atof(argv[4]);

    int nthreads = atoi(argv[5]);
    int seed = atoi(argv[6]);

    omp_set_num_threads(nthreads);
    int nloops = (int)floor(n*loop_density);

    if (nthreads == 1)
    {

        cout << endl;
        cout << "This is the sequential version of the code for 2D Ising model instance generation and solution." << endl << endl;

        cout << "Initializing " << n << "-by-" << m << " Ising lattice..." << endl << endl;
        Ising lattice(n, m, nthreads, 0);
        lattice.rand_spins();

        cout << "Generating instance with " << nloops << " frustrated loops..." << endl;
        double time_gen = lattice.add_loops(nloops,seed);
        cout << "Time used for generation: " << time_gen << " seconds." << endl;
        int eplant = lattice.energy_sol;
        cout << "The planted ground state energy is " << -eplant << "." << endl << endl;

        cout << "Solving the instance with SA with " << nsweep << " sweeps..." << endl;
        double time_sol = lattice.sweeps(nsweep, 0, true);
        cout << "Time used for solution: " << time_sol << " seconds." << endl;

        lattice.update_energy();
        int esol = lattice.energy;
        cout << "The obtained minimum energy is " << -esol << "." << endl;
        cout << "Fractional energy deviation: " << setprecision(2) << 100.0*(eplant - esol)/eplant << "%.";

        cout << endl;

    }
    else
    {

        cout << endl;
        cout << "This is the parallelized version of the code for 2D Ising model instance generation and solution." << endl;
        cout << "Number of threads: " << nthreads << endl << endl;

        cout << "Initializing " << n << "-by-" << m << " Ising lattice..." << endl << endl;
        Ising lattice(n, m, nthreads, 0);
        lattice.rand_spins();

        cout << "Generating instance with " << nloops << " frustrated loops..." << endl;
        double time_gen = lattice.add_loops_parallel(nloops,seed);
        cout << "Time used for generation: " << time_gen << " seconds." << endl;
        int eplant = lattice.energy_sol;
        cout << "The planted ground state energy is " << -eplant << "." << endl << endl;

        Row part = getpart(nthreads);
        cout << "Partitioning the lattice with a " << part[0] << "-by-" << part[1] << " scheme..." << endl << endl;
        subIsing sublattice(lattice);

        cout << "Solving the instance with SA with " << nsweep << " sweeps..." << endl;
        double time_sol = sublattice.sweeps_parallel(nsweep, 0);
        cout << "Time used for solution: " << time_sol << " seconds." << endl;

        sublattice.update_energy();
        int esol = sublattice.energy;
        cout << "The obtained minimum energy is " << -esol << "." << endl;
        cout << "Fractional energy deviation: " << setprecision(2) << 100.0*(eplant - esol)/eplant << "%.";

        cout << endl;

    }

    return 0;

}

