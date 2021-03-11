#include "subIsing.h"

typedef vector< vector<int> > Matrix;
typedef vector<int> Row;

int main(int argc, char** argv)
{

    if (argc != 4)
    {
        cout << "Input parameter: <n> <nsweep> <# of runs>" << endl;
        cout << "See README for details." << endl;
        return 0;
    }

    int n = atoi(argv[1]);
    int m = n;

    int nsweep = atoi(argv[2]);

    double loop_density = 0.2;
    int nloops = (int)floor(n*loop_density);

    int runs = atoi(argv[3]);

    vector<vector<double>> time;
    vector<vector<double>> engdev;

    cout << endl;
    cout << "This example tests the parallel efficiency of the SA solver." << endl;
    cout << "The test is done on a " << n << "-by-" << m << " lattice with " << nsweep << " SA sweeps." << endl << endl;

    for (int i=1; i<=24; i++)
    {

        int nthreads = i;
        omp_set_num_threads(nthreads);

        Ising lattice(n,m,nthreads,0);
        lattice.rand_spins();

        vector<double> time_row;
        vector<double> engdev_row;

        cout << "Solving " << runs << " instances with " << nthreads << " threads." << endl;

        for (int seed=0; seed<runs; seed++)

        {

            lattice.add_loops(nloops, seed);
            int eplant = lattice.energy_sol;
            int esol = 0;
            double t = 0;

            if (nthreads == 1)
            {
                t = lattice.sweeps(nsweep,seed,true);
                lattice.update_energy();
                esol = lattice.energy;
            }
            else
            {
                subIsing sublattice(lattice);
                t = sublattice.sweeps_parallel(nsweep,seed);
                sublattice.update_energy();
                esol = sublattice.energy;
            }

            double ed = (eplant - esol)/(double)eplant;
            lattice.reset(seed+1);
            time_row.push_back(t);
            engdev_row.push_back(ed);
        }

        double t_tot = accumulate(time_row.begin(), time_row.end(), 0.0);
        double ed_avg = accumulate(engdev_row.begin(), engdev_row.end(), 0.0)/engdev_row.size();

        cout << "Total time used: " << t_tot << "." << endl;
        cout << "Average fractional energy deviation: " << ed_avg << "." << endl << endl;

        time.push_back(time_row);
        engdev.push_back(engdev_row);

    }

    ofstream file_t;
    ofstream file_eng;
    file_t.open("test_sol_t.csv");
    file_eng.open("test_sol_edev.csv");

    for (int i=0; i<24; i++)
    {
        for (int j=0; j<runs; j++)
        {
            file_t << time[i][j] << ",";
            file_eng << engdev[i][j] << ",";
        }
        file_t << "\n";
        file_eng << "\n";
    }

    return 0;

}

