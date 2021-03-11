#include "Ising.h"

Ising::Ising(int n0, int m0, int nthreads0, int seed0)
{

    n = n0;
    m = m0;
    energy = 0;
    nthreads = nthreads0;
    seed = seed0;

    spins = mat_empty(n,m);
    fields = mat_empty(n,m);
    bonds = ten_empty(n,m,2);

}

Ising::Ising(Matrix &spins0, Tensor &bonds0, Matrix &fields0)
{

    spins = spins0; bonds = bonds0; fields = fields0;
    n = spins.size();
    m = spins[0].size();

}

void Ising::eye_spins()
{

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            spins[i][j] = 1;
        }
    }

}

void Ising::rand_spins()
{

    mt19937 rng(seed);
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            spins[i][j] = -1 + 2*uid(rng);
        }
    }

}

void Ising::ferro_bonds()
{

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            int i1 = (i+n-1)%n; int i2 = (i+1)%n;
            int j1 = (j+m-1)%m; int j2 = (j+1)%m;
            bonds[i][j][0] = 1; bonds[i][j][1] = 1;
            fields[i][j] = spins[i1][j] + spins[i2][j] + spins[i][j1] + spins[i][j2];
        }
    }

    energy_sol = 2*n*m;

}

double Ising::add_loops(int nloops, int seed0)
{

    double start_time = omp_get_wtime();

    Cycle cycle(n,m,seed0);
    for (int i=1; i<nloops; i++)
    {
        cycle.walk();
        Matrix c = cycle.close;
        int seed1 = seed0*nloops + i;
        cycle.reset(seed1);
        energy_sol += add_loop(c);
    }

    return omp_get_wtime() - start_time;

}

double Ising::add_loops_parallel(int nloops, int seed0)
{

    double start_time = omp_get_wtime();

    int temp = 0;
    Cycle cycle(n,m,0);

    # pragma omp parallel for firstprivate(cycle) reduction(+:temp)
    for (int i=0; i<nloops; i++)
    {
        int seed1 = seed0*nloops + i;
        cycle.reset(seed1);
        cycle.walk();
        Matrix c = cycle.close;
        temp += add_loop(c);
    }
    energy_sol += temp;

    return omp_get_wtime() - start_time;

}

void Ising::separate_bonds()
{

    bonds_right = getmat(bonds,0);
    bonds_up = getmat(bonds,1);

    b_left = getrow(bonds_right,1);
    b_right = getrow(bonds_right,3);
    b_down = getrow(bonds_up,2);
    b_up = getrow(bonds_up,4);

}

void Ising::update_fields()
{

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            int i1 = (i+n-1)%n; int i2 = (i+1)%n;
            int j1 = (j+m-1)%m; int j2 = (j+1)%m;
            fields[i][j] = spins[i1][j]*bonds[i1][j][0] + spins[i2][j]*bonds[i][j][0]
                         + spins[i][j1]*bonds[i][j1][1] + spins[i][j2]*bonds[i][j][1];
        }
    }

}

void Ising::update_energy()
{

    int temp = 0;
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            temp += spins[i][j]*fields[i][j];
        }
    }
    energy = temp/2;

}

void Ising::sweep_bound(double beta, int seed0, bool corner)
{

    mt19937 rng(seed0);

    Matrix mat = randmat(n,m);

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            int x = mat[i][j];
            int ii = x/m; int jj = x%m;
            flip_bound(ii,jj,beta,rng,corner);
        }
    }

    in_left = getrow(spins, 1);
    in_right = getrow(spins, 3);
    in_down = getrow(spins, 2);
    in_up = getrow(spins, 4);

}

double Ising::sweeps(int nsweep, int seed0, bool r)
{

   double start_time = omp_get_wtime();

    mt19937 rng(seed0);
    double beta = 0;
    double iter = 0;
    for (int ns=0; ns<nsweep; ns++)
    {
        iter = ns/(double) nsweep;
        beta = 0.01 + (log(n)-0.01)*iter;
        if (r) {sweep_rand(beta,rng);}
        else {sweep(beta,rng);}
    }

    return omp_get_wtime() - start_time;

}

void Ising::reset(int seed0)
{
    seed = seed0;
    rand_spins();
    update_fields();
}

void Ising::print_spins()
{
    print_mat_order(spins);
}

void Ising::print_fields()
{
    print_mat_order(fields);
}

void Ising::print_bonds()
{
    print_ten(bonds);
}

void Ising::add_bond(Row point1, Row point2, int bond)
{

    int x1 = point1[0]; int y1 = point1[1];
    int x2 = point2[0]; int y2 = point2[1];
    int dx = x2-x1; int dy = y2-y1;

    if (dx>0)
    {
        bonds[x1][y1][0] += bond;
    }
    else if (dx<0)
    {
        bonds[x2][y1][0] += bond;
    }
    else
    {
        if (dy>0)
        {
            bonds[x1][y1][1] += bond;
        }
        else
        {
            bonds[x1][y2][1] += bond;
        }
    }

}

int Ising::add_loop(Matrix cycle)
{

    int len = cycle.size()-2;
    Row lastpos = cycle[len]; Row pos = cycle[0]; Row nextpos = cycle[1];
    for (int edge=0; edge<=len; edge++)
    {
        fields[pos[0]][pos[1]] += spins[lastpos[0]][lastpos[1]] + spins[nextpos[0]][nextpos[1]];
        add_bond(pos,nextpos,1);
        lastpos = pos; pos = nextpos; nextpos = cycle[(edge+2)%(len+1)];
    }

    int index1 = rand() % (len+1);
    int index2 = (index1+1) % (len+1);
    Row pos1 = cycle[index1]; Row pos2 = cycle[index2];
    fields[pos1[0]][pos1[1]] += -2*spins[pos2[0]][pos2[1]];
    fields[pos2[0]][pos2[1]] += -2*spins[pos1[0]][pos1[1]];
    add_bond(pos1,pos2,-2);

    return len-1;

}

void Ising::flip(int i, int j, double beta, mt19937& rng)
{

    int spin = spins[i][j];
    int field = fields[i][j];
    int diff = -2*field*spin;
    double prob = exp(beta*diff);
    prob = prob/(1+prob);
    double x = u01(rng);

    if (prob + x > 1)
    {
        spins[i][j] = -spin;
        int i1 = (i+n-1)%n; int i2 = (i+1)%n;
        int j1 = (j+m-1)%m; int j2 = (j+1)%m;
        fields[i1][j] += -2*spin*bonds[i1][j][0];
        fields[i2][j] += -2*spin*bonds[i][j][0];
        fields[i][j1] += -2*spin*bonds[i][j1][1];
        fields[i][j2] += -2*spin*bonds[i][j][1];
    }

}

void Ising::flip_bound(int i, int j, double beta, mt19937& rng, bool corner)
{

    int spin = spins[i][j];
    int field = fields[i][j];
    int diff = -2*field*spin;
    double prob = exp(beta*diff);
    prob = prob/(1+prob);
    double x = u01(rng);

    if (prob + x > 1)
    {
        spins[i][j] = -spin;
        int i1 = (i-1+n)%n; int i2 = (i+1)%n;
        int j1 = j-1; int j2 = j+1;

        if ((i != 0) || (!corner))
        {
            fields[i1][j] += -2*spin*bonds[i1][j][0];
        }
        if ((i != n-1) || (!corner))
        {
            fields[i2][j] += -2*spin*bonds[i][j][0];
        }
        if (j != 0)
        {
            fields[i][j1] += -2*spin*bonds[i][j1][1];
        }
        if (j != m-1)
        {
            fields[i][j2] += -2*spin*bonds[i][j][1];
        }

    }

}

void Ising::sweep_row(int i, double beta, mt19937& rng)
{

    for (int j=0; j<m; j++)
    {
        flip(i,j,beta,rng);
    }

}

void Ising::sweep(double beta, mt19937& rng)
{

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            flip(i,j,beta,rng);
        }
    }

}

void Ising::sweep_rand(double beta, mt19937& rng)
{

    Matrix mat = randmat(n,m);

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            int x = mat[i][j];
            int ii = x/m; int jj = x%m;
            flip(ii,jj,beta,rng);
        }
    }

}
