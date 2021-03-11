#include "subIsing.h"

subIsing::subIsing(Ising &lattice)
{

    int nthreads = lattice.nthreads;

    Row npart = getpart(nthreads);
    k1 = npart[0]; k2 = npart[1];
    if (k1 > 1) {corner = true;}

    n = lattice.n; m = lattice.m;

    Row nk1 = part(n,k1);
    Row nk2 = part(m,k2);

    for (int i=0; i<k1; i++)
    {
        Rowlattice rowlattice;
        Row nrow; Row mrow;
        for (int j=0; j<k2; j++)
        {
            int i1 = 0; int j1 = 0;
            if (i != 0) {i1 = nk1[i-1]+1;}
            if (j != 0) {j1 = nk2[j-1]+1;}
            int i2 = nk1[i]+1; int j2 = nk2[j]+1;
            Matrix subspins = submat(lattice.spins, i1, i2, j1, j2);
            Tensor subbonds = subten(lattice.bonds, i1, i2, j1, j2);
            Matrix subfields = submat(lattice.fields, i1, i2, j1, j2);
            rowlattice.push_back(Ising(subspins,subbonds,subfields));
            nrow.push_back(i2-i1); mrow.push_back(j2-j1);
        }
        nmat.push_back(nrow); mmat.push_back(mrow);
        matlattice.push_back(rowlattice);
    }

}

void subIsing::update_energy()
{

    auto_exchange();
    int temp = 0;
    for (int i=0; i<k1; i++)
    {
        for (int j=0; j<k2; j++)
        {
            matlattice[i][j].update_energy();
            temp += matlattice[i][j].energy;
        }
    }
    energy = temp;

}

void subIsing::print_spins(int i, int j)
{
    matlattice[i][j].print_spins();
}

void subIsing::print_fields(int i, int j)
{
    matlattice[i][j].print_fields();
}

void subIsing::print_bonds(int i, int j)
{
    matlattice[i][j].print_bonds();
}

void subIsing::print_inout(int i, int j)
{

    Ising& lattice = matlattice[i][j];
    Matrix in = {lattice.in_left,lattice.in_right,lattice.in_down,lattice.in_up};
    Matrix out = {lattice.out_left,lattice.out_right,lattice.out_down,lattice.out_up};
    print_mat(in);
    print_mat(out);

}

void subIsing::combine()
{

    spins = mat_empty(n,m);
    fields = mat_empty(n,m);
    bonds = ten_empty(n,m,2);

    int ni = 0;
    int n0 = 0;
    int m0 = 0;
    for (int ki=0; ki<k1; ki++)
    {
        int nj = 0;
        for (int kj=0; kj<k2; kj++)
        {
            Ising& lattice = matlattice[ki][kj];
            Matrix spins0 = lattice.spins;
            Matrix fields0 = lattice.fields;
            Tensor bonds0 = lattice.bonds;
            n0 = lattice.n; m0 = lattice.m;
            for (int i=0; i<n0; i++)
            {
                int ii = ni+i;
                for (int j=0; j<m0; j++)
                {
                    int jj = nj+j;
                    spins[ii][jj] = spins0[i][j];
                    fields[ii][jj] = fields0[i][j];
                    bonds[ii][jj][0] = bonds0[i][j][0];
                    bonds[ii][jj][1] = bonds0[i][j][1];
                }
            }
            nj += m0;
        }
        ni += n0;
    }

}

void subIsing::print_combine_spins()
{
    combine();
    print_mat_order(spins);
}

void subIsing::print_combine_fields()
{
    combine();
    print_mat_order(fields);
}

void subIsing::separate_bonds()
{

    for (int i=0; i<k1; i++)
    {
        for (int j=0; j<k2; j++)
        {
            matlattice[i][j].separate_bonds();
        }
    }

}

void subIsing::contact(Ising &mid, Ising &down, Ising &up)
{

    mid.in_down = getrow(mid.spins, 2);
    mid.in_up = getrow(mid.spins, 4);

    mid.out_down = getrow(down.spins, 4);
    mid.out_up = getrow(up.spins, 2);

}

void subIsing::contact(Ising &mid, Ising &left, Ising &right, Ising &down, Ising &up)
{

    mid.in_left = getrow(mid.spins, 1);
    mid.in_right = getrow(mid.spins, 3);
    mid.in_down = getrow(mid.spins, 2);
    mid.in_up = getrow(mid.spins, 4);

    mid.out_left = getrow(left.spins, 3);
    mid.out_right = getrow(right.spins, 1);
    mid.out_down = getrow(down.spins, 4);
    mid.out_up = getrow(up.spins, 2);

}

void subIsing::auto_contact()
{

    if (k1 == 1)
    {
        for (int j=0; j<k2; j++)
        {
            int j1 = (j-1+k2)%k2; int j2 = (j+1)%k2;
            contact(matlattice[0][j],matlattice[0][j1],matlattice[0][j2]);
        }
    }
    else
    {
        for (int i=0; i<k1; i++)
        {
            for (int j=0; j<k2; j++)
            {
                int i1 = (i-1+k1)%k1; int i2 = (i+1)%k1;
                int j1 = (j-1+k2)%k2; int j2 = (j+1)%k2;
                contact(matlattice[i][j],matlattice[i1][j],matlattice[i2][j],matlattice[i][j1],matlattice[i][j2]);
            }
        }
    }

}

void subIsing::exchange_fields(Ising &mid, Ising &down, Ising &up)
{

    int n = mid.n; int m = mid.m;

    Row downdiff = rowdiff(down.in_up, mid.out_down);
    Row updiff = rowdiff(up.in_down, mid.out_up);

    for (int i=0; i<n; i++)
    {
        mid.fields[i][0] += downdiff[i]*down.b_up[i];
        mid.fields[i][m-1] += updiff[i]*mid.b_up[i];
    }

    mid.out_down = down.in_up;
    mid.out_up = up.in_down;


}

void subIsing::exchange_fields(Ising &mid, Ising &left, Ising &right, Ising &down, Ising &up)
{

    int n = mid.n; int m = mid.m;

    Row leftdiff = rowdiff(left.in_right, mid.out_left);
    Row rightdiff = rowdiff(right.in_left, mid.out_right);
    Row downdiff = rowdiff(down.in_up, mid.out_down);
    Row updiff = rowdiff(up.in_down, mid.out_up);

    for (int j=0; j<m; j++)
    {
        mid.fields[0][j] += leftdiff[j]*left.b_right[j];
        mid.fields[n-1][j] += rightdiff[j]*mid.b_right[j];
    }

    for (int i=0; i<n; i++)
    {
        mid.fields[i][0] += downdiff[i]*down.b_up[i];
        mid.fields[i][m-1] += updiff[i]*mid.b_up[i];
    }

    mid.out_left = left.in_right;
    mid.out_right = right.in_left;
    mid.out_down = down.in_up;
    mid.out_up = up.in_down;

}

void subIsing::auto_exchange()
{

    if (k1 == 1)
    {
        for (int j=0; j<k2; j++)
        {
            int j1 = (j-1+k2)%k2; int j2 = (j+1)%k2;
            exchange_fields(matlattice[0][j],matlattice[0][j1],matlattice[0][j2]);
        }
    }
    else
    {
        for (int i=0; i<k1; i++)
        {
            for (int j=0; j<k2; j++)
            {
                int i1 = (i-1+k1)%k1; int i2 = (i+1)%k1;
                int j1 = (j-1+k2)%k2; int j2 = (j+1)%k2;
                exchange_fields(matlattice[i][j],matlattice[i1][j],matlattice[i2][j],matlattice[i][j1],matlattice[i][j2]);
            }
        }
    }

}

double subIsing::sweeps(int nsweep, int seed1)
{

    double start_time = omp_get_wtime();

    separate_bonds();
    auto_contact();

    double beta, iter;
    for (int ns=0; ns<nsweep; ns++)
    {
        iter = (double)(ns)/(nsweep-1);
        beta = 0.01 + (log(n)-0.01)*iter;

        for (int i=0; i<k1; i++)
        {
            for (int j=0; j<k2; j++)
            {
                int seed0 = (ns*k1*k2 + i*k2 + j) % 2^32;
                matlattice[i][j].sweep_bound(beta, seed0+seed1, corner);
            }
        }

        auto_exchange();
    }

    return omp_get_wtime() - start_time;

}

double subIsing::sweeps_parallel(int nsweep, int seed1)
{

    double start_time = omp_get_wtime();

    separate_bonds();
    auto_contact();

    double beta, iter;
    for (int ns=0; ns<nsweep; ns++)
    {
        iter = (double)(ns)/(nsweep-1);
        beta = 0.01 + (log(n)-0.01)*iter;

        #pragma omp parallel
        {
            int id = omp_get_thread_num();
            int i = id/k2; int j = id%k2;
            int seed0 = (ns*k1*k2 + id) % 2^32;
            matlattice[i][j].sweep_bound(beta, seed0+seed1, corner);
        }

        auto_exchange();
    }

    return omp_get_wtime() - start_time;

}
