#ifndef CYCLE_H
#define CYCLE_H

#include "helper.h"
using namespace std;

typedef vector< vector<int> > Matrix;
typedef vector<int> Row;

class Cycle
{

    public:

        Cycle(int n0, int m0, int seed0);

        Matrix close;

        void walk();

        void reset(int seed0);

    private:

        int n, m, seed;
        struct Position
        {
            int x;
            int y;
            int dir;
        };
        Position pos;

        Matrix trail;

        void step(int dir);
        void steprand(mt19937& rng);

};

#endif // CYCLE_H
