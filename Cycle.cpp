#include "Cycle.h"

Cycle::Cycle(int n0, int m0, int seed0)
{

    n = n0; m = m0;
    int x = randnum(0,n-1); pos.x = x;
    int y = randnum(0,n-1); pos.y = y;
    trail.push_back({x,y});
    pos.dir = 0;
    seed = seed0;

}

void Cycle::walk()
{

    mt19937 rng(seed);
    bool intersect = false;
    Row loc = {0,0};
    Matrix::iterator it;

    while (!intersect)
    {
        steprand(rng);
        loc = {pos.x,pos.y};
        it = find(trail.begin(), trail.end()-1, loc);
        if (it != trail.end()-1)
        {
            intersect = true;
        }

    }

    Matrix temp(it,trail.end());
    close = temp;

}

void Cycle::reset(int seed0)
{

    int x = randnum(0,n-1); pos.x = x;
    int y = randnum(0,n-1); pos.y = y;
    pos.dir = 0;
    trail.clear();
    trail.push_back({x,y});
    close.clear();
    seed = seed0;

}

void Cycle::step(int dir)
{

    if (dir == 1)
    {
        pos.x = pos.x - 1;
        trail.push_back({pos.x,pos.y});
        pos.dir = 1;
    }
    else if (dir == 2)
    {
        pos.y = pos.y - 1;
        trail.push_back({pos.x,pos.y});
        pos.dir = 2;
    }
    else if (dir == 3)
    {
        pos.x = pos.x + 1;
        trail.push_back({pos.x,pos.y});
        pos.dir = 3;
    }
    else if (dir == 4)
    {
        pos.y = pos.y + 1;
        trail.push_back({pos.x,pos.y});
        pos.dir = 4;
    }

}

void Cycle::steprand(mt19937& rng)
{

    vector<int> dirs;
    int dir1=0, dir2=0, c=0;
    int dir3 = (pos.dir+2)%4;

    if (pos.x == 0)
    {
        dir1 = 1;
        c++;
    }
    else if (pos.x == n-1)
    {
        dir1 = 3;
        c++;
    }

    if (pos.y == 0)
    {
        dir2 = 2;
        c++;
    }
    else if (pos.y == m-1)
    {
        dir2 = 4;
        c++;
    }

    for (int i=1; i<=4; i++)
    {
        if ( (i != dir1) && (i != dir2) && (i != dir3) )
        {
            dirs.push_back(i);
        }
    }

    int ind = uid2(rng)%(3-c);
    int dir = dirs[ind];
    step(dir);

}
