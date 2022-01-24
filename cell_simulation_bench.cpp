/* This is the bench for simulation of growing dense active matter. [1]
author: Chengjie Luo 
email: c.luo@tue.nl
date: 2021/8/26
[1]Tjhung, Elsen, and Ludovic Berthier. "Analogies between growing dense active matter and soft driven glasses." Physical Review Research 2.4 (2020): 043334.
*/

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <chrono>
#include <random>
#include "cell_simulation.hpp"

using namespace std;

//class of 2d particles with r={x,y} and sr is the diameter.
class Particle
{
public:
    Particle()
    {
    }
    Particle(int pid, std::vector<double> pr, double psr, std::vector<double> pv = std::vector<double>(2, 0), double psv = 0, int ptype = 0)
    {
        r = pr;
        sr = psr;
        v = pv;
        sv = psv;
        id = pid;
        type = ptype;
    }

    void setPos(const double px, const double py, const double pz)
    {
        r[0] = px;
        r[1] = py;
        r[2] = pz;
    }

    void setPos(std::vector<double> pr)
    {
        r = pr;
    }

    void setId(const int pid)
    {
        id = pid;
    }

    void setType(const int ptype)
    {
        type = ptype;
    }
    std::vector<double> getPos()
    {
        return r;
    }
    int getId()
    {
        return id;
    }
    int getType()
    {
        return type;
    }

    std::vector<double> r = std::vector<double>(2, 0);
    std::vector<double> v = std::vector<double>(2, 0);
    double sr = 0;
    double sv = 0;
    int id = -1;
    int type = -1;
    double alpha = 0.01;
};

vector<Particle> particles;
double epsilon, uc, alpha_max, dt, tmax, ksi;
double tnow = 0;
long long int timestep = 0;
int Np;
double usum;

unsigned seed;
default_random_engine generator;
uniform_real_distribution<double> distribution(0.0, 2. * M_PI);

unsigned seed_alpha;
default_random_engine generator_alpha;
uniform_real_distribution<double> distribution_alpha(0.0, 1.0);

ofstream ofile_config;

double vmax = 0;
double svmax = 0;
int main(int argc, char **argv)
{
    if (argc == 1)
    {
        cout << " input should be : uc alpha_max dt tmax configfile seed seedalpha " << endl;
        return 0;
    }
    auto start = chrono::steady_clock::now();

    int iarg = 1;

    //// parse arguments and set parameters
    uc = stod(argv[iarg++]);
    cout << "uc=" << uc << endl;
    alpha_max = stod(argv[iarg++]);
    cout << "alpha_max=" << alpha_max << endl;
    dt = stod(argv[iarg++]);
    cout << "dt=" << dt << endl;
    tmax = stod(argv[iarg++]);
    cout << "tmax=" << tmax << endl;
    string configfile = argv[iarg++];
    cout << "configfile is " << configfile << endl;
    seed = stoi(argv[iarg++]);
    cout << "seed=" << seed << endl;
    seed_alpha = stoi(argv[iarg++]);
    cout << "seed_alpha=" << seed_alpha << endl;

    epsilon = 1.;
    ksi = 1.;

    generator.seed(seed);
    generator_alpha.seed(seed_alpha);
    //// initialize
    vector<double> tmpr{0, 0};
    double tmpsr = 1;
    Particle p(0, tmpr, tmpsr);
    particles.push_back(p);
    Np = 1;

    //// run simulation
    ofile_config.open(configfile);

    while (tnow < tmax)
    {
        SingleStep();
    }
    ofile_config.close();

    cout << "vmax=" << vmax << endl;
    cout << "svmax=" << svmax << endl;
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time : "
         << chrono::duration_cast<chrono::seconds>(end - start).count()
         << " sec" << endl;
}

void SingleStep()
{
    // if (timestep < 10)
    //     Printconfig();
    // else if (timestep < 100)
    // {
    //     if (timestep % 1 == 0)
    //         Printconfig();
    // }
    // else if (timestep < 1000)
    // {
    //     if (timestep % 10 == 0)
    //         Printconfig();
    // }
    // else if (timestep < 10000)
    // {
    //     if (timestep % 100 == 0)
    //         Printconfig();
    // }
    // else if (timestep < 100000)
    // {
    //     if (timestep % 1000 == 0)
    //         Printconfig();
    // }
    // else if (timestep < 1000000)
    // {
    //     if (timestep % 10000 == 0)
    //         Printconfig();
    // }
    // else if (timestep < 10000000)
    // {
    //     if (timestep % 100000 == 0)
    //         Printconfig();
    // }


    if (timestep < 10)
        Printconfig();
    else if (timestep < 100)
    {
        if (timestep % 5 == 0)
            Printconfig();
    }
    else if (timestep < 1000)
    {
        if (timestep % 50 == 0)
            Printconfig();
    }
    else if (timestep < 10000)
    {
        if (timestep % 500 == 0)
            Printconfig();
    }
    else if (timestep < 100000)
    {
        if (timestep % 5000 == 0)
            Printconfig();
    }
    else if (timestep < 1000000)
    {
        if (timestep % 50000 == 0)
            Printconfig();
    }
    else if (timestep < 10000000)
    {
        if (timestep % 500000 == 0)
            Printconfig();
    }
    else if (timestep < 100000000)
    {
        if (timestep % 5000000 == 0)
            Printconfig();
    }
    else if (timestep < 1000000000)
    {
        if (timestep % 50000000 == 0)
            Printconfig();
    }
    ++timestep;
    tnow = timestep * dt;
    Cal_v_sv();

    Cal_r_sr();
}

void Printconfig()
{
    ofile_config << "timestep=" << tnow << endl;
    ofile_config << "Np=" << Np << endl;
    for (int ip = 0; ip < Np; ip++)
    {
        ofile_config << particles[ip].r[0] << ' ' << particles[ip].r[1] << ' ' << particles[ip].sr << endl;
    }
}
void Cal_v_sv()
{

    for (int ip = 0; ip < Np; ip++)
    {
        particles[ip].v[0] = 0;
        particles[ip].v[1] = 0;
        particles[ip].sv = 0;
    }
    usum = 0;
    for (int ip = 0; ip < Np - 1; ip++)
    {
        for (int jp = ip + 1; jp < Np; jp++)
        {
            double dx = particles[ip].r[0] - particles[jp].r[0];
            double dy = particles[ip].r[1] - particles[jp].r[1];
            double smean = 0.5 * (particles[ip].sr + particles[jp].sr);
            double dr = sqrt(dx * dx + dy * dy);
            if (dr <= smean)
            {
                double prefactor = epsilon * 0.5 * (particles[ip].sr * particles[ip].sr + particles[jp].sr * particles[jp].sr);
                double fac = prefactor * (1. / smean * (1. / dr - 1. / smean));

                particles[ip].v[0] += fac * dx;
                particles[ip].v[1] += fac * dy;

                particles[jp].v[0] += (-fac * dx);
                particles[jp].v[1] += (-fac * dy);

                usum += prefactor * 0.5 * (1. - dr / smean) * (1. - dr / smean);
            }
        }
    }
    usum *= 4. / Np / M_PI;
    //ofile_config << setprecision(15) << "usum" << usum << endl;

    if (usum > uc)
    {
        for (int ip = 0; ip < Np; ip++)
        {
            particles[ip].sv = 0;
        }
    }
    else
    {
        for (int ip = 0; ip < Np; ip++)
        {
            particles[ip].sv = particles[ip].alpha * (uc - usum) / uc;
            //cout << "sv=" << particles[ip].sv << endl;
        }
    }

    for (int ip = 0; ip < Np; ip++)
    {
        if (vmax < fabs(particles[ip].v[0]))
            vmax = fabs(particles[ip].v[0]);
        if (vmax < fabs(particles[ip].v[1]))
            vmax = fabs(particles[ip].v[1]);
        if (svmax < fabs(particles[ip].sv))
            svmax = fabs(particles[ip].sv);
    }
}

void Cal_r_sr()
{
    //update x,y first
    for (int ip = 0; ip < Np; ip++)
    {
        particles[ip].r[0] += dt * particles[ip].v[0] / ksi;
        particles[ip].r[1] += dt * particles[ip].v[1] / ksi;
    }
    for (int ip = 0; ip < Np; ip++)
    {
        particles[ip].sr += dt * particles[ip].sv;
        if (particles[ip].sr > sqrt(2))
        {
            double tmpx, tmpy, newx1, newy1, newx2, newy2;
            tmpx = particles[ip].r[0];
            tmpy = particles[ip].r[1];
            double theta = 0.; //should be random number or pressure related
            theta = distribution(generator);
            //cout << "theta=" << theta << endl;
            newx1 = tmpx - (sqrt(2.) - 1.) * 0.5 * cos(theta);
            newy1 = tmpy - (sqrt(2.) - 1.) * 0.5 * sin(theta);
            newx2 = tmpx + (sqrt(2.) - 1.) * 0.5 * cos(theta);
            newy2 = tmpy + (sqrt(2.) - 1.) * 0.5 * sin(theta);

            vector<double> tmpr1{newx1, newy1}, tmpr2{newx2, newy2};
            // update left particle
            particles[ip].r = tmpr1;
            particles[ip].sr = 1;
            // update a new right particle
            Particle newp(Np, tmpr2, 1);
            newp.alpha = distribution_alpha(generator_alpha) * alpha_max;
            //cout << "alpha=" << newp.alpha << endl;
            particles.push_back(newp);
            Np = particles.size();
        }
    }
}