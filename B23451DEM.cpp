// B23451
// Mukul
// HPSC_2026 Assignment-1
// ME 522
// To compile go to the folder through CLI in which this file is stored and type:
// gfortran -O3 -fopenmp B23451DEM.cpp -o B23451DEM.exe
// for running the file on 1 thread
//.\B23451DEM.exe 1
// for running the file on 2 threads
//.\B23451DEM.exe 2
// And so on.
// for fastest results, run as: 
// .\B23451DEM.exe

#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>


using namespace std;


// -------------------- Parameters ---------------------------
const double g = 9.81;        // gravity (meter/second^2)
const double dt = 1e-4;       // timestep (seconds)
const int steps = 1000000;    // number of timesteps
const double k = 1;           // spring constant (newton/meter)
const double gamma = 5e-3;    // damping coefficient(newton*second/meter)
const double radius = 1e-3;   // particle radius(meter)
const double box_x = 1.0;     // container walls with x-axis as normal, x = 0 and x = 1 meter 
const double box_y = 1.0;     // container walls with y-axis as normal, y = 0 and y = 1 meter 
const double box_z = 1.0;     // container walls with z-axis as normal, z = 0 and z = 1 meter 


// -------------------- Particle Struct --------------------
struct Particle {
    double x, y, z;     // position
    double vx, vy, vz;  // velocity
    double m;           // mass
};


// -------------------- Main --------------------------------
int main(int argc, char* argv[]) {
    int num_threads=1;

    if(argc > 1) num_threads=atoi(argv[1]);
        
    // Manage Input-Output
    #ifndef ONLINE_JUDGE
    freopen("input.txt","r",stdin); // Write the path where the input file is stored.
    #endif
    
    
    if(num_threads>1) omp_set_num_threads(num_threads); // adjust threads
    
    // Initialize particles
    int N;cin>>N;
    std::vector<Particle> particles(N);
    for(int i=0;i<N;i++){
        cin>>particles[i].x;
        cin>>particles[i].y;
        cin>>particles[i].z;
        cin>>particles[i].vx;
        cin>>particles[i].vy;
        cin>>particles[i].vz;
        cin>>particles[i].m;
    }

    // Simulation loop
    vector<double> dvxv(N,0.),dvyv(N,0.),dvzv(N,0.);
    if(num_threads==1) goto single_thread;
    for (int step = 0; step < steps; step++) { 
        
    #pragma omp parallel for 
    for (int i = 0; i < N; i++) {
        Particle& p = particles[i];

        // Gravity
        double fx = 0.0, fy = 0.0, fz = -p.m * g;

        // Wall contact (bottom only)
        if (p.z - radius < 0.0) {
            double overlap = radius - p.z;
            fz += k * overlap - gamma * p.vz;
        }

        // Wall contact (top only)
        else if(p.z + radius - box_z > 0.0) {
            double overlap = p.z + radius - box_z;
            fz-= k * overlap - gamma * p.vz;
        }

        // Wall contact (wall with positive normal as x-axis)
        if (p.x - radius < 0.0) {
            double overlap = radius - p.x;
            fx += k * overlap - gamma * p.vx;
        }

        // Wall contact (wall with normal as -ve x-axis)
        else if(p.x + radius - box_x > 0.0) {
            double overlap = p.x + radius - box_x;
            fx-= k * overlap - gamma * p.vx;
        }

        // Wall contact (wall with positive normal as y-axis)
        if (p.y - radius < 0.0) {
            double overlap = radius - p.y;
            fy += k * overlap - gamma * p.vy;
        }

        // Wall contact (wall with normal as -ve y-axis)
        else if(p.y + radius - box_y > 0.0) {
            double overlap = p.y + radius - box_y;
            fy-= k * overlap - gamma * p.vy;
        }

        // Particle–particle contacts
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            Particle& q = particles[j];

            double dx = p.x - q.x;
            double dy = p.y - q.y;
            double dz = p.z - q.z;
            double dist_sq = dx*dx + dy*dy + dz*dz;

            if (dist_sq < 4*radius*radius) {
                double dist=sqrt(dist_sq);
                double overlap = 2*radius - dist;
                double nx = dx / dist;
                double ny = dy / dist;
                double nz = dz / dist;

                double dvx = p.vx - q.vx;
                double dvy = p.vy - q.vy;
                double dvz = p.vz - q.vz;
                double vn = dvx*nx + dvy*ny + dvz*nz;

                double f_contact = k*overlap - gamma*vn;
                fx += f_contact * nx;
                fy += f_contact * ny;
                fz += f_contact * nz;
            }
        }

        // Update velocity (semi-implicit Euler)
        dvxv[i] = (fx / p.m) * dt;
        dvyv[i] = (fy / p.m) * dt;
        dvzv[i] = (fz / p.m) * dt;
    }
    
    for (int i = 0; i < N; i++){
        particles[i].vx +=dvxv[i];
        particles[i].vy +=dvyv[i];
        particles[i].vz +=dvzv[i];
        particles[i].x += particles[i].vx * dt;
        particles[i].y += particles[i].vy * dt;
        particles[i].z += particles[i].vz * dt;
    }  
    }
    goto endofcode;
    single_thread:
    for (int step = 0; step < steps; step++) { 
        
    for (int i = 0; i < N; i++) {
        Particle& p = particles[i];

        // Gravity
        double fx = 0.0, fy = 0.0, fz = -p.m * g;

        // Wall contact (bottom only)
        if (p.z - radius < 0.0) {
            double overlap = radius - p.z;
            fz += k * overlap - gamma * p.vz;
        }

        // Wall contact (top only)
        else if(p.z + radius - box_z > 0.0) {
            double overlap = p.z + radius - box_z;
            fz-= k * overlap - gamma * p.vz;
        }

        // Wall contact (wall with positive normal as x-axis)
        if (p.x - radius < 0.0) {
            double overlap = radius - p.x;
            fx += k * overlap - gamma * p.vx;
        }

        // Wall contact (wall with normal as -ve x-axis)
        else if(p.x + radius - box_x > 0.0) {
            double overlap = p.x + radius - box_x;
            fx-= k * overlap - gamma * p.vx;
        }

        // Wall contact (wall with positive normal as y-axis)
        if (p.y - radius < 0.0) {
            double overlap = radius - p.y;
            fy += k * overlap - gamma * p.vy;
        }

        // Wall contact (wall with normal as -ve y-axis)
        else if(p.y + radius - box_y > 0.0) {
            double overlap = p.y + radius - box_y;
            fy-= k * overlap - gamma * p.vy;
        }

        // Particle–particle contacts
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            Particle& q = particles[j];

            double dx = p.x - q.x;
            double dy = p.y - q.y;
            double dz = p.z - q.z;
            double dist_sq = dx*dx + dy*dy + dz*dz;

            if (dist_sq < 4*radius*radius) {
                double dist=sqrt(dist_sq);
                double overlap = 2*radius - dist;
                double nx = dx / dist;
                double ny = dy / dist;
                double nz = dz / dist;

                double dvx = p.vx - q.vx;
                double dvy = p.vy - q.vy;
                double dvz = p.vz - q.vz;
                double vn = dvx*nx + dvy*ny + dvz*nz;

                double f_contact = k*overlap - gamma*vn;
                fx += f_contact * nx;
                fy += f_contact * ny;
                fz += f_contact * nz;
            }
        }

        // Update velocity (semi-implicit Euler)
        dvxv[i] = (fx / p.m) * dt;
        dvyv[i] = (fy / p.m) * dt;
        dvzv[i] = (fz / p.m) * dt;
    }
    
    for (int i = 0; i < N; i++){
        particles[i].vx +=dvxv[i];
        particles[i].vy +=dvyv[i];
        particles[i].vz +=dvzv[i];
        particles[i].x += particles[i].vx * dt;
        particles[i].y += particles[i].vy * dt;
        particles[i].z += particles[i].vz * dt;
    }  
    }
    endofcode:
    return 0;
}
