#include <stdio.h>      // For input and output operations
#include <stdlib.h>     // For memory allocation and random number functions
#include <math.h>       // For mathematical functions like sqrt
#include <time.h>       // For time function (used to seed random number generator)

// Define simulation constants
#define NX 30           // Number of grid points in x-direction
#define NY 30           // Number of grid points in y-direction
#define NT 4000         // Number of time steps
#define IP 8000         // Total number of particles
#define Dd 1.0          // Grid spacing
#define Dt 0.05         // Time step size
#define TAU 0.000755    // Buoyancy term scaling factor
#define WB 0.35         // Bubble velocity
#define QG 0.001        // Gas flow rate

// Function to generate random number between -1 and 1
double rand_uniform() {
    return 2.0 * ((double)rand() / RAND_MAX) - 1.0;  // Scale and shift rand() to [-1, 1]
}

// Function to allocate 2D array and initialize with a given value
double** allocate_2d_array(int rows, int cols, double init_val) {
    double** array = (double**)malloc(rows * sizeof(double*));  // Allocate row pointers
    for (int i = 0; i < rows; i++) {
        array[i] = (double*)malloc(cols * sizeof(double));  // Allocate each row
        for (int j = 0; j < cols; j++) {
            array[i][j] = init_val;  // Initialize each element
        }
    }
    return array;  // Return pointer to the array
}

// Function to free a 2D array
void free_2d_array(double** array, int rows) {
    for (int i = 0; i < rows; i++) free(array[i]);  // Free each row
    free(array);  // Free row pointers
}

// Function to write concentration data to CSV file
void write_concentration_to_csv(double** c) {
    FILE* f = fopen("concentration.csv", "w");  // Open file for writing
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            fprintf(f, "%f", c[i][j]);  // Write concentration value
            if (i < NX - 1) fprintf(f, ",");  // Add comma if not last column
        }
        fprintf(f, "\n");  // New line after each row
    }
    fclose(f);  // Close file
}

// Function to write velocity field to CSV file
void write_velocity_to_csv(double** u, double** v) {
    FILE* f = fopen("velocity.csv", "w");  // Open file
    fprintf(f, "i,j,u,v\n");  // Write header
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            fprintf(f, "%d,%d,%f,%f\n", i, j, u[i][j], v[i][j]);  // Write grid and velocity
        }
    }
    fclose(f);  // Close file
}

// Function to write kinetic energy over time to CSV
void write_kinetic_energy_to_csv(double* KinE) {
    FILE* f = fopen("kinetic_energy.csv", "w");  // Open file
    for (int i = 0; i < NT; i++) {
        fprintf(f, "%f\n", KinE[i]);  // Write energy value for each time step
    }
    fclose(f);  // Close file
}

// Main simulation function
void simulate(double cp, double ed, double sources[][2], int num_sources, double* KinE) {
    int ipp = IP / NT;  // Particles per time step
    int max_particles = IP * num_sources;  // Max total particles, num_source means number of gas source point 
    double vol = QG * Dt / ipp;  // Volume of gas added per particle per time step

    // Allocate fields
    double **h = allocate_2d_array(NX, NY, 1.0);      // pressure head
    double **hn = allocate_2d_array(NX, NY, 1.0);     // new pressure head
    double **u = allocate_2d_array(NX, NY, 0.0);      // Velocity x
    double **un = allocate_2d_array(NX, NY, 0.0);     // New velocity x
    double **v = allocate_2d_array(NX, NY, 0.0);      // Velocity y
    double **vn = allocate_2d_array(NX, NY, 0.0);     // New velocity y
    double **c = allocate_2d_array(NX, NY, 0.0);      // Concentration
    double **co = allocate_2d_array(NX, NY, 0.0);     // Old concentration

    // Allocate particle positions
    double* x = (double*)malloc(max_particles * sizeof(double));  // Particle x positions
    double* y = (double*)malloc(max_particles * sizeof(double));  // Particle y positions
    int particle_count = 0;  // Initial particle count

    // Time loop
    for (int k = 0; k < NT; k++) {
        int ipt = (ipp * (k + 1) < max_particles) ? ipp * (k + 1) : max_particles;  // Active particles
        double D = sqrt(6.0 * ed / Dt);  // Diffusion term

        // Inject new particles at sources
        for (int s = 0; s < num_sources; s++) {
            if (k * ipp < IP) {
                for (int i = 0; i < ipp; i++) {
                    x[particle_count] = sources[s][0];  // Initial x
                    y[particle_count] = sources[s][1];  // Initial y
                    particle_count++;
                }
            }
        }

        // Update pressure 
        for (int j = 1; j < NY - 1; j++)
            for (int i = 1; i < NX - 1; i++)
                hn[i][j] = h[i][j]
                         - Dt / Dd / 2.0 * (u[i + 1][j] - u[i][j] + v[i][j + 1] - v[i][j])
                         + (c[i][j] - co[i][j]) * vol / (Dd * Dd);  // Update formula

        // Apply boundary conditions to pressure
        for (int j = 0; j < NY; j++) {
            hn[0][j] = hn[1][j]; //Sets the leftmost column to the same value as its neighbor to the right → ∂hn/∂x = 0 at left boundary.
            hn[NX - 1][j] = hn[NX - 2][j];//Sets the rightmost column to the same value as the one just left of it → ∂hn/∂x = 0 at right boundary.
        }
        for (int i = 0; i < NX; i++) {
            hn[i][0] = hn[i][1];//Sets the bottom row equal to the one above it → ∂hn/∂y = 0 at bottom boundary.
            hn[i][NY - 1] = hn[i][NY - 2];//Sets the top row equal to the one below it → ∂hn/∂y = 0 at top boundary.
        }

        // Copy updated height/pressure and concentration
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++) {
                h[i][j] = hn[i][j];
                co[i][j] = c[i][j];
            }

        // Update x-velocity
        for (int j = 1; j < NY - 1; j++)
            for (int i = 2; i < NX - 1; i++)
                un[i][j] = u[i][j]
                         - Dt / Dd * cp * cp * (h[i][j] - h[i - 1][j])
                         + ed * 5.0 * Dt / (Dd * Dd) *
                         (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - 4.0 * u[i][j]);

        // x-velocity boundary conditions
        for (int j = 0; j < NY; j++) {
            un[1][j] = -(h[1][j] - 1.0) * cp / h[1][j];
            un[NX - 2][j] = (h[NX - 3][j] - 1.0) * cp / h[NX - 3][j];
        }
        for (int i = 0; i < NX; i++) {
            un[i][0] = un[i][1];
            un[i][NY - 1] = un[i][NY - 2];
        }

        // Update y-velocity
        for (int j = 2; j < NY - 1; j++)
            for (int i = 1; i < NX - 1; i++)
                vn[i][j] = v[i][j]
                         - Dt / Dd * cp * cp * (h[i][j] - h[i][j - 1])
                         + ed * 5.0 * Dt / (Dd * Dd) *
                         (v[i][j + 1] + v[i][j - 1] + v[i - 1][j] + v[i + 1][j] - 4.0 * v[i][j])
                         + Dt * (c[i][j] + c[i][j - 1]) / 2.0 / (Dd * Dd) * TAU;

        // y-velocity boundary conditions
        for (int j = 0; j < NY; j++) {
            vn[0][j] = vn[1][j];
            vn[NX - 1][j] = vn[NX - 2][j];
        }

        // Copy new velocities
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++) {
                u[i][j] = un[i][j];
                v[i][j] = vn[i][j];
            }

        // Move particles based on velocity and diffusion
        for (int i = 0; i < ipt; i++) {
            int X = (int)(x[i] / Dd);
            int Y = (int)(y[i] / Dd);
            double uu = 0.0, vv = 0.0;
            if (X >= 1 && X < NX - 1 && Y >= 1 && Y < NY - 1) {
                uu = u[X][Y] + (u[X + 1][Y] - u[X][Y]) * (x[i] - X * Dd) / Dd;
                vv = v[X][Y] + (v[X][Y + 1] - v[X][Y]) * (y[i] - Y * Dd) / Dd;
            }
            x[i] += Dt * (uu + rand_uniform() * D);
            y[i] += Dt * (vv + rand_uniform() * D + WB);
        }

        // Reset concentration
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
                c[i][j] = 0.0;

        // Recalculate concentration from particles
        for (int i = 0; i < ipt; i++) {
            int X = (int)(x[i] / Dd);
            int Y = (int)(y[i] / Dd);
            if (X >= 1 && X < NX - 1 && Y >= 1 && Y < NY - 1) {
                c[X][Y] += 1.0;
            }
        }

        // Compute kinetic energy
        double ke = 0.0;
        for (int i = 1; i < NX - 1; i++)
            for (int j = 1; j < NY - 1; j++)
                ke += u[i][j] * u[i][j] + v[i][j] * v[i][j];
        KinE[k] = ke;
    }

    // Write results to files
    write_concentration_to_csv(c);
    write_velocity_to_csv(u, v);
    write_kinetic_energy_to_csv(KinE);

    // Free memory
    free_2d_array(h, NX);
    free_2d_array(hn, NX);
    free_2d_array(u, NX);
    free_2d_array(un, NX);
    free_2d_array(v, NX);
    free_2d_array(vn, NX);
    free_2d_array(c, NX);
    free_2d_array(co, NX);
    free(x);
    free(y);
}

// Main function
int main() {
    srand(time(NULL));  // Seed random number generator

    double KinE[NT];  // Array to store kinetic energy
    double sources[1][2] = {{14.5, 1.5}};  // Define particle source

    simulate(4.0, 0.1, sources, 1, KinE);  // Run the simulation this is for base case

    printf("Simulation complete. CSV files written.\n");  // Notify user
    return 0;  // Exit program
}
