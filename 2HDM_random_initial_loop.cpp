//==============================================================================
// 2HDM_single.cpp
// Description: Serial (non-MPI) simulation of domain wall dynamics in the
//              Z₂-symmetric two-Higgs-doublet model (2HDM). Uses 4th-order
//              finite differences and sponge layer absorbing boundaries.
//==============================================================================

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <vector>
#include <random>
#include <sstream>
#include <string>

using namespace std;

const double pi = 4.0 * atan(1.0);

// Simulation parameters (can be adjusted)
const int nts = 2;                  // Number of time levels stored in field arrays (not used)
const long long int nx = 128;      // Grid points in x-direction
const long long int ny = 128;      // Grid points in y-direction
const long long int nz = 1;        // Set to 1 for 2D simulation
const long long int nPos = nx * ny * nz;   // Total number of grid points

const int radius2 = nx / 2 - 2;    // Radius of the active simulation region

const double dx = 0.5;            // Grid spacing in x [1/m_h]
const double dy = 0.5;            // Grid spacing in y
const double dz = 1;            // Grid spacing in z (only used in 3D)
const double dt = 0.1;             // Time step [1/(c m_h)]

const int nt = static_cast<int>(15 * (nx * dx / (2 * dt)));   // Number of time steps (light crossing time)

const bool count_global = true;    // Compute global energy and wall area
const bool makeGif = false;         // Output snapshots for GIF generation
const int saveFreq = 20;            // Frequency of snapshot output

const string dir_path = "Data/";   // Output directory
const int countRate = 20;          // Progress reporting interval

const int nb_fields = 8;           // Number of real scalar fields (two doublets)

// 2HDM parameters (see Appendix 2 of the paper)
const long double m_h = 125;       // Higgs mass [GeV]
const long double V_sm = 246;      // Vacuum expectation value [GeV]
const long double m_H = 200;       // Heavy Higgs mass [GeV]
const long double m_A = 200;       // Pseudoscalar mass [GeV]
const long double m_H_pm = 200;    // Charged Higgs mass [GeV]

// Dimensionless parameters (set by the above masses and mixing angles)
const long double M_h = m_h / m_h;          // 1
const long double v_sm = V_sm / V_sm;       // 1
const long double M_H = m_H / m_h;
const long double M_A = m_A / m_h;
const long double M_H_pm = m_H_pm / m_h;

const long double a = atan(0.85);              // Mixing angle α = π/4
const long double b = atan(0.85);              // Mixing angle β = π/4
const long double s_a = sin(a), c_a = cos(a);
const long double s_b = sin(b), c_b = cos(b);
const long double t_b = tan(b), ct_b = 1.0 / t_b;

// Dimensionless potential coefficients (scaled by m_h and v_sm)
const long double mu_1_sq = (1.0 / (M_h*M_h)) * 0.5 * ((M_h*M_h * c_a*c_a) + (M_H*M_H * s_a*s_a) + ((M_h*M_h - M_H*M_H) * c_a * s_a * t_b));
const long double mu_2_sq = (1.0 / (M_h*M_h)) * 0.5 * ((M_h*M_h * s_a*s_a) + (M_H*M_H * c_a*c_a) + ((M_h*M_h - M_H*M_H) * c_a * s_a * ct_b));

const double lambda_1 = (v_sm*v_sm / (M_h*M_h)) * (M_h*M_h * c_a*c_a + M_H*M_H * s_a*s_a) / (2 * c_b*c_b * v_sm*v_sm);
const double lambda_2 = (v_sm*v_sm / (M_h*M_h)) * (M_h*M_h * s_a*s_a + M_H*M_H * c_a*c_a) / (2 * s_b*s_b * v_sm*v_sm);
const double lambda_3 = (v_sm*v_sm / (M_h*M_h)) * ((M_h*M_h - M_H*M_H) * c_a * s_a + 2 * M_H_pm*M_H_pm * c_b * s_b) / (c_b * s_b * v_sm*v_sm);
const double l4_m_l5 = (v_sm*v_sm / (M_h*M_h)) * (-2 * M_H_pm*M_H_pm) / (v_sm*v_sm);
const double l4_p_l5 = (v_sm*v_sm / (M_h*M_h)) * (2 * (M_A*M_A - M_H_pm*M_H_pm)) / (v_sm*v_sm);


// Damping Paramaters:
const int damped_nt = 300; // Number of time steps for which damping is imposed. Useful for random initial conditions
const double dampFac = 1.0; // Magnitude of damping term, unclear how strong to make this

//-----------------------------------------------------------------------------
// Main simulation driver
//-----------------------------------------------------------------------------
int main() {
    // Read oscillon profile data (produced by function_shape_2HDM.py)
    vector<vector<double>> data;

    // Output tag for file names
    int seed;
    for (seed = 1; seed < 101; seed++) {
        const string outTag = "_2d_2hdm_nx_0_" + to_string(nx) + "_nt" + to_string(nt) + "_seed" + to_string(seed);

    // Allocate field arrays: two time levels, plus derivatives
    vector<vector<double>> fields_now(nb_fields, vector<double>(nPos, 0.0));
    vector<vector<double>> fields_past(nb_fields, vector<double>(nPos, 0.0));
    vector<vector<double>> fields_dx(nb_fields, vector<double>(nPos, 0.0));
    vector<vector<double>> fields_dy(nb_fields, vector<double>(nPos, 0.0));
    vector<vector<double>> fields_dz(nb_fields, vector<double>(nPos, 0.0));
    vector<vector<double>> fields_dt(nb_fields, vector<double>(nPos, 0.0));

    // Independent random number generators for each field
    mt19937 gens[8];
    for (int i = 0; i < 8; ++i) gens[i] = mt19937(seed + i);
    uniform_real_distribution<double> dist(-1.0, 1.0);

    // Initialize fields with random values in [-1, 1]
    for (long long int idx = 0; idx < nPos; ++idx) {
        fields_now[0][idx] = dist(gens[0]);
        fields_now[1][idx] = dist(gens[1]);
        fields_now[2][idx] = dist(gens[2]);
        fields_now[3][idx] = dist(gens[3]);
        fields_now[4][idx] = dist(gens[4]);
        fields_now[5][idx] = dist(gens[5]);
        fields_now[6][idx] = dist(gens[6]);
        fields_now[7][idx] = dist(gens[7]);

        // Past time level equals the current (zero velocity)
        for (int f = 0; f < 8; ++f) fields_past[f][idx] = fields_now[f][idx];
    }

    struct timeval start, end;
    gettimeofday(&start, NULL);
    cout << "Initial data generated. Starting time evolution..." << endl;

    // Open output file for global quantities
    ofstream valsPerLoop;
    if (count_global) {
        string valsPerLoopPath = dir_path + "valsPerLoop" + outTag + ".txt";
        valsPerLoop.open(valsPerLoopPath.c_str());
        if (!valsPerLoop.is_open()) cerr << "Error: Could not open " << valsPerLoopPath << endl;
    }

    // Time evolution loop
    for (int TimeStep = 0; TimeStep <= nt; ++TimeStep) {
        // Friction coefficient (currently zero)
        double fric;
        if (TimeStep < damped_nt) {fric=dampFac;}
        else{fric = 0.0;}

        // Accumulators for energy and wall area
        double totalEnergy = 0.0, totalKE = 0.0, totalPE = 0.0, totalTE = 0.0;
        double totalADW = 0.0;

        // Evolve all points
        for (long long int idx = 0; idx < nPos; ++idx) {
            // Coordinates for absorbing boundary
            double xval = (idx / (ny * nz)) - nx/2 + 0.5;
            double yval = ((idx % (ny * nz)) / nz) - ny/2 + 0.5;
            double zval = (idx % nz) - nz/2 + 0.5;
            double r_curr = sqrt(xval*xval + yval*yval + zval*zval);

            // Periodic neighbour indices (only used if inside radius2)
            auto wrap = [&](long long int i) -> long long int {
                return (i + nPos) % nPos;
            };

            // x direction (neighbours along x)
            long long int imx = wrap(idx - ny * nz);
            long long int ipx = wrap(idx + ny * nz);
            long long int imxx = wrap(idx - 2 * ny * nz);
            long long int ipxx = wrap(idx + 2 * ny * nz);

            // y direction (neighbours along y)
            long long int imy = wrap(idx - nz);
            long long int ipy = wrap(idx + nz);
            long long int imyy = wrap(idx - 2 * nz);
            long long int ipyy = wrap(idx + 2 * nz);

            // z direction: in 2D, all z neighbours point to the same point
            long long int imz, ipz, imzz, ipzz;
            if (nz == 1) {
                imz = ipz = imzz = ipzz = idx;
            } else {
                imz = wrap(idx - 1);
                ipz = wrap(idx + 1);
                imzz = wrap(idx - 2);
                ipzz = wrap(idx + 2);
            }

            // Current field values
            double f0 = fields_now[0][idx];
            double f1 = fields_now[1][idx];
            double f2 = fields_now[2][idx];
            double f3 = fields_now[3][idx];
            double f4 = fields_now[4][idx];
            double f5 = fields_now[5][idx];
            double f6 = fields_now[6][idx];
            double f7 = fields_now[7][idx];

            // Temporary storage for derivatives
            vector<double> fieldt(nb_fields), fieldx(nb_fields), fieldy(nb_fields), fieldz(nb_fields);
            vector<double> fieldtt(nb_fields);

            // Compute derivatives for all components (will be used only if inside radius2)
            for (int comp = 0; comp < nb_fields; ++comp) {
                // 4th-order Laplacian
                double fieldxx = (16.0 * (fields_now[comp][ipx] + fields_now[comp][imx])
                                - 30.0 * fields_now[comp][idx]
                                - fields_now[comp][ipxx] - fields_now[comp][imxx]) / (12.0 * dx * dx);
                double fieldyy = (16.0 * (fields_now[comp][ipy] + fields_now[comp][imy])
                                - 30.0 * fields_now[comp][idx]
                                - fields_now[comp][ipyy] - fields_now[comp][imyy]) / (12.0 * dy * dy);
                double fieldzz = (16.0 * (fields_now[comp][ipz] + fields_now[comp][imz])
                                - 30.0 * fields_now[comp][idx]
                                - fields_now[comp][ipzz] - fields_now[comp][imzz]) / (12.0 * dz * dz);

                // First time derivative
                fieldt[comp] = (fields_now[comp][idx] - fields_past[comp][idx]) / dt;

                // Spatial derivatives for energy
                fieldx[comp] = (fields_now[comp][ipx] - fields_now[comp][imx]) / (2.0*dx);
                fieldy[comp] = (fields_now[comp][ipy] - fields_now[comp][imy]) / (2.0*dy);
                fieldz[comp] = (fields_now[comp][ipz] - fields_now[comp][imz]) / (2.0*dz);

                // Equation of motion (second time derivative) – base part common to all
                double base1 = fieldxx + fieldyy + fieldzz - fric * fieldt[comp]
                            - ( - mu_1_sq * fields_now[comp][idx]
                              + 2.0 * lambda_1 * fields_now[comp][idx] * (f0*f0 + f1*f1 + f2*f2 + f3*f3)
                            + lambda_3 * fields_now[comp][idx] * (f4*f4 + f5*f5 + f6*f6 + f7*f7) );
                double base2 = fieldxx + fieldyy + fieldzz - fric * fieldt[comp]
                                - ( - mu_2_sq * fields_now[comp][idx]
                                + lambda_3 * fields_now[comp][idx] * (f0*f0 + f1*f1 + f2*f2 + f3*f3)
                                + 2.0 * lambda_2 * fields_now[comp][idx] * (f4*f4 + f5*f5 + f6*f6 + f7*f7) );

                // Additional interaction terms (component-specific)
                if (comp == 0) {
                    fieldtt[comp] = base1
                                  - l4_m_l5 * f4 * (f0*f4 + f1*f5 + f2*f6 + f3*f7)
                                  - l4_p_l5 * f5 * (f0*f5 - f1*f4 + f2*f7 - f3*f6);
                } else if (comp == 1) {
                    fieldtt[comp] = base1
                                  - l4_m_l5 * f5 * (f0*f4 + f1*f5 + f2*f6 + f3*f7)
                                  - l4_p_l5 * (-f4) * (f0*f5 - f1*f4 + f2*f7 - f3*f6);
                } else if (comp == 2) {
                    fieldtt[comp] = base1
                                  - l4_m_l5 * f6 * (f0*f4 + f1*f5 + f2*f6 + f3*f7)
                                  - l4_p_l5 * f7 * (f0*f5 - f1*f4 + f2*f7 - f3*f6);
                } else if (comp == 3) {
                    fieldtt[comp] = base1
                                  - l4_m_l5 * f7 * (f0*f4 + f1*f5 + f2*f6 + f3*f7)
                                  - l4_p_l5 * (-f6) * (f0*f5 - f1*f4 + f2*f7 - f3*f6);
                } else if (comp == 4) {
                    fieldtt[comp] = base2
                                  - l4_m_l5 * f0 * (f0*f4 + f1*f5 + f2*f6 + f3*f7)
                                  - l4_p_l5 * (-f1) * (f0*f5 - f1*f4 + f2*f7 - f3*f6);
                } else if (comp == 5) {
                    fieldtt[comp] = base2
                                  - l4_m_l5 * f1 * (f0*f4 + f1*f5 + f2*f6 + f3*f7)
                                  - l4_p_l5 * f0 * (f0*f5 - f1*f4 + f2*f7 - f3*f6);
                } else if (comp == 6) {
                    fieldtt[comp] = base2
                                  - l4_m_l5 * f2 * (f0*f4 + f1*f5 + f2*f6 + f3*f7)
                                  - l4_p_l5 * (-f3) * (f0*f5 - f1*f4 + f2*f7 - f3*f6);
                } else if (comp == 7) {
                    fieldtt[comp] = base2
                                  - l4_m_l5 * f3 * (f0*f4 + f1*f5 + f2*f6 + f3*f7)
                                  - l4_p_l5 * f2 * (f0*f5 - f1*f4 + f2*f7 - f3*f6);
                }
            }

            // Update fields
            for (int comp = 0; comp < nb_fields; ++comp) {
                double new_val = 2.0 * fields_now[comp][idx] - fields_past[comp][idx] + dt * dt * fieldtt[comp];
                fields_past[comp][idx] = new_val;
                fields_dt[comp][idx] = fieldt[comp];
                fields_dx[comp][idx] = fieldx[comp];
                fields_dy[comp][idx] = fieldy[comp];
                fields_dz[comp][idx] = fieldz[comp];
                }
            

            // Energy and wall area accumulation
            for (int comp = 0; comp < nb_fields; ++comp) {
                double k = fieldt[comp];
                double dx_ = fieldx[comp];
                double dy_ = fieldy[comp];
                double dz_ = fieldz[comp];
                totalKE += (k*k)*dx*dy*dz;
                totalTE += (dx_*dx_ + dy_*dy_ + dz_*dz_)*dx*dy*dz;
                totalEnergy += (k*k + dx_*dx_ + dy_*dy_ + dz_*dz_)*dx*dy*dz;
            }
            double pot = (-mu_1_sq * (f0*f0 + f1*f1 + f2*f2 + f3*f3)
                        - mu_2_sq * (f4*f4 + f5*f5 + f6*f6 + f7*f7)
                        + lambda_1 * pow(f0*f0 + f1*f1 + f2*f2 + f3*f3, 2)
                        + lambda_2 * pow(f4*f4 + f5*f5 + f6*f6 + f7*f7, 2)
                        + lambda_3 * (f0*f0 + f1*f1 + f2*f2 + f3*f3) * (f4*f4 + f5*f5 + f6*f6 + f7*f7)
                        + l4_m_l5 * pow(f0*f4 + f1*f5 + f2*f6 + f3*f7, 2)
                        + l4_p_l5 * pow(f0*f5 - f1*f4 + f2*f7 - f3*f6, 2) + 0.125)*dx*dy*dz;
            totalPE += pot;
            totalEnergy += pot;

            // Wall area (using R1 sign change)
            double R1_i = 2.0 * (f0*f4 + f1*f5 + f2*f6 + f3*f7);
            double R1_ipx = 2.0 * (fields_now[0][ipx]*fields_now[4][ipx] +
                                    fields_now[1][ipx]*fields_now[5][ipx] +
                                    fields_now[2][ipx]*fields_now[6][ipx] +
                                    fields_now[3][ipx]*fields_now[7][ipx]);
            double R1_ipy = 2.0 * (fields_now[0][ipy]*fields_now[4][ipy] +
                                    fields_now[1][ipy]*fields_now[5][ipy] +
                                    fields_now[2][ipy]*fields_now[6][ipy] +
                                    fields_now[3][ipy]*fields_now[7][ipy]);
            double R1_ipz = 2.0 * (fields_now[0][ipz]*fields_now[4][ipz] +
                                    fields_now[1][ipz]*fields_now[5][ipz] +
                                    fields_now[2][ipz]*fields_now[6][ipz] +
                                    fields_now[3][ipz]*fields_now[7][ipz]);
            double R1_imx = 2.0 * (fields_now[0][imx]*fields_now[4][imx] +
                                    fields_now[1][imx]*fields_now[5][imx] +
                                    fields_now[2][imx]*fields_now[6][imx] +
                                    fields_now[3][imx]*fields_now[7][imx]);
            double R1_imy = 2.0 * (fields_now[0][imy]*fields_now[4][imy] +
                                    fields_now[1][imy]*fields_now[5][imy] +
                                    fields_now[2][imy]*fields_now[6][imy] +
                                    fields_now[3][imy]*fields_now[7][imy]);
            double R1_imz = 2.0 * (fields_now[0][imz]*fields_now[4][imz] +
                                    fields_now[1][imz]*fields_now[5][imz] +
                                    fields_now[2][imz]*fields_now[6][imz] +
                                    fields_now[3][imz]*fields_now[7][imz]);

            double R1x = (R1_ipx - R1_imx) / (2.0*dx);
            double R1y = (R1_ipy - R1_imy) / (2.0*dy);
            double R1z = (R1_ipz - R1_imz) / (2.0*dz);
            double norm = sqrt(R1x*R1x + R1y*R1y + R1z*R1z);
            double denom = abs(R1x) + abs(R1y) + abs(R1z);

            if (R1_i * R1_ipx < 0)
                totalADW += dy * dz * norm / denom;
            if (R1_i * R1_ipy < 0)
                totalADW += dx * dz * norm / denom;
            if (R1_i * R1_ipz < 0)
                totalADW += dx * dy * norm / denom;
        }
    

        // Write global quantities to file
        if (count_global && TimeStep == 0 && valsPerLoop.is_open())
            valsPerLoop << "Energy Energy-KE Energy-tention Energy-potential ADW_Full" << endl;
        if (count_global && valsPerLoop.is_open())
            valsPerLoop << totalEnergy << " " << totalKE << " " << totalTE << " " << totalPE << " " << totalADW << endl;

        // Swap time levels for next iteration
        fields_now.swap(fields_past);

        // Output snapshots for GIF
        if (makeGif && TimeStep % saveFreq == 0) {
            string TimeStepPath = dir_path + "fields_timestep=" + to_string(TimeStep) + outTag + ".txt";
            ofstream Gif(TimeStepPath.c_str());
            Gif << "R0 R1 R2 R3 R4 R5 n1 n2 n3 psi1 psi2 psi3 psi4 psi5 psi6 psi7 psi8 spare1 spare2 spare3 E KE PE" << endl;

            for (long long int j = 0; j < nPos; ++j) {
                double psi1 = fields_now[0][j];
                double psi2 = fields_now[1][j];
                double psi3 = fields_now[2][j];
                double psi4 = fields_now[3][j];
                double psi5 = fields_now[4][j];
                double psi6 = fields_now[5][j];
                double psi7 = fields_now[6][j];
                double psi8 = fields_now[7][j];

                double psi1t = fields_dt[0][j];
                double psi2t = fields_dt[1][j];
                double psi3t = fields_dt[2][j];
                double psi4t = fields_dt[3][j];
                double psi5t = fields_dt[4][j];
                double psi6t = fields_dt[5][j];
                double psi7t = fields_dt[6][j];
                double psi8t = fields_dt[7][j];

                double psi1x = fields_dx[0][j];
                double psi2x = fields_dx[1][j];
                double psi3x = fields_dx[2][j];
                double psi4x = fields_dx[3][j];
                double psi5x = fields_dx[4][j];
                double psi6x = fields_dx[5][j];
                double psi7x = fields_dx[6][j];
                double psi8x = fields_dx[7][j];

                double psi1y = fields_dy[0][j];
                double psi2y = fields_dy[1][j];
                double psi3y = fields_dy[2][j];
                double psi4y = fields_dy[3][j];
                double psi5y = fields_dy[4][j];
                double psi6y = fields_dy[5][j];
                double psi7y = fields_dy[6][j];
                double psi8y = fields_dy[7][j];

                double psi1z = fields_dz[0][j];
                double psi2z = fields_dz[1][j];
                double psi3z = fields_dz[2][j];
                double psi4z = fields_dz[3][j];
                double psi5z = fields_dz[4][j];
                double psi6z = fields_dz[5][j];
                double psi7z = fields_dz[6][j];
                double psi8z = fields_dz[7][j];

                double R0 = psi1*psi1 + psi2*psi2 + psi3*psi3 + psi4*psi4 + psi5*psi5 + psi6*psi6 + psi7*psi7 + psi8*psi8;
                double R1 = 2.0 * (psi1*psi5 + psi2*psi6 + psi3*psi7 + psi4*psi8);
                double R2 = 2.0 * (psi1*psi6 + psi3*psi8 - psi2*psi5 - psi4*psi7);
                double R3 = (psi1*psi1+psi2*psi2+psi3*psi3+psi4*psi4) - (psi5*psi5+psi6*psi6+psi7*psi7+psi8*psi8);
                double R4 = 2.0 * (psi1*psi7 - psi2*psi8 - psi3*psi5 + psi4*psi6);
                double R5 = 2.0 * (psi1*psi8 + psi2*psi7 - psi3*psi6 - psi4*psi5);

                double n1 = -2.0 * (psi1*psi3 + psi2*psi4 + psi5*psi7 + psi6*psi8);
                double n2 = -2.0 * (psi1*psi4 - psi2*psi3 + psi5*psi8 - psi6*psi7);
                double n3 = -(psi1*psi1 + psi2*psi2 - psi3*psi3 - psi4*psi4 + psi5*psi5 + psi6*psi6 - psi7*psi7 - psi8*psi8);

                double KE = psi1t*psi1t + psi2t*psi2t + psi3t*psi3t + psi4t*psi4t +
                            psi5t*psi5t + psi6t*psi6t + psi7t*psi7t + psi8t*psi8t;

                double TE = psi1x*psi1x + psi2x*psi2x + psi3x*psi3x + psi4x*psi4x +
                            psi5x*psi5x + psi6x*psi6x + psi7x*psi7x + psi8x*psi8x +
                            psi1y*psi1y + psi2y*psi2y + psi3y*psi3y + psi4y*psi4y +
                            psi5y*psi5y + psi6y*psi6y + psi7y*psi7y + psi8y*psi8y +
                            psi1z*psi1z + psi2z*psi2z + psi3z*psi3z + psi4z*psi4z +
                            psi5z*psi5z + psi6z*psi6z + psi7z*psi7z + psi8z*psi8z;

                double PE = 0.125 - mu_1_sq*(psi1*psi1+psi2*psi2+psi3*psi3+psi4*psi4)
                                  - mu_2_sq*(psi5*psi5+psi6*psi6+psi7*psi7+psi8*psi8)
                                  + lambda_1 * pow(psi1*psi1+psi2*psi2+psi3*psi3+psi4*psi4,2)
                                  + lambda_2 * pow(psi5*psi5+psi6*psi6+psi7*psi7+psi8*psi8,2)
                                  + lambda_3 * (psi1*psi1+psi2*psi2+psi3*psi3+psi4*psi4)*(psi5*psi5+psi6*psi6+psi7*psi7+psi8*psi8)
                                  + l4_m_l5 * pow(psi1*psi5+psi2*psi6+psi3*psi7+psi4*psi8,2)
                                  + l4_p_l5 * pow(psi1*psi6 - psi2*psi5 + psi3*psi8 - psi4*psi7,2);
                double E = KE + PE + TE;

                double spare1, spare2, spare3 =0;
                Gif << R0 << " " << R1 << " " << R2 << " " << R3 << " " << R4 << " "
                    << R5 << " " << n1 << " " << n2 << " " << n3 << " "
                    << psi1 << " " << psi2 << " " << psi3 << " " << psi4 << " "
                    << psi5 << " " << psi6 << " " << psi7 << " " << psi8 << " "
                    << spare1 << " " << spare2 << " " << spare3 << " " << E << " " << KE << " " << PE << endl;
            }
            Gif.close();
        }

        // Progress indicator
        if (TimeStep % countRate == 0)
            cout << "\rFile:"<< outTag<<", Timestep " << TimeStep << " completed." << flush;
    }

    cout << "\rFile:"<< outTag<<", Timestep " << TimeStep << " completed." << endl;
    gettimeofday(&end, NULL);
    cout << "Time taken: " << end.tv_sec - start.tv_sec << "s" << endl;
    }
    return 0;
}