#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <vector>
#include <mpi.h>
#include <random>
#include <sstream>
#include <string>
#include <cassert>

using namespace std;
const double pi = 4.0 * atan(1.0);





//Simulation paramaters (adjustable):

const int nts = 2; // Number of time steps saved in data arrays

const long long int nx = 256; // Grid Dimensions
const long long int ny = 256;
const long long int nz = 1; // Set nz = 1 for 2D.
const long long int nPos = nx * ny * nz;

int radius=nx/16;
int radius2=nx/2-2;

const double dx = 0.25; //Grid Spacings
const double dy = 0.25;
const double dz = 0.25;
const double dt = 0.05; 

const int nt = (nx * dx / (2 * dt)); // nt required for sim to end at light crossing time is nx*dx/(2*dt)

const int seed = 8;

const string outTag = "_2d_2hdm_nx_0_" + to_string(nx) + "_nt" + to_string(nt) + "_seed" + to_string(seed) + "_Z2_bubbles";

const bool calcEnergy = true; // Output Choices
const bool wallDetect = true;
const bool finalOut = false;

const bool makeGif = true;
const int saveFreq = 5;

const string file_path = __FILE__;
const string dir_path = "Data/"; // Data Directory Location (currently wherever code is located)

const int countRate = 20; // Increments for simulation progress status output.

const string ic_type = "random";

const int nb_fields = 8; // Number of fields in simulation





// 2HDM Z_2 Symmetric Potential Set-Up:
 
// Mass and Energy Paramaters (CAN be chosen)
const long double m_h = 125;
const long double V_sm = 246;
const long double m_H = 200;
const long double m_A = 200;
const long double m_H_pm = 200;

// Scaled Mass and Energy Paramaters (NOT to be edited)
const long double M_h = m_h / m_h;
const long double v_sm = V_sm / V_sm;
const long double M_H = m_H / m_h;
const long double M_A = m_A / m_h;
const long double M_H_pm = m_H_pm / m_h;

// Mixing Angle Paramaters
const long double a = atan(1); 
const long double b = atan(1);

const long double s_a = sin(a);
const long double c_a = cos(a);

const long double s_b = sin(b);
const long double c_b = cos(b);
const long double t_b = tan(b);
const long double ct_b = pow(tan(b), -1);

// Dimensionless Potential Paramaters
const long double mu_1_sq = (1 / pow(M_h, 2)) * 0.5 * ((pow(M_h, 2) * pow(c_a, 2)) + (pow(M_H, 2) * pow(s_a, 2)) + ((pow(M_h, 2) - pow(M_H, 2)) * c_a * s_a * t_b)); // mu_sq paramaters scaled by 1/M_h^2 (sets length scale)
const long double mu_2_sq = (1 / pow(M_h, 2)) * 0.5 * ((pow(M_h, 2) * pow(s_a, 2)) + (pow(M_H, 2) * pow(c_a, 2)) + ((pow(M_h, 2) - pow(M_H, 2)) * c_a * s_a * ct_b));

const double lambda_1 = (pow(v_sm, 2) / pow(M_h, 2)) * (pow(M_h, 2) * pow(c_a, 2) + pow(M_H, 2) * pow(s_a, 2)) / (2 * pow(c_b, 2) * pow(v_sm, 2)); // lambda paramaters scaled by v_SM^2/M_h^2 (sets energy scale)
const double lambda_2 = (pow(v_sm, 2) / pow(M_h, 2)) * (pow(M_h, 2) * pow(s_a, 2) + pow(M_H, 2) * pow(c_a, 2)) / (2 * pow(s_b, 2) * pow(v_sm, 2));
const double lambda_3 = (pow(v_sm, 2) / pow(M_h, 2)) * ((pow(M_h, 2) - pow(M_H, 2)) * c_a * s_a + 2 * pow(M_H_pm, 2) * c_b * s_b) / (c_b * s_b * pow(v_sm, 2));
const double l4_m_l5 = (pow(v_sm, 2) / pow(M_h, 2)) * (-2 * pow(M_H_pm, 2)) / (pow(v_sm, 2));
const double l4_p_l5 = (pow(v_sm, 2) / pow(M_h, 2)) * (2 * (pow(M_A, 2) - pow(M_H_pm, 2))) / (pow(v_sm, 2));

//Doublet VEVs
const long double v1 = c_b * v_sm;
const long double v2 = s_b * v_sm;

const long double minval1 = pow(((4 * lambda_2 * mu_1_sq - 2 * (lambda_3 + l4_m_l5) * mu_2_sq)/(4 * lambda_1 * lambda_2 - (lambda_3 +l4_m_l5))),0.5);
const long double minval2 = pow(((4 * lambda_1 * mu_2_sq - 2 * (lambda_3 + l4_m_l5) * mu_1_sq)/(4 * lambda_1 * lambda_2 - (lambda_3 +l4_m_l5))),0.5);



// Damping Paramaters:
const int damped_nt = 0; // Number of time steps for which damping is imposed. Useful for random initial conditions
const double dampFac = 1.0; // Magnitude of damping term, unclear how strong to make this
const int ntHeld = 1; // Hold fields fixed (but effectively continue expansion) for this number of timesteps. Attempting to get the network into the scaling regime. Not sure how useful this is...
const bool expandDamp = false; // If true then the universe expands during the damping regime.

// Expansion Paramaters:
const double alpha = 2; // Factor multiplying hubble damping term for use in PRS algorithm. alpha = #dims has been claimed to give similar dynamics without changing string width. alpha = #dims - 1 is the usual factor.
const double be_ta = 2; // Scale factor^be_ta is the factor that multiplies the potential contribution to the EoMs. Standard is 2, PRS is 0.
const double scaling = 0; // Power law scaling of the scale factor wrt tau. Using conformal time so rad dom is scaling=1 while matter dom is scaling=2. scaling=0 returns a static universe.




void txt_to_vectordouble(vector<vector<double>>& res, const string& pathname) {
    ifstream infile(pathname);
    
    // 检查文件是否成功打开
    if (!infile.is_open()) {
        cerr << "ERROR: Cannot open file: " << pathname << endl;
        return;
    }
    
    string line;
    int line_count = 0;
    
    while (getline(infile, line)) {
        // 跳过空行
        if (line.empty()) {
            continue;
        }
        
        // 跳过可能存在的注释行（以#开头）
        if (line[0] == '#') {
            continue;
        }
        
        vector<double> row;
        stringstream ss(line);
        double value;
        
        // 使用 stringstream 读取每一行的数值
        while (ss >> value) {
            row.push_back(value);
        }
        
        // 确保行不为空再添加到结果中
        if (!row.empty()) {
            res.push_back(row);
            line_count++;
        }
    }
    
    infile.close();
    
    // 输出读取状态信息
    // cout << "Successfully read " << line_count << " lines from " << pathname << endl;
    
    // 检查数据质量
    if (res.empty()) {
        cerr << "WARNING: No data was read from the file!" << endl;
        return;
    }
    
    // 检查数据维度
    size_t expected_columns = 2; // 根据您的数据格式调整
    for (size_t i = 0; i < res.size(); ++i) {
        if (res[i].size() != expected_columns) {
            cerr << "WARNING: Line " << i + 1 << " has " << res[i].size() 
                 << " columns, expected " << expected_columns << endl;
        }
    }
}
    
double fun(int radius,vector<vector<double>> data,int typ) {
    // cout << data[radius][typ] <<endl;

    //cout << radius<<"_"<<typ<< "_"<< data[radius][typ] <<endl;

    
    
    return data[radius][typ];
}







double draw_circle_0(int area, int position, vector<vector<double>> data){
    int sidelen = nx;
    int zval=position/pow(nx,2);
    int pos2=position%(nx*nx);

    int xval=pos2/nx;
    int yval=pos2%nx;

    long double radi= pow((xval-nx/2+0.5),2)+pow((yval-ny/2+0.5),2);
    double val=fun(0.5*sqrt(radi)*dx/0.1 ,data, 0);
    return val;


}





double draw_circle_1(int area, int position, vector<vector<double>> data){
    int sidelen = nx;
    int zval=position/pow(nx,2);
    int pos2=position%(nx*nx);

    int xval=pos2/nx;
    int yval=pos2%nx;

    long double radi= pow((xval-nx/2+0.5),2)+pow((yval-ny/2+0.5),2);
    double val=fun(0.5*sqrt(radi)*dx/0.1 ,data, 1);
    return val;

}





double draw_circle2(int area, int position, int radius){
    int sidelen = nx;
    int zval=position/pow(nx,2);
    int pos2=position%(nx*nx);

    int xval=pos2/nx;
    int yval=pos2%nx;
    int radi= pow((zval-nz/2),2)+pow((xval-nx/2),2)+pow((yval-ny/2),2);
    double val= exp(-(radi/radius));
    //double val=  pow((pow((radi/radius),2)-1),2);

    //cout <<position<<"-"<<radi<<"-" <<xval <<"-"<< yval <<"-"<< zval << endl;

    // if (radi>pow(radius,2))
    // {return 0;}
    // else{return 1;}
    return val;

}


// Begginning of Simulation:
int main(int argc, char** argv) {
    cout << mu_1_sq <<endl;
    cout << mu_2_sq <<endl;
    cout << lambda_1 <<endl;
    cout << lambda_2 <<endl;
    cout << lambda_3 <<endl;
    cout << l4_m_l5 <<endl;
    cout << l4_p_l5 <<endl;
    vector<vector<double>> data;
    txt_to_vectordouble(data, "expression_2D_polar.txt");

    // Initialize MPI

    // Init MPI
    MPI_Init(&argc, &argv);

    // Get the rank and size
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    long long int chunk = nPos / size;
    long long int chunkRem = nPos - size * chunk;

    long long int coreSize;
    if (rank >= chunkRem) { coreSize = chunk; }
    else { coreSize = chunk + 1; }

    // Calculate the position of the start of the chunk in the full array
    long long int coreStart, coreEnd;
    if (rank < chunkRem) { coreStart = rank * (chunk + 1); coreEnd = (rank + 1) * (chunk + 1); }
    else { coreStart = rank * chunk + chunkRem; coreEnd = (rank + 1) * chunk + chunkRem; }


    // Calculate the halo sizes (all data up to the previous x row at start of chunk and all data up to the next x row at the end of the chunk)
    long long int frontHaloSize, backHaloSize, nbrFrontHaloSize, nbrBackHaloSize, remFront, remBack;
    remFront = coreStart % (ny * nz);
    remBack = coreEnd % (ny * nz);
    if (remFront == 0) { // Smallest possible halo size

        frontHaloSize = 2 * ny * nz;
        nbrBackHaloSize = 2 * ny * nz;

    }
    else {

        // The two sum to 3*ny*nz rather than 2*ny*nz. This is inefficient and should be avoided if possible.

        frontHaloSize = 2 * ny * nz + remFront;
        nbrBackHaloSize = 4 * ny * nz - remFront;

    }

    if (remBack == 0) {

        backHaloSize = 2 * ny * nz;
        nbrFrontHaloSize = 2 * ny * nz;

    }
    else {

        backHaloSize = 4 * ny * nz - remBack;
        nbrFrontHaloSize = 2 * ny * nz + remBack;

    }

    // Size the array needs to be to hold the core and the two halos.
    long long int totSize = frontHaloSize + coreSize + backHaloSize;

    // Calculate the position of the start of the local array (including the haloes) in the full array. This quantity wraps around (i.e -ve numbers mean the other side of array)
    long long int dataStart = coreStart - frontHaloSize;
    long long int dataEnd = coreEnd + backHaloSize;

    // Warnings
    if (rank == 0) {

        if (size == 1) { cout << "Warning: Only one processor being used. This code is not designed for only one processor and may not work." << endl; }
        if (chunk < ny * nz) { cout << "Warning: Chunk size is less than the minimum halo size (i.e chunk neighbour data). Code currently assumes this is not the case so it probably won't work." << endl; }

    }

    // Define variables for the fields
    vector<vector<double>> fields(nb_fields, vector<double>(2 * totSize, 0.0)); 
    vector<vector<double>> fields_dx(nb_fields, vector<double>(2 * totSize, 0.0)); 
    vector<vector<double>> fields_dy(nb_fields, vector<double>(2 * totSize, 0.0)); 
    vector<vector<double>> fields_dz(nb_fields, vector<double>(2 * totSize, 0.0)); 
    vector<vector<double>> fields_dt(nb_fields, vector<double>(2 * totSize, 0.0)); // Define and initialize vector of field vectors.
    vector<double> fieldx(nb_fields, 0.0), fieldy(nb_fields, 0.0), fieldz(nb_fields, 0.0), fieldt(nb_fields, 0.0), fieldtt(nb_fields, 0.0), localKinEnergy(nb_fields, 0.0); // Needed for calculation of energy. nb_fields components.
    double fieldxx, fieldyy, fieldzz; // Only need them to calculate second time derivative of each field individually, thus can be reused.
    double totalLocalEnergy, localNDW, localADW_simple, localADW_full, x0, y0, z0;
    long long int i, j, k, TimeStep, tNow, tPast, comp, imx, ipx, imy, ipy, imz, ipz, ipxmy, ipxmz, imxpy, ipymz, imxpz, imypz, imxx, ipxx, imyy, ipyy, imzz, ipzz;





    struct timeval start, end;
    if (rank == 0) { gettimeofday(&start, NULL); }

    stringstream ss;

    MPI_Barrier(MPI_COMM_WORLD);
    

    //Creates Output Files if required
    string icPath = dir_path + "ic.txt";
    ifstream ic(icPath.c_str());


    string finalFieldPath = dir_path + "finalFields" + outTag + ".txt";
    ofstream finalFields(finalFieldPath.c_str());


    string valsPerLoopPath = dir_path + "valsPerLoop" + outTag + ".txt";
    ofstream valsPerLoop(valsPerLoopPath.c_str());



    // Index values (not neccessarily on grid and hence not integers) of the zero coordinate.
    x0 = 0.5 * (nx - 1);
    y0 = 0.5 * (ny - 1);
    z0 = 0.5 * (nz - 1);



    if (ic_type == "random") {

        // Creates and assigns RIC for each of the 8 fields independantly.

        // Use the seed to generate the data
        mt19937 generator_1(seed);
        mt19937 generator_2(seed + 1);
        mt19937 generator_3(seed + 2);
        mt19937 generator_4(seed + 3);
        mt19937 generator_5(seed + 4);
        mt19937 generator_6(seed + 5);
        mt19937 generator_7(seed + 6);
        mt19937 generator_8(seed + 7);

        uniform_real_distribution<double> distribution(-0.0, 0.0); // Uniform distribution for the phase of the strings


        double phi1Assign;
        double phi2Assign;
        double phi3Assign;
        double phi4Assign;
        double phi5Assign;
        double phi6Assign;
        double phi7Assign;
        double phi8Assign;



        // Skip the random numbers ahead to the appropriate point.
        for (i = 0; i < coreStart; i++) {
            phi1Assign = distribution(generator_1);
            phi2Assign = distribution(generator_2);
            phi3Assign = distribution(generator_3);
            phi4Assign = distribution(generator_4);
            phi5Assign = distribution(generator_5);
            phi6Assign = distribution(generator_6);
            phi7Assign = distribution(generator_7);
            phi8Assign = distribution(generator_8);
        }



        for (i = frontHaloSize; i < coreSize + frontHaloSize; i++) {
        //cout <<"init-"<<i<<"-"<<i+rank*chunk<<"-"<<rank<<"-"<<size<< endl;
            phi1Assign = distribution(generator_1);
            phi2Assign = distribution(generator_2);
            phi3Assign = distribution(generator_3);
            phi4Assign = distribution(generator_4);
            phi5Assign = distribution(generator_5);
            phi6Assign = distribution(generator_6);
            phi7Assign = distribution(generator_7);
            phi8Assign = distribution(generator_8);

            //Assign values to fields
            fields[0][i] = phi1Assign;
            fields[1][i] = phi2Assign;
            fields[2][i] = phi3Assign;
            fields[3][i] = phi4Assign;
            fields[4][i] = phi5Assign;
            fields[5][i] = phi6Assign;
            fields[6][i] = phi7Assign;
            fields[7][i] = phi8Assign;




            fields[0][i] =fields[0][i]+ draw_circle_0( nPos,  i+rank*chunk-frontHaloSize,  data)+ draw_circle_1( nPos,  i+rank*chunk-frontHaloSize,  data);
            fields[2][i] =fields[2][i]- draw_circle_0( nPos,  i+rank*chunk-frontHaloSize,  data);
            fields[4][i] =fields[4][i]- draw_circle_0( nPos,  i+rank*chunk-frontHaloSize,  data)+ draw_circle_1( nPos,  i+rank*chunk-frontHaloSize,  data);
            fields[6][i] =fields[6][i]+ draw_circle_0( nPos,  i+rank*chunk-frontHaloSize,  data);


            fields[1][i] =0;
            fields[3][i] =0;
            fields[5][i] =0;
            fields[7][i] =0;






            fields[0][totSize + i] = + draw_circle_0( nPos,  i+rank*chunk-frontHaloSize,  data) + draw_circle_1( nPos,  i+rank*chunk-frontHaloSize,  data)*cos(0.6*dt);
            fields[1][totSize + i] = 0;
            fields[2][totSize + i] = - draw_circle_0( nPos,  i+rank*chunk-frontHaloSize,  data) - draw_circle_1( nPos,  i+rank*chunk-frontHaloSize,  data)*sin(0.6*dt);
            fields[3][totSize + i] = 0;
            fields[4][totSize + i] = - draw_circle_0( nPos,  i+rank*chunk-frontHaloSize,  data) + draw_circle_1( nPos,  i+rank*chunk-frontHaloSize,  data)*cos(0.6*dt);
            fields[5][totSize + i] = 0;
            fields[6][totSize + i] = + draw_circle_0( nPos,  i+rank*chunk-frontHaloSize,  data) - draw_circle_1( nPos,  i+rank*chunk-frontHaloSize,  data)*sin(0.6*dt);
            fields[7][totSize + i] = 0;
            
            // cout << "." << endl; 
                }
        

        // Now that the core data has been generated, need to communicate the haloes between processes:

        // Loop over the different fields (the nb_fields components of the vector of fields)
        for (comp = 0; comp < nb_fields; comp++) {

            MPI_Sendrecv(&fields[comp][frontHaloSize], nbrBackHaloSize, MPI_DOUBLE, (rank - 1 + size) % size, comp, // Send this
                &fields[comp][coreSize + frontHaloSize], backHaloSize, MPI_DOUBLE, (rank + 1) % size, comp, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Receive this

            MPI_Sendrecv(&fields[comp][coreSize + frontHaloSize - nbrFrontHaloSize], nbrFrontHaloSize, MPI_DOUBLE, (rank + 1) % size, comp,
                &fields[comp][0], frontHaloSize, MPI_DOUBLE, (rank - 1 + size) % size, comp, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Sendrecv(&fields[comp][totSize + frontHaloSize], nbrBackHaloSize, MPI_DOUBLE, (rank - 1 + size) % size, comp + nb_fields,
                &fields[comp][totSize + coreSize + frontHaloSize], backHaloSize, MPI_DOUBLE, (rank + 1) % size, comp + nb_fields, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Sendrecv(&fields[comp][totSize + coreSize + frontHaloSize - nbrFrontHaloSize], nbrFrontHaloSize, MPI_DOUBLE, (rank + 1) % size, comp + nb_fields,
                &fields[comp][totSize], frontHaloSize, MPI_DOUBLE, (rank - 1 + size) % size, comp + nb_fields, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


        }


    }




    gettimeofday(&end, NULL);

    if (rank == 0) { cout << "Initial data loaded/generated in: " << end.tv_sec - start.tv_sec << "s" << endl; }




    // Main for loop that evolves the fields:
    for (TimeStep = 0; TimeStep < nt; TimeStep++) {
        
            //cout << TimeStep << endl;

        double fric, tau;


        // Expansion during damping check:
        if (expandDamp) { tau = (1 + (ntHeld + TimeStep) * dt)/(1 + (ntHeld) * dt); }
        else { tau = (1 + (ntHeld + TimeStep - damped_nt) * dt)/(1 + (ntHeld - damped_nt) * dt); }


        // Is damping switched on or not?
        if (TimeStep < damped_nt) {

            if (expandDamp) { fric = dampFac + alpha * scaling / tau; } // denominator is conformal time
            else { fric = dampFac; }

        }
        else {

            if (expandDamp) { fric = alpha * scaling / tau; } // Time needs to have moved along during the damped phase
            else { fric = alpha * scaling / tau; } // Time was not progressing during the damped phase

        }

        tNow = (TimeStep + 1) % 2;
        tPast = TimeStep % 2;





        // Main calculations and evolutions section, using the EoM:
        totalLocalEnergy = 0;
        for (comp = 0; comp < nb_fields; comp++) { localKinEnergy[comp] = 0; }
        localNDW = 0;
        localADW_simple = 0;
        localADW_full = 0;

        //Loops over and evolves all assigned data points
        for (i = frontHaloSize; i < coreSize + frontHaloSize; i++) { 
            
            
            //Encoding the periodic boundary conditions for the evolving data point:
            
            // No need to worry about periodicity with the x neighbours because halo is designed to contain them
            imx = i - ny * nz;
            ipx = i + ny * nz;
            imxx = i - 2 * ny * nz;
            ipxx = i + 2 * ny * nz;


            // Need to account for the periodicity of the space for the other two directions:
            // Convert to global position in array to do modulo arithmetic. The second to last term gives ny*nz*floor((i+dataStart)/(ny*nz)). The last term converts back to the position in the local array 
            imy = (i + dataStart - nz + ny * nz) % (ny * nz) + ((i + dataStart) / (ny * nz)) * ny * nz - dataStart;
            ipy = (i + dataStart + nz) % (ny * nz) + ((i + dataStart) / (ny * nz)) * ny * nz - dataStart;
            imyy = (i + dataStart - 2 * nz + ny * nz) % (ny * nz) + ((i + dataStart) / (ny * nz)) * ny * nz - dataStart;
            ipyy = (i + dataStart + 2 * nz) % (ny * nz) + ((i + dataStart) / (ny * nz)) * ny * nz - dataStart;

            imz = (i + dataStart - 1 + nz) % nz + ((i + dataStart) / nz) * nz - dataStart;
            ipz = (i + dataStart + 1) % nz + ((i + dataStart) / nz) * nz - dataStart;
            imzz = (i + dataStart - 2 + nz) % nz + ((i + dataStart) / nz) * nz - dataStart;
            ipzz = (i + dataStart + 2) % nz + ((i + dataStart) / nz) * nz - dataStart;

            // Additionally needed for wilson loop calculations. Avoid using x shifted points first as this makes the calculations more complicated and some of these points aren't in the correct positions
            ipxmy = imy + ny * nz;
            ipxmz = imz + ny * nz;
            imxpy = ipy - ny * nz;
            ipymz = (ipy + dataStart - 1 + nz) % nz + ((ipy + dataStart) / nz) * nz - dataStart;
            imxpz = ipz - ny * nz;
            imypz = (imy + dataStart + 1) % nz + ((imy + dataStart) / nz) * nz - dataStart;





            // Define in Shorthand the current value of each field at the evolving point:
            double f0 = fields[0][totSize * tNow + i];
            double f1 = fields[1][totSize * tNow + i];
            double f2 = fields[2][totSize * tNow + i];
            double f3 = fields[3][totSize * tNow + i];
            double f4 = fields[4][totSize * tNow + i];
            double f5 = fields[5][totSize * tNow + i];
            double f6 = fields[6][totSize * tNow + i];
            double f7 = fields[7][totSize * tNow + i];





            // Calculate the second order time derivatives and update the fields:


            //First Doublet:
            for (comp = 0; comp <= 3; comp++) {

                //2nd Spatial Deriviitive Calculations (central finite difference method) - to 4th Order:
                fieldxx = (16 * (fields[comp][totSize * tNow + ipx] + fields[comp][totSize * tNow + imx]) - 30 * fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + ipxx] - fields[comp][totSize * tNow + imxx]) / (12 * dx * dx * a * a);
                fieldyy = (16 * (fields[comp][totSize * tNow + ipy] + fields[comp][totSize * tNow + imy]) - 30 * fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + ipyy] - fields[comp][totSize * tNow + imyy]) / (12 * dy * dy * a * a);
                fieldzz = (16 * (fields[comp][totSize * tNow + ipz] + fields[comp][totSize * tNow + imz]) - 30 * fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + ipzz] - fields[comp][totSize * tNow + imzz]) / (12 * dz * dz * a * a);

                //1st Temporal Derivitive Calculation (only required for damping):
                fieldt[comp] = (fields[comp][totSize * tNow + i] - fields[comp][totSize * tPast + i]) / dt;
                fieldx[comp] = (fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + imx]) / (dx * a);
                fieldy[comp] = (fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + imy]) / (dy * a);
                fieldz[comp] = (fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + imz]) / (dz * a);
                


                //2nd Temporal Derivitive Calculation - Using EoM:
                fieldtt[comp] = fieldxx + fieldyy + fieldzz - fric * fieldt[comp] - pow(pow(tau, scaling), be_ta) * (-mu_1_sq * fields[comp][totSize * tNow + i] + 2 * lambda_1 * fields[comp][totSize * tNow + i] * (pow(f0, 2) + pow(f1, 2) + pow(f2, 2) + pow(f3, 2)) + lambda_3 * fields[comp][totSize * tNow + i] * (pow(f4, 2) + pow(f5, 2) + pow(f6, 2) + pow(f7, 2)));

                //Different potential contributions dependant upon the component of the field vector:
                if (comp == 0) {

                    fieldtt[comp] += -pow(pow(tau, scaling), be_ta) * ((l4_m_l5)*f4 * (f0 * f4 + f1 * f5 + f2 * f6 + f3 * f7));
                    fieldtt[comp] += -pow(pow(tau, scaling), be_ta) * ((l4_p_l5)*f5 * (f0 * f5 - f1 * f4 + f2 * f7 - f3 * f6));

                }

                if (comp == 1) {
                    fieldtt[comp] =0;

                }

                if (comp == 2) {

                    fieldtt[comp] += -pow(pow(tau, scaling), be_ta) * ((l4_m_l5)*f6 * (f0 * f4 + f1 * f5 + f2 * f6 + f3 * f7));
                    fieldtt[comp] += -pow(pow(tau, scaling), be_ta) * ((l4_p_l5)*f7 * (f0 * f5 - f1 * f4 + f2 * f7 - f3 * f6));

                }

                if (comp == 3) {

                    fieldtt[comp] =0;

                }

                // Updates value of the field (2nd Order centralfinite difference method):
                double val = 2 * fields[comp][totSize * tNow + i] - fields[comp][totSize * tPast + i] + dt * dt * fieldtt[comp];

                double xval1 = (i+rank*chunk-frontHaloSize)/nx-nx/2+0.5;
                double yval1 = (i+rank*chunk-frontHaloSize)%nx-ny/2+0.5;
                
                // cout << xval1<< "_" << yval1<<  "_" << radius<< endl;
                if (pow(xval1,2)+pow(yval1,2)<=pow(radius2,2)) {fields[comp][totSize * tPast + i] = val;}
                // fields[comp][totSize * tPast + i] = val;
                if (pow(xval1,2)+pow(yval1,2)>pow(radius2,2) and xval1>0 and yval1>0) {fields[comp][totSize * tPast + i] = (fields[comp][totSize * tNow + imx]*xval1*xval1 +fields[comp][totSize * tNow + imy]*yval1*yval1 )/((pow(xval1,2)+pow(yval1,2)));}
                if (pow(xval1,2)+pow(yval1,2)>pow(radius2,2) and xval1>0 and yval1<0) {fields[comp][totSize * tPast + i] = (fields[comp][totSize * tNow + imx]*xval1*xval1 +fields[comp][totSize * tNow + ipy]*yval1*yval1 )/((pow(xval1,2)+pow(yval1,2)));}
                if (pow(xval1,2)+pow(yval1,2)>pow(radius2,2) and xval1<0 and yval1>0) {fields[comp][totSize * tPast + i] = (fields[comp][totSize * tNow + ipx]*xval1*xval1 +fields[comp][totSize * tNow + imy]*yval1*yval1 )/((pow(xval1,2)+pow(yval1,2)));}
                if (pow(xval1,2)+pow(yval1,2)>pow(radius2,2) and xval1<0 and yval1<0) {fields[comp][totSize * tPast + i] = (fields[comp][totSize * tNow + ipx]*xval1*xval1 +fields[comp][totSize * tNow + ipy]*yval1*yval1 )/((pow(xval1,2)+pow(yval1,2)));}
                fields_dt[comp][totSize * tPast + i] = fieldt[comp];
                fields_dx[comp][totSize * tPast + i] = fieldx[comp];
                fields_dy[comp][totSize * tPast + i] = fieldy[comp];
                fields_dz[comp][totSize * tPast + i] = fieldz[comp];

            }

            //Second Doublet:
            for (comp = 4; comp <= 7; comp++) {

                //2nd Spatial Deriviitive Calculations (central finite difference method) - to 4th Order:
                fieldxx = (16 * (fields[comp][totSize * tNow + ipx] + fields[comp][totSize * tNow + imx]) - 30 * fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + ipxx] - fields[comp][totSize * tNow + imxx]) / (12 * dx * dx * a * a);
                fieldyy = (16 * (fields[comp][totSize * tNow + ipy] + fields[comp][totSize * tNow + imy]) - 30 * fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + ipyy] - fields[comp][totSize * tNow + imyy]) / (12 * dy * dy * a * a);
                fieldzz = (16 * (fields[comp][totSize * tNow + ipz] + fields[comp][totSize * tNow + imz]) - 30 * fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + ipzz] - fields[comp][totSize * tNow + imzz]) / (12 * dz * dz * a * a);

                //1st Temporal Derivitive Calculation (only required for damping):
                fieldt[comp] = (fields[comp][totSize * tNow + i] - fields[comp][totSize * tPast + i]) / dt;
                fieldx[comp] = (fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + imx]) / (dx * a);
                fieldy[comp] = (fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + imy]) / (dy * a);
                fieldz[comp] = (fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + imz]) / (dz * a);

                //2nd Temporal Derivitive Calculation - Using EoM:
                fieldtt[comp] = fieldxx + fieldyy + fieldzz - fric * fieldt[comp] - pow(pow(tau, scaling), be_ta) * (-mu_2_sq * fields[comp][totSize * tNow + i] + lambda_3 * fields[comp][totSize * tNow + i] * (pow(f0, 2) + pow(f1, 2) + pow(f2, 2) + pow(f3, 2)) + 2 * lambda_2 * fields[comp][totSize * tNow + i] * (pow(f4, 2) + pow(f5, 2) + pow(f6, 2) + pow(f7, 2)));

                //Different potential contributions dependant upon the component of the field vector:
                if (comp == 4) {

                    fieldtt[comp] += -pow(pow(tau, scaling), be_ta) * ((l4_m_l5)*f0 * (f0 * f4 + f1 * f5 + f2 * f6 + f3 * f7));
                    fieldtt[comp] += -pow(pow(tau, scaling), be_ta) * ((l4_p_l5) * -f1 * (f0 * f5 - f1 * f4 + f2 * f7 - f3 * f6));

                }

                if (comp == 5) {

                    fieldtt[comp] =0;

                }

                if (comp == 6) {

                    fieldtt[comp] += -pow(pow(tau, scaling), be_ta) * ((l4_m_l5)*f2 * (f0 * f4 + f1 * f5 + f2 * f6 + f3 * f7));
                    fieldtt[comp] += -pow(pow(tau, scaling), be_ta) * ((l4_p_l5) * -f3 * (f0 * f5 - f1 * f4 + f2 * f7 - f3 * f6));

                }

                if (comp == 7) {

                    fieldtt[comp] =0;
                }

                // Updates value of the field (2nd Order centralfinite difference method):
                double val = 2 * fields[comp][totSize * tNow + i] - fields[comp][totSize * tPast + i] + dt * dt * fieldtt[comp];

                double xval1 = (i+rank*chunk-frontHaloSize)/nx-nx/2+0.5;
                double yval1 = (i+rank*chunk-frontHaloSize)%nx-ny/2+0.5;
                
                // cout << xval1<< "_" << yval1<<  "_" << radius<< endl;
                if (pow(xval1,2)+pow(yval1,2)<=pow(radius2,2)) {fields[comp][totSize * tPast + i] = val;}
                // fields[comp][totSize * tPast + i] = val;
                if (pow(xval1,2)+pow(yval1,2)>pow(radius2,2) and xval1>0 and yval1>0) {fields[comp][totSize * tPast + i] = (fields[comp][totSize * tNow + imx]*xval1*xval1 +fields[comp][totSize * tNow + imy]*yval1*yval1 )/((pow(xval1,2)+pow(yval1,2)));}
                if (pow(xval1,2)+pow(yval1,2)>pow(radius2,2) and xval1>0 and yval1<0) {fields[comp][totSize * tPast + i] = (fields[comp][totSize * tNow + imx]*xval1*xval1 +fields[comp][totSize * tNow + ipy]*yval1*yval1 )/((pow(xval1,2)+pow(yval1,2)));}
                if (pow(xval1,2)+pow(yval1,2)>pow(radius2,2) and xval1<0 and yval1>0) {fields[comp][totSize * tPast + i] = (fields[comp][totSize * tNow + ipx]*xval1*xval1 +fields[comp][totSize * tNow + imy]*yval1*yval1 )/((pow(xval1,2)+pow(yval1,2)));}
                if (pow(xval1,2)+pow(yval1,2)>pow(radius2,2) and xval1<0 and yval1<0) {fields[comp][totSize * tPast + i] = (fields[comp][totSize * tNow + ipx]*xval1*xval1 +fields[comp][totSize * tNow + ipy]*yval1*yval1 )/((pow(xval1,2)+pow(yval1,2)));}
                fields_dt[comp][totSize * tPast + i] = fieldt[comp];
                fields_dx[comp][totSize * tPast + i] = fieldx[comp];
                fields_dy[comp][totSize * tPast + i] = fieldy[comp];
                fields_dz[comp][totSize * tPast + i] = fieldz[comp];

            }





            // Calculate the energy contained in this process's domain
            if (calcEnergy) {

                for (comp = 0; comp < nb_fields; comp++) {
                    fieldx[comp] = (fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + imx]) / dx;
                    fieldy[comp] = (fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + imy]) / dy;
                    fieldz[comp] = (fields[comp][totSize * tNow + i] - fields[comp][totSize * tNow + imz]) / dz;
                    // cout <<fieldz[comp]<< endl;

                    localKinEnergy[comp] += (pow(fieldt[comp], 2) + pow(fieldx[comp], 2) + pow(fieldy[comp], 2) + pow(fieldz[comp], 2)) * dx * dy * dz; // Kinetic energy

                    totalLocalEnergy += localKinEnergy[comp];
                }

                //Potential Terms
                totalLocalEnergy += -mu_1_sq * (pow(f0, 2) + pow(f1, 2) + pow(f2, 2) + pow(f3, 2)) - mu_2_sq * (pow(f4, 2) + pow(f5, 2) + pow(f6, 2) + pow(f7, 2));
                totalLocalEnergy += lambda_1 * pow((pow(f0, 2) + pow(f1, 2) + pow(f2, 2) + pow(f3, 2)), 2);
                totalLocalEnergy += lambda_2 * pow((pow(f4, 2) + pow(f5, 2) + pow(f6, 2) + pow(f7, 2)), 2);
                totalLocalEnergy += lambda_3 * ((pow(f0, 2) + pow(f1, 2) + pow(f2, 2) + pow(f3, 2)) * (pow(f4, 2) + pow(f5, 2) + pow(f6, 2) + pow(f7, 2)));
                totalLocalEnergy += (l4_m_l5)*pow((f0 * f4 + f1 * f5 + f2 * f6 + f3 * f7), 2);
                totalLocalEnergy += (l4_p_l5)*pow((f0 * f5 - f1 * f4 + f2 * f7 - f3 * f6), 2);

            }

            // If the sign of phi flips between any two neighbours, consider that as a wall detection. Sum this up. Only look at forward neighbours so I'm not double counting.
            if (wallDetect) {

                double R1_i = 2 * (fields[0][totSize * tNow + i] * fields[4][totSize * tNow + i] + fields[1][totSize * tNow + i] * fields[5][totSize * tNow + i] + fields[2][totSize * tNow + i] * fields[6][totSize * tNow + i] + fields[3][totSize * tNow + i] * fields[7][totSize * tNow + i]);
                double R1_ipx = 2 * (fields[0][totSize * tNow + ipx] * fields[4][totSize * tNow + ipx] + fields[1][totSize * tNow + ipx] * fields[5][totSize * tNow + ipx] + fields[2][totSize * tNow + ipx] * fields[6][totSize * tNow + ipx] + fields[3][totSize * tNow + ipx] * fields[7][totSize * tNow + ipx]);
                double R1_ipy = 2 * (fields[0][totSize * tNow + ipy] * fields[4][totSize * tNow + ipy] + fields[1][totSize * tNow + ipy] * fields[5][totSize * tNow + ipy] + fields[2][totSize * tNow + ipy] * fields[6][totSize * tNow + ipy] + fields[3][totSize * tNow + ipy] * fields[7][totSize * tNow + ipy]);
                double R1_ipz = 2 * (fields[0][totSize * tNow + ipz] * fields[4][totSize * tNow + ipz] + fields[1][totSize * tNow + ipz] * fields[5][totSize * tNow + ipz] + fields[2][totSize * tNow + ipz] * fields[6][totSize * tNow + ipz] + fields[3][totSize * tNow + ipz] * fields[7][totSize * tNow + ipz]);

                double R1_imx = 2 * (fields[0][totSize * tNow + imx] * fields[4][totSize * tNow + imx] + fields[1][totSize * tNow + imx] * fields[5][totSize * tNow + imx] + fields[2][totSize * tNow + imx] * fields[6][totSize * tNow + imx] + fields[3][totSize * tNow + imx] * fields[7][totSize * tNow + imx]);
                double R1_imy = 2 * (fields[0][totSize * tNow + imy] * fields[4][totSize * tNow + imy] + fields[1][totSize * tNow + imy] * fields[5][totSize * tNow + imy] + fields[2][totSize * tNow + imy] * fields[6][totSize * tNow + imy] + fields[3][totSize * tNow + imy] * fields[7][totSize * tNow + imy]);
                double R1_imz = 2 * (fields[0][totSize * tNow + imz] * fields[4][totSize * tNow + imz] + fields[1][totSize * tNow + imz] * fields[5][totSize * tNow + imz] + fields[2][totSize * tNow + imz] * fields[6][totSize * tNow + imz] + fields[3][totSize * tNow + imz] * fields[7][totSize * tNow + imz]);



                // x neighbour
                if (R1_i * R1_ipx < 0) {

                    localNDW += 1;
                    localADW_simple += 2.0 * dy * dz / 3.0;

                    double R1x = (R1_i - R1_imx) / dx;
                    double R1y = (R1_i - R1_imy) / dy;
                    double R1z = (R1_i - R1_imz) / dz;
                    localADW_full += dy * dz * sqrt(pow(R1x, 2) + pow(R1y, 2) + pow(R1z, 2)) / (abs(R1x) + abs(R1y) + abs(R1z));
                }

                // y neighbour
                if (R1_i * R1_ipy < 0) {

                    localNDW += 1;
                    localADW_simple += 2.0 * dy * dz / 3.0;

                    double R1x = (R1_i - R1_imx) / dx;
                    double R1y = (R1_i - R1_imy) / dy;
                    double R1z = (R1_i - R1_imz) / dz;
                    localADW_full += dx * dz * sqrt(pow(R1x, 2) + pow(R1y, 2) + pow(R1z, 2)) / (abs(R1x) + abs(R1y) + abs(R1z));
                }

                // z neighbour
                if (R1_i * R1_ipz < 0) {

                    localNDW += 1;
                    localADW_simple += 2.0 * dy * dz / 3.0;

                    double R1x = (R1_i - R1_imx) / dx;
                    double R1y = (R1_i - R1_imy) / dy;
                    double R1z = (R1_i - R1_imz) / dz;
                    localADW_full += dx * dy * sqrt(pow(R1x, 2) + pow(R1y, 2) + pow(R1z, 2)) / (abs(R1x) + abs(R1y) + abs(R1z));
                }

            }

        }

        // Puts required headers on valsPerLoop output file:
        if (TimeStep == 0 and rank == 0) {
            if (calcEnergy and wallDetect) { valsPerLoop << "Energy" << " " << "NDW" << " " << "ADW_Simple" << " " << "ADW_Full" << endl; }
            else {
                if (calcEnergy) { valsPerLoop << "Energy" << endl; }
                if (wallDetect) { valsPerLoop << "NDW" << " " << "ADW_Simple" << " " << "ADW_Full" << endl; }
            }
        }



        // If calculating the energy, add it all up and output to text
        if (calcEnergy) {

            if (rank == 0) {

                double energy = totalLocalEnergy; // Initialise the energy as the energy in the domain of this process. Then add the energy in the regions of the other processes.

                for (i = 1; i < size; i++) { MPI_Recv(&totalLocalEnergy, 1, MPI_DOUBLE, i, 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  energy += totalLocalEnergy; }

                valsPerLoop << energy << " ";

            }
            else { MPI_Send(&totalLocalEnergy, 1, MPI_DOUBLE, 0, 20, MPI_COMM_WORLD); }

        }


        // Sum up the locally detected walls and output to text
        if (wallDetect) {

            if (rank == 0) {

                double NDW = localNDW;
                double ADW_simple = localADW_simple;
                double ADW_full = localADW_full;

                for (i = 1; i < size; i++) {

                    MPI_Recv(&localNDW, 1, MPI_DOUBLE, i, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    NDW += localNDW;

                    MPI_Recv(&localADW_simple, 1, MPI_DOUBLE, i, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    ADW_simple += localADW_simple;

                    MPI_Recv(&localADW_full, 1, MPI_DOUBLE, i, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    ADW_full += localADW_full;

                }

                valsPerLoop << NDW << " " << ADW_simple << " " << ADW_full;

            }
            else {

                MPI_Send(&localNDW, 1, MPI_DOUBLE, 0, 21, MPI_COMM_WORLD);
                MPI_Send(&localADW_simple, 1, MPI_DOUBLE, 0, 21, MPI_COMM_WORLD);
                MPI_Send(&localADW_full, 1, MPI_DOUBLE, 0, 21, MPI_COMM_WORLD);

            }

        }


        if (rank == 0 and (calcEnergy or wallDetect)) { valsPerLoop << endl; }


        // Update the core
        // Send sections of the core that are haloes for the other processes across to the relevant process. Then receive data for the halo of this process.
        for (comp = 0; comp < nb_fields; comp++) {

            MPI_Sendrecv(&fields[comp][totSize * tPast + frontHaloSize], nbrBackHaloSize, MPI_DOUBLE, (rank - 1 + size) % size, comp, // Send this
                &fields[comp][totSize * tPast + coreSize + frontHaloSize], backHaloSize, MPI_DOUBLE, (rank + 1) % size, comp, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Receive this

            MPI_Sendrecv(&fields[comp][totSize * tPast + coreSize + frontHaloSize - nbrFrontHaloSize], nbrFrontHaloSize, MPI_DOUBLE, (rank + 1) % size, comp,
                &fields[comp][totSize * tPast], frontHaloSize, MPI_DOUBLE, (rank - 1 + size) % size, comp, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }





        //Output the final fields.
        if (finalOut and TimeStep == nt - 1) {

            if (rank == 0) {
                finalFields << "R0" << " " << "R1" << " " << "R2" << " " << "R3" << " " << "R4" << " " << "R5" << " " << "n1" << " " << "n2" << " " << "n3" << " " << "psi1" << " " << "psi2" << " " << "psi3" << " " << "psi4" << " " << "psi5" << " " << "psi6" << " " << "psi7" << " " << "psi8" << " " << "Jemx" << " " <<"Jemy" << " " << "Jemz" << " " << "E" <<  " " << "KE" <<  " " << "PE" << endl;

                vector<vector<double>> fieldsOut(nb_fields, vector<double>(nPos, 0.0));
                vector<vector<double>> fieldsOut_dx(nb_fields, vector<double>(nPos, 0.0));
                vector<vector<double>> fieldsOut_dy(nb_fields, vector<double>(nPos, 0.0));
                vector<vector<double>> fieldsOut_dz(nb_fields, vector<double>(nPos, 0.0));
                vector<vector<double>> fieldsOut_dt(nb_fields, vector<double>(nPos, 0.0));

                double R0;
                double R1;
                double R2;
                double R3;
                double R4;
                double R5;

                double n1;
                double n2;
                double n3;

                double p1;
                double p2;
                double p3;
                double p4;
                double p5;
                double p6;
                double p7;
                double p8;

                double p1x;
                double p2x;
                double p3x;
                double p4x;
                double p5x;
                double p6x;
                double p7x;
                double p8x;

                double p1y;
                double p2y;
                double p3y;
                double p4y;
                double p5y;
                double p6y;
                double p7y;
                double p8y;

                double p1z;
                double p2z;
                double p3z;
                double p4z;
                double p5z;
                double p6z;
                double p7z;
                double p8z;

                double J0x;
                double J0y;
                double J0z;
                
                double J1x;
                double J1y;
                double J1z;
                
                double J2x;
                double J2y;
                double J2z;
                
                double J3x;
                double J3y;
                double J3z;
                
                double Jemx;
                double Jemy;
                double Jemz;

                double energy;
                double ke;
                double pe;
                vector<double> energyv(nPos, 0.0);




                for (comp = 0; comp < nb_fields; comp++) {

                    for (i = 0; i < coreSize; i++) {
                        fieldsOut[comp][i] = fields[comp][frontHaloSize + i]; 



                        fieldsOut_dx[comp][i] = fields_dx[comp][frontHaloSize + i]; 
                        fieldsOut_dy[comp][i] = fields_dy[comp][frontHaloSize + i]; 
                        fieldsOut_dz[comp][i] = fields_dz[comp][frontHaloSize + i]; 
                        fieldsOut_dt[comp][i] = fields_dt[comp][frontHaloSize + i]; 
                        energyv[i] =+ (pow(fieldsOut_dx[comp][i], 2) + pow(fieldsOut_dy[comp][i], 2) + pow(fieldsOut_dz[comp][i], 2) + pow(fieldsOut_dt[comp][i], 2)) * dx * dy * dz; // Kinetic energy

                    }


                    for (i = 1; i < size; i++) {

                        int localCoreStart;
                        int localCoreSize;
                        if (i < chunkRem) { localCoreStart = i * (chunk + 1); localCoreSize = chunk + 1; }
                        else { localCoreStart = i * chunk + chunkRem; localCoreSize = chunk; }

                        MPI_Recv(&fieldsOut[comp][localCoreStart], localCoreSize, MPI_DOUBLE, i, comp, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    }

                }

                //cout <<fieldsOut[0][0]<<"-"<<fieldsOut[1][0]<<"-" <<fieldsOut[2][0]<<"-"<< fieldsOut[3][0]<<"-"<< fieldsOut[4][0]<<"-"<< fieldsOut[5][0]<<"-"<< fieldsOut[6][0]<<"-"<< fieldsOut[7][0]<< endl;

                for (i = 0; i < nPos; i++) {

                    R0 = pow(fieldsOut[0][i], 2) + pow(fieldsOut[1][i], 2) + pow(fieldsOut[2][i], 2) + pow(fieldsOut[3][i], 2) + pow(fieldsOut[4][i], 2) + pow(fieldsOut[5][i], 2) + pow(fieldsOut[6][i], 2) + pow(fieldsOut[7][i], 2);
                    R1 = 2 * (fieldsOut[0][i] * fieldsOut[4][i] + fieldsOut[1][i] * fieldsOut[5][i] + fieldsOut[2][i] * fieldsOut[6][i] + fieldsOut[3][i] * fieldsOut[7][i]);
                    R2 = 2 * (fieldsOut[0][i] * fieldsOut[5][i] - fieldsOut[1][i] * fieldsOut[4][i] + fieldsOut[2][i] * fieldsOut[7][i] - fieldsOut[3][i] * fieldsOut[6][i]);
                    R3 = pow(fieldsOut[0][i], 2) + pow(fieldsOut[1][i], 2) + pow(fieldsOut[2][i], 2) + pow(fieldsOut[3][i], 2) - pow(fieldsOut[4][i], 2) - pow(fieldsOut[5][i], 2) - pow(fieldsOut[6][i], 2) - pow(fieldsOut[7][i], 2);
                    R4 = 2 * (fieldsOut[0][i] * fieldsOut[6][i] - fieldsOut[1][i] * fieldsOut[7][i] - fieldsOut[2][i] * fieldsOut[4][i] + fieldsOut[3][i] * fieldsOut[5][i]);
                    R5 = 2 * (fieldsOut[0][i] * fieldsOut[7][i] + fieldsOut[1][i] * fieldsOut[6][i] - fieldsOut[2][i] * fieldsOut[5][i] - fieldsOut[3][i] * fieldsOut[4][i]);

                    n1 = -2 * (fieldsOut[0][i] * fieldsOut[2][i] + fieldsOut[1][i] * fieldsOut[3][i] + fieldsOut[4][i] * fieldsOut[6][i] + fieldsOut[5][i] * fieldsOut[7][i]);
                    n2 = -2 * (fieldsOut[0][i] * fieldsOut[3][i] - fieldsOut[1][i] * fieldsOut[2][i] + fieldsOut[4][i] * fieldsOut[7][i] - fieldsOut[5][i] * fieldsOut[6][i]);
                    n3 = -1 * (pow(fieldsOut[0][i], 2) + pow(fieldsOut[1][i], 2) - pow(fieldsOut[2][i], 2) - pow(fieldsOut[3][i], 2) + pow(fieldsOut[4][i], 2) + pow(fieldsOut[5][i], 2) - pow(fieldsOut[6][i], 2) - pow(fieldsOut[7][i], 2));

                    p1 = fieldsOut[0][i];
                    p2 = fieldsOut[1][i];
                    p3 = fieldsOut[2][i];
                    p4 = fieldsOut[3][i];
                    p5 = fieldsOut[4][i];
                    p6 = fieldsOut[5][i];
                    p7 = fieldsOut[6][i];
                    p8 = fieldsOut[7][i];

                    p1x = fieldsOut_dx[0][i];
                    p2x = fieldsOut_dx[1][i];
                    p3x = fieldsOut_dx[2][i];
                    p4x = fieldsOut_dx[3][i];
                    p5x = fieldsOut_dx[4][i];
                    p6x = fieldsOut_dx[5][i];
                    p7x = fieldsOut_dx[6][i];
                    p8x = fieldsOut_dx[7][i];

                    p1y = fieldsOut_dy[0][i];
                    p2y = fieldsOut_dy[1][i];
                    p3y = fieldsOut_dy[2][i];
                    p4y = fieldsOut_dy[3][i];
                    p5y = fieldsOut_dy[4][i];
                    p6y = fieldsOut_dy[5][i];
                    p7y = fieldsOut_dy[6][i];
                    p8y = fieldsOut_dy[7][i];

                    p1z = fieldsOut_dz[0][i];
                    p2z = fieldsOut_dz[1][i];
                    p3z = fieldsOut_dz[2][i];
                    p4z = fieldsOut_dz[3][i];
                    p5z = fieldsOut_dz[4][i];
                    p6z = fieldsOut_dz[5][i];
                    p7z = fieldsOut_dz[6][i];
                    p8z = fieldsOut_dz[7][i];



                    J0x = - p1 * p2x + p2 * p1x - p3 * p4x + p4 * p3x - p5 * p6x + p6 * p5x - p7 * p8x + p8 * p7x ;
                    J0y = - p1 * p2y + p2 * p1y - p3 * p4y + p4 * p3y - p5 * p6y + p6 * p5y - p7 * p8y + p8 * p7y ;
                    J0z = - p1 * p2z + p2 * p1z - p3 * p4z + p4 * p3z - p5 * p6z + p6 * p5z - p7 * p8z + p8 * p7z ;

                    J1x = - p1 * p4x + p2 * p3x - p3 * p2x + p4 * p1x - p5 * p8x + p6 * p7x - p7 * p6x + p8 * p5x ;
                    J1y = - p1 * p4y + p2 * p3y - p3 * p2y + p4 * p1y - p5 * p8y + p6 * p7y - p7 * p6y + p8 * p5y ;
                    J1z = - p1 * p4z + p2 * p3z - p3 * p2z + p4 * p1z - p5 * p8z + p6 * p7z - p7 * p6z + p8 * p5z ;

                    J2x = p1 * p3x + p2 * p4x - p3 * p1x - p4 * p2x + p5 * p7x + p6 * p8x - p7 * p5x - p8 * p6x ;
                    J2y = p1 * p3y + p2 * p4y - p3 * p1y - p4 * p2y + p5 * p7y + p6 * p8y - p7 * p5y - p8 * p6y ;
                    J2z = p1 * p3z + p2 * p4z - p3 * p1z - p4 * p2z + p5 * p7z + p6 * p8z - p7 * p5z - p8 * p6z ;

                    J3x = - p1 * p2x + p2 * p1x + p3 * p4x - p4 * p3x - p5 * p6x + p6 * p5x + p7 * p8x - p8 * p7x ;
                    J3y = - p1 * p2y + p2 * p1y + p3 * p4y - p4 * p3y - p5 * p6y + p6 * p5y + p7 * p8y - p8 * p7y ;
                    J3z = - p1 * p2z + p2 * p1z + p3 * p4z - p4 * p3z - p5 * p6z + p6 * p5z + p7 * p8z - p8 * p7z ;

                    Jemx=J0x + n1 * J1x + n2 * J2x + n3 * J3x;
                    Jemy=J0y + n1 * J1y + n2 * J2y + n3 * J3y;
                    Jemz=J0z + n1 * J1z + n2 * J2z + n3 * J3z;



                //Potential Terms
                    ke = energyv[i];
                    pe = -mu_1_sq * (pow(p1, 2) + pow(p2, 2) + pow(p3, 2) + pow(p4, 2)) - mu_2_sq * (pow(p5, 2) + pow(p6, 2) + pow(p7, 2) + pow(p8, 2)) + lambda_1 * pow((pow(p1, 2) + pow(p2, 2) + pow(p3, 2) + pow(p4, 2)), 2) + lambda_2 * pow((pow(p5, 2) + pow(p6, 2) + pow(p7, 2) + pow(p8, 2)), 2) + lambda_3 * ((pow(p1, 2) + pow(p2, 2) + pow(p3, 2) + pow(p4, 2)) * (pow(p5, 2) + pow(p6, 2) + pow(p7, 2) + pow(p8, 2))) + (l4_m_l5)*pow((p1 * p5 + p2 * p6 + p3 * p7 + p4 * p8), 2) + (l4_p_l5)*pow((p1 * p6 - p2 * p5 + p3 * p8 - p4 * p7), 2);

                    energy = ke + pe;






                    
                    //cout <<p1z<<"-"<<p2z<<"-" <<p3z<<"-" <<p4z<<"-"<<p5z<<"-" <<p6z<<"-"<<p7z<<"-" <<p8z<< endl;


                    finalFields << R0 << " " << R1 << " " << R2 << " " << R3 << " " << R4 << " " << R5 << " " << n1 << " " << n2 << " " << n3 <<" " << p1 << " " << p2 << " " << p3 << " " << p4 << " " << p5 << " " << p6 << " " << p7 << " " << p8 << " " << Jemx << " " <<Jemy << " " << Jemz << " " << energy <<  " " << ke <<  " " << pe << endl;

                }



            }

            else {



                for (comp = 0; comp < nb_fields; comp++) {

                    MPI_Send(&fields[comp][frontHaloSize], coreSize, MPI_DOUBLE, 0, comp, MPI_COMM_WORLD);

                }
            }





        }

        //Gif Output

        if (makeGif and TimeStep % saveFreq == 0) {


            if (rank == 0) {
                // Create files
                string TimeStepPath = dir_path + "fields_timestep=" + to_string(TimeStep) + outTag + ".txt";
                ofstream Gif(TimeStepPath.c_str());
                Gif << "R0" << " " << "R1" << " " << "R2" << " " << "R3" << " " << "R4" << " " << "R5" << " " << "n1" << " " << "n2" << " " << "n3" << " " << "psi1" << " " << "psi2" << " " << "psi3" << " " << "psi4" << " " << "psi5" << " " << "psi6" << " " << "psi7" << " " << "psi8" << " " << "Jemx" << " " <<"Jemy" << " " << "Jemz" << " " << "E" <<  " " << "KE" <<  " " << "PE" <<endl;

                vector<vector<double>> fieldsOutnt(nb_fields, vector<double>(nPos, 0.0));
                vector<vector<double>> fieldsOutnt_dx(nb_fields, vector<double>(nPos, 0.0));
                vector<vector<double>> fieldsOutnt_dy(nb_fields, vector<double>(nPos, 0.0));
                vector<vector<double>> fieldsOutnt_dz(nb_fields, vector<double>(nPos, 0.0));
                vector<vector<double>> fieldsOutnt_dt(nb_fields, vector<double>(nPos, 0.0));
                double R0nt;
                double R1nt;
                double R2nt;
                double R3nt;
                double R4nt;
                double R5nt;

                double n1nt;
                double n2nt;
                double n3nt;

                double psi1;
                double psi2;
                double psi3;
                double psi4;
                double psi5;
                double psi6;
                double psi7;
                double psi8;
                
                double psi1x;
                double psi2x;
                double psi3x;
                double psi4x;
                double psi5x;
                double psi6x;
                double psi7x;
                double psi8x;

                double psi1y;
                double psi2y;
                double psi3y;
                double psi4y;
                double psi5y;
                double psi6y;
                double psi7y;
                double psi8y;

                double psi1z;
                double psi2z;
                double psi3z;
                double psi4z;
                double psi5z;
                double psi6z;
                double psi7z;
                double psi8z;
                
                double psi1t;
                double psi2t;
                double psi3t;
                double psi4t;
                double psi5t;
                double psi6t;
                double psi7t;
                double psi8t;

                
                double J0_x;
                double J0_y;
                double J0_z;
                
                double J1_x;
                double J1_y;
                double J1_z;
                
                double J2_x;
                double J2_y;
                double J2_z;
                
                double J3_x;
                double J3_y;
                double J3_z;
                
                double Jem_x;
                double Jem_y;
                double Jem_z;

                double ken;
                double pen;
                double energyn=0;
                vector<double> energynt(nPos, 0.0);

                for (comp = 0; comp < nb_fields; comp++) {

                    for (j = 0; j < coreSize; j++) { 
                        fieldsOutnt[comp][j] = fields[comp][frontHaloSize + j]; 
                        
                        fieldsOutnt_dx[comp][j] = fields_dx[comp][frontHaloSize + j]; 
                        fieldsOutnt_dy[comp][j] = fields_dy[comp][frontHaloSize + j]; 
                        fieldsOutnt_dz[comp][j] = fields_dz[comp][frontHaloSize + j]; 
                        fieldsOutnt_dt[comp][j] = fields_dt[comp][frontHaloSize + j]; 
                        //energynt[j] =+ (pow(fieldsOutnt_dx[comp][j], 2) + pow(fieldsOutnt_dy[comp][j], 2) + pow(fieldsOutnt_dz[comp][j], 2) + pow(fieldsOutnt_dt[comp][j], 2)) * dx * dy * dz; // Kinetic energy
                        //cout <<j<<" "<<rank<<" "<< size<<" " <<fieldsOutnt[1][j]<<" " << fieldsOutnt_dt[1][j]<<" " << fieldsOutnt_dx[1][j]<<" " << fieldsOutnt_dy[1][j]<<" " << fieldsOutnt_dz[1][j]<<" " <<energynt[j]<< endl;


                        }

                    for (j = 1; j < size; j++) {

                        int localCoreStartnt;
                        int localCoreSizent;
                        if (j < chunkRem) { localCoreStartnt = j * (chunk + 1); localCoreSizent = chunk + 1; }
                        else { localCoreStartnt = j * chunk + chunkRem; localCoreSizent = chunk; }

                        MPI_Recv(&fieldsOutnt[comp][localCoreStartnt], localCoreSizent, MPI_DOUBLE, j, comp, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        
                        MPI_Recv(&fieldsOutnt_dx[comp][localCoreStartnt], localCoreSizent, MPI_DOUBLE, j, comp, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&fieldsOutnt_dy[comp][localCoreStartnt], localCoreSizent, MPI_DOUBLE, j, comp, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&fieldsOutnt_dz[comp][localCoreStartnt], localCoreSizent, MPI_DOUBLE, j, comp, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&fieldsOutnt_dt[comp][localCoreStartnt], localCoreSizent, MPI_DOUBLE, j, comp, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        

                    }
                }
                //cout <<fieldsOutnt[0][0]<<"-"<<fieldsOutnt[1][0]<<"-" <<fieldsOutnt[2][0]<<"-"<< fieldsOutnt[3][0]<<"-"<< fieldsOutnt[4][0]<<"-"<< fieldsOutnt[5][0]<<"-"<< fieldsOutnt[6][0]<<"-"<< fieldsOutnt[7][0]<< endl;

                for (j = 0; j < nPos; j++) {
                    //cout <<frontHaloSize + j<<" " <<fields[comp][frontHaloSize + j]<<" " << fields_dt[comp][frontHaloSize + j]<<" " << fields_dx[comp][frontHaloSize + j]<<" " << fields_dy[comp][frontHaloSize + j]<<" " << fields_dz[comp][frontHaloSize + j]<<" " <<energynt[j]<< endl;
                    //cout <<j<<" "<<rank<<" "<< size<<" " <<fieldsOutnt[1][j]<<" " << fieldsOutnt_dt[1][j]<<" " << fieldsOutnt_dx[1][j]<<" " << fieldsOutnt_dy[1][j]<<" " << fieldsOutnt_dz[1][j]<<" " <<energynt[j]<< endl;


                    R0nt = pow(fieldsOutnt[0][j], 2) + pow(fieldsOutnt[1][j], 2) + pow(fieldsOutnt[2][j], 2) + pow(fieldsOutnt[3][j], 2) + pow(fieldsOutnt[4][j], 2) + pow(fieldsOutnt[5][j], 2) + pow(fieldsOutnt[6][j], 2) + pow(fieldsOutnt[7][j], 2);
                    R1nt = 2 * (fieldsOutnt[0][j] * fieldsOutnt[4][j] + fieldsOutnt[1][j] * fieldsOutnt[5][j] + fieldsOutnt[2][j] * fieldsOutnt[6][j] + fieldsOutnt[3][j] * fieldsOutnt[7][j]);
                    R2nt = 2 * (fieldsOutnt[0][j] * fieldsOutnt[5][j] + fieldsOutnt[2][j] * fieldsOutnt[7][j] - fieldsOutnt[1][j] * fieldsOutnt[4][j] - fieldsOutnt[3][j] * fieldsOutnt[6][j]);
                    R3nt = pow(fieldsOutnt[0][j], 2) + pow(fieldsOutnt[1][j], 2) + pow(fieldsOutnt[2][j], 2) + pow(fieldsOutnt[3][j], 2) - pow(fieldsOutnt[4][j], 2) - pow(fieldsOutnt[5][j], 2) - pow(fieldsOutnt[6][j], 2) - pow(fieldsOutnt[7][j], 2);
                    R4nt = 2 * (fieldsOutnt[0][j] * fieldsOutnt[6][j] - fieldsOutnt[1][j] * fieldsOutnt[7][j] - fieldsOutnt[2][j] * fieldsOutnt[4][j] + fieldsOutnt[3][j] * fieldsOutnt[5][j]);
                    R5nt = 2 * (fieldsOutnt[0][j] * fieldsOutnt[7][j] + fieldsOutnt[1][j] * fieldsOutnt[6][j] - fieldsOutnt[2][j] * fieldsOutnt[5][j] - fieldsOutnt[3][j] * fieldsOutnt[4][j]);

                    n1nt = -2 * (fieldsOutnt[0][j] * fieldsOutnt[2][j] + fieldsOutnt[1][j] * fieldsOutnt[3][j] + fieldsOutnt[4][j] * fieldsOutnt[6][j] + fieldsOutnt[5][j] * fieldsOutnt[7][j]);
                    n2nt = -2 * (fieldsOutnt[0][j] * fieldsOutnt[3][j] - fieldsOutnt[1][j] * fieldsOutnt[2][j] + fieldsOutnt[4][j] * fieldsOutnt[7][j] - fieldsOutnt[5][j] * fieldsOutnt[6][j]);
                    n3nt = -1 * (pow(fieldsOutnt[0][j], 2) + pow(fieldsOutnt[1][j], 2) - pow(fieldsOutnt[2][j], 2) - pow(fieldsOutnt[3][j], 2) + pow(fieldsOutnt[4][j], 2) + pow(fieldsOutnt[5][j], 2) - pow(fieldsOutnt[6][j], 2) - pow(fieldsOutnt[7][j], 2));

                    psi1 = fieldsOutnt[0][j];
                    psi2 = fieldsOutnt[1][j];
                    psi3 = fieldsOutnt[2][j];
                    psi4 = fieldsOutnt[3][j];
                    psi5 = fieldsOutnt[4][j];

                    psi6 = fieldsOutnt[5][j];
                    psi7 = fieldsOutnt[6][j];
                    psi8 = fieldsOutnt[7][j];

                    psi1x = fieldsOutnt_dx[0][j];
                    psi2x = fieldsOutnt_dx[1][j];
                    psi3x = fieldsOutnt_dx[2][j];
                    psi4x = fieldsOutnt_dx[3][j];
                    psi5x = fieldsOutnt_dx[4][j];
                    psi6x = fieldsOutnt_dx[5][j];
                    psi7x = fieldsOutnt_dx[6][j];
                    psi8x = fieldsOutnt_dx[7][j];

                    psi1y = fieldsOutnt_dy[0][j];
                    psi2y = fieldsOutnt_dy[1][j];
                    psi3y = fieldsOutnt_dy[2][j];
                    psi4y = fieldsOutnt_dy[3][j];
                    psi5y = fieldsOutnt_dy[4][j];
                    psi6y = fieldsOutnt_dy[5][j];
                    psi7y = fieldsOutnt_dy[6][j];
                    psi8y = fieldsOutnt_dy[7][j];

                    psi1z = fieldsOutnt_dz[0][j];
                    psi2z = fieldsOutnt_dz[1][j];
                    psi3z = fieldsOutnt_dz[2][j];
                    psi4z = fieldsOutnt_dz[3][j];
                    psi5z = fieldsOutnt_dz[4][j];
                    psi6z = fieldsOutnt_dz[5][j];
                    psi7z = fieldsOutnt_dz[6][j];
                    psi8z = fieldsOutnt_dz[7][j];

                    psi1t = fieldsOutnt_dt[0][j];
                    psi2t = fieldsOutnt_dt[1][j];
                    psi3t = fieldsOutnt_dt[2][j];
                    psi4t = fieldsOutnt_dt[3][j];
                    psi5t = fieldsOutnt_dt[4][j];
                    psi6t = fieldsOutnt_dt[5][j];
                    psi7t = fieldsOutnt_dt[6][j];
                    psi8t = fieldsOutnt_dt[7][j];

                    J0_x = - psi1 * psi2x + psi2 * psi1x - psi3 * psi4x + psi4 * psi3x - psi5 * psi6x + psi6 * psi5x - psi7 * psi8x + psi8 * psi7x ;
                    J0_y = - psi1 * psi2y + psi2 * psi1y - psi3 * psi4y + psi4 * psi3y - psi5 * psi6y + psi6 * psi5y - psi7 * psi8y + psi8 * psi7y ;
                    J0_z = - psi1 * psi2z + psi2 * psi1z - psi3 * psi4z + psi4 * psi3z - psi5 * psi6z + psi6 * psi5z - psi7 * psi8z + psi8 * psi7z ;

                    J1_x = - psi1 * psi4x + psi2 * psi3x - psi3 * psi2x + psi4 * psi1x - psi5 * psi8x + psi6 * psi7x - psi7 * psi6x + psi8 * psi5x ;
                    J1_y = - psi1 * psi4y + psi2 * psi3y - psi3 * psi2y + psi4 * psi1y - psi5 * psi8y + psi6 * psi7y - psi7 * psi6y + psi8 * psi5y ;
                    J1_z = - psi1 * psi4z + psi2 * psi3z - psi3 * psi2z + psi4 * psi1z - psi5 * psi8z + psi6 * psi7z - psi7 * psi6z + psi8 * psi5z ;

                    J2_x = psi1 * psi3x + psi2 * psi4x - psi3 * psi1x - psi4 * psi2x + psi5 * psi7x + psi6 * psi8x - psi7 * psi5x - psi8 * psi6x ;
                    J2_y = psi1 * psi3y + psi2 * psi4y - psi3 * psi1y - psi4 * psi2y + psi5 * psi7y + psi6 * psi8y - psi7 * psi5y - psi8 * psi6y ;
                    J2_z = psi1 * psi3z + psi2 * psi4z - psi3 * psi1z - psi4 * psi2z + psi5 * psi7z + psi6 * psi8z - psi7 * psi5z - psi8 * psi6z ;

                    J3_x = - psi1 * psi2x + psi2 * psi1x + psi3 * psi4x - psi4 * psi3x - psi5 * psi6x + psi6 * psi5x + psi7 * psi8x - psi8 * psi7x ;
                    J3_y = - psi1 * psi2y + psi2 * psi1y + psi3 * psi4y - psi4 * psi3y - psi5 * psi6y + psi6 * psi5y + psi7 * psi8y - psi8 * psi7y ;
                    J3_z = - psi1 * psi2z + psi2 * psi1z + psi3 * psi4z - psi4 * psi3z - psi5 * psi6z + psi6 * psi5z + psi7 * psi8z - psi8 * psi7z ;

                    // Jem_x=J0_x + n1nt * J1_x + n2nt * J2_x + n3nt * J3_x;
                    // Jem_y=J0_y + n1nt * J1_y + n2nt * J2_y + n3nt * J3_y;
                    // Jem_z=J0_z + n1nt * J1_z + n2nt * J2_z + n3nt * J3_z;

                    Jem_x=(psi8x+psi7x+psi6x+psi5x+psi4x+psi3x+psi2x+psi1x)* dx * dy * dz;
                    Jem_y=(psi8y+psi7y+psi6y+psi5y+psi4y+psi3y+psi2y+psi1y)* dx * dy * dz;
                    Jem_z=(psi8t*psi8t+psi7t*psi7t+psi6t*psi6t+psi5t*psi5t+psi4t*psi4t+psi3t*psi3t+psi2t*psi2t+psi1t*psi1t)* dx * dy * dz;

                    ken = (psi8x*psi8x+psi8y*psi8y+psi8z*psi8z+psi8t*psi8t+psi7x*psi7x+psi7y*psi7y+psi7z*psi7z+psi7t*psi7t+psi6x*psi6x+psi6y*psi6y+psi6z*psi6z+psi6t*psi6t+psi5x*psi5x+psi5y*psi5y+psi5z*psi5z+psi5t*psi5t+psi4x*psi4x+psi4y*psi4y+psi4z*psi4z+psi4t*psi4t+psi3x*psi3x+psi3y*psi3y+psi3z*psi3z+psi3t*psi3t+psi2x*psi2x+psi2y*psi2y+psi2z*psi2z+psi2t*psi2t+psi1x*psi1x+psi1y*psi1y+psi1z*psi1z+psi1t*psi1t)* dx * dy * dz;

                    pen = -mu_1_sq * (pow(psi1, 2) + pow(psi2, 2) + pow(psi3, 2) + pow(psi4, 2)) - mu_2_sq * (pow(psi5, 2) + pow(psi6, 2) + pow(psi7, 2) + pow(psi8, 2)) + lambda_1 * pow((pow(psi1, 2) + pow(psi2, 2) + pow(psi3, 2) + pow(psi4, 2)), 2) + lambda_2 * pow((pow(psi5, 2) + pow(psi6, 2) + pow(psi7, 2) + pow(psi8, 2)), 2) + lambda_3 * ((pow(psi1, 2) + pow(psi2, 2) + pow(psi3, 2) + pow(psi4, 2)) * (pow(psi5, 2) + pow(psi6, 2) + pow(psi7, 2) + pow(psi8, 2))) + (l4_m_l5)*pow((psi1 * psi5 + psi2 * psi6 + psi3 * psi7 + psi4 * psi8), 2) + (l4_p_l5)*pow((psi1 * psi6 - psi2 * psi5 + psi3 * psi8 - psi4 * psi7), 2);
                    energyn = ken + pen;



                    //cout <<psi1z<<"-"<<psi2z<<"-" <<psi3z<<"-" <<psi4z<<"-"<<psi5z<<"-" <<psi6z<<"-"<<psi7z<<"-" <<psi8z<< endl;
                    Gif << R0nt << " " << R1nt << " " << R2nt << " " << R3nt << " " << R4nt << " " << R5nt << " " << n1nt << " " << n2nt << " " << n3nt  << " " << psi1 << " " << psi2 << " " << psi3 << " " << psi4 << " " << psi5 << " " << psi6 << " " << psi7 << " " << psi8 << " " << Jem_x << " " <<Jem_y << " " << Jem_z << " " << energyn <<  " " << ken <<  " " << pen << endl;

                }

            }

            else {

                for (comp = 0; comp < nb_fields; comp++) {

                    MPI_Send(&fields[comp][frontHaloSize], coreSize, MPI_DOUBLE, 0, comp, MPI_COMM_WORLD);
                        
                    MPI_Send(&fields_dx[comp][frontHaloSize], coreSize, MPI_DOUBLE, 0, comp, MPI_COMM_WORLD);
                    MPI_Send(&fields_dy[comp][frontHaloSize], coreSize, MPI_DOUBLE, 0, comp, MPI_COMM_WORLD);
                    MPI_Send(&fields_dz[comp][frontHaloSize], coreSize, MPI_DOUBLE, 0, comp, MPI_COMM_WORLD);
                    MPI_Send(&fields_dt[comp][frontHaloSize], coreSize, MPI_DOUBLE, 0, comp, MPI_COMM_WORLD);

                }
            }

        }

	
	// Simulation Progress Output
        if (rank == 0 and TimeStep % countRate == 0) {
        
            cout << "\rTimestep " << TimeStep << " completed.";
        
        }

        // Barrier before going to the next timestep.
        MPI_Barrier(MPI_COMM_WORLD);

    }


    if (rank == 0) {

        cout << "\rTimestep " << nt << " completed." << endl;

        gettimeofday(&end, NULL);

        cout << "Time taken: " << end.tv_sec - start.tv_sec << "s" << endl;

    }


    // Deletes redundent outpur files if not used:
    if (rank == 0) {

        if (!finalOut) {
            finalFields.close();
            remove(finalFieldPath.c_str());
        }

        if (!calcEnergy and !wallDetect) { 
            valsPerLoop.close();
            remove(valsPerLoopPath.c_str());
        }
    }


    MPI_Finalize();


    return 0;
}
