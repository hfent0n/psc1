// Translate this file with
//
// g++ -O3 assignment-code.cpp -o assignment-code
//
// Run it with
//
// ./demo-code
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2020 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>

#include <cmath>
#include <vector>
#include <iterator>
using namespace std;

double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double* x0;
double* x1;
double* x2;

/**
 * Equivalent to x storing the velocities.
 */
double* v0;
double* v1;
double* v2;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;

double C;
/**
 * Set up scenario from the command line.
 *
 * If you need additional helper data structures, you can
 * initialise them here. Alternatively, you can introduce a
 * totally new function to initialise additional data fields and
 * call this new function from main after setUp(). Either way is
 * fine.
 *
 * This operation's semantics is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-4) / 7;

  C = 0.01 / NumberOfBodies; 
 
  x0    = new double [NumberOfBodies];
  x1    = new double [NumberOfBodies];
  x2    = new double [NumberOfBodies];

  v0    = new double [NumberOfBodies];
  v1    = new double [NumberOfBodies];
  v2    = new double [NumberOfBodies];

  mass  = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    
    x0[i] = std::stof(argv[readArgument]); readArgument++;
    x1[i] = std::stof(argv[readArgument]); readArgument++;
    x2[i] = std::stof(argv[readArgument]); readArgument++;

    v0[i] = std::stof(argv[readArgument]); readArgument++;
    v1[i] = std::stof(argv[readArgument]); readArgument++;
    v2[i] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  
  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x0[i]
        << " "
        << x1[i]
        << " "
        << x2[i]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}

void printMass(){
  cout << "[";
  for (int i=0; i<NumberOfBodies;i++){
    if (i + 1 != NumberOfBodies) cout << mass[i] << ", ";
    else cout << mass[i];
  }
  cout << ']' << endl;
}
void printPos(){
  for (int i = 0; i < NumberOfBodies; i++){
    cout << "[";
    cout << x0[i] << ", " << x1[i] << ", " << x2[i];
    cout << ']' << endl;
  }
 
}

void printVel(){
  
  for (int i = 0; i < NumberOfBodies; i++){
    cout << "[";
    cout << v0[i] << ", " << v1[i] << ", " << v2[i];
    cout << ']' << endl;
  }
 
}

/**
 * This is the main operation you should change in the assignment. You might
 * want to add a few more variables or helper functions, but this is where the
 * magic happens.
 */

void updateBody() {
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();
  bool touched = false;
  // force0 = force along x direction
  // force1 = force along y direction
  // force2 = force along z direction
  double* force0 = new double[NumberOfBodies];
  double* force1 = new double[NumberOfBodies];
  double* force2 = new double[NumberOfBodies];

  for (int i = 0; i < NumberOfBodies; i++ ){
    force0[i] = 0.0;
    force1[i] = 0.0;
    force2[i] = 0.0;
  }

  
  int mergeCount=0;
  int toBeMerged[NumberOfBodies][NumberOfBodies];
  for (int i = 0; i < NumberOfBodies; i++){
    for (int j = 0; j < NumberOfBodies; j++){
      toBeMerged[i][j] = 0;
    }
  }

  //create copies of positions  
  double x0_c[NumberOfBodies];
  double x1_c[NumberOfBodies];
  double x2_c[NumberOfBodies];

  memcpy(x0_c, x0, NumberOfBodies*sizeof(double));
  memcpy(x1_c, x1, NumberOfBodies*sizeof(double));
  memcpy(x2_c, x2, NumberOfBodies*sizeof(double));

  //create copies of the velocities
  double v0_c[NumberOfBodies];
  double v1_c[NumberOfBodies];
  double v2_c[NumberOfBodies];
  
  memcpy(v0_c, v0, NumberOfBodies*sizeof(double));
  memcpy(v1_c, v1, NumberOfBodies*sizeof(double));
  memcpy(v2_c, v2, NumberOfBodies*sizeof(double));

  

  //calculate force at t
  for (int i = 0; i < NumberOfBodies; i++){
    #pragma omp simd reduction(min:minDx)
    for (int j = 0; j < NumberOfBodies; j++){
      
      const double distance = sqrt(
        (x0_c[i]-x0_c[j]) * (x0_c[i]-x0_c[j]) +
        (x1_c[i]-x1_c[j]) * (x1_c[i]-x1_c[j]) +
        (x2_c[i]-x2_c[j]) * (x2_c[i]-x2_c[j])
      );

      
      // x,y,z forces acting on particle i
      if (distance){
        force0[i] += (x0_c[j]-x0_c[i]) * mass[i]*mass[j] / (distance*distance*distance);
        force1[i] += (x1_c[j]-x1_c[i]) * mass[i]*mass[j] / (distance*distance*distance);
        force2[i] += (x2_c[j]-x2_c[i]) * mass[i]*mass[j] / (distance*distance*distance);
      }
      minDx = std::min( minDx, distance );
    }
  }

  // update positions for t + 0.5
  #pragma ivdep 
  for ( int i = 0; i<NumberOfBodies; i++){
    x0_c[i] = x0_c[i] + (timeStepSize * 0.5) * v0_c[i];
    x1_c[i] = x1_c[i] + (timeStepSize * 0.5) * v1_c[i];
    x2_c[i] = x2_c[i] + (timeStepSize * 0.5) * v2_c[i];
  }
  
  // update velocities for t + 0.5 
  #pragma ivdep 
  for ( int i = 0; i < NumberOfBodies; i++){
    v0_c[i] = v0_c[i] + (timeStepSize * 0.5) * force0[i] / mass[i];
    v1_c[i] = v1_c[i] + (timeStepSize * 0.5) * force1[i] / mass[i];
    v2_c[i] = v2_c[i] + (timeStepSize * 0.5) * force2[i] / mass[i];
  }


  //calculate forces at the t + 0.5
  for (int i = 0; i < NumberOfBodies; i++){
    #pragma omp simd reduction(min:minDx)
    for (int j = 0; j < NumberOfBodies; j++){
      
      const double distance = sqrt(
        (x0_c[i]-x0_c[j]) * (x0_c[i]-x0_c[j]) +
        (x1_c[i]-x1_c[j]) * (x1_c[i]-x1_c[j]) +
        (x2_c[i]-x2_c[j]) * (x2_c[i]-x2_c[j])
      );

      if (distance <= C * (mass[i]*mass[j])) {
        if (i<j){
          toBeMerged[i][j] = 1;
          touched = true;
          mergeCount++;
        }
        
      }
      // x,y,z forces acting on particle i
      if (distance){
        force0[i] += (x0_c[j]-x0_c[i]) * mass[i]*mass[j] / (distance*distance*distance);
        force1[i] += (x1_c[j]-x1_c[i]) * mass[i]*mass[j] / (distance*distance*distance);
        force2[i] += (x2_c[j]-x2_c[i]) * mass[i]*mass[j] / (distance*distance*distance);
      }
      minDx = std::min( minDx, distance );
    }
  }

  // update positions for t + 1
  #pragma ivdep 
  for ( int i = 0; i<NumberOfBodies; i++){
    x0[i] = x0[i] + timeStepSize * v0_c[i];
    x1[i] = x1[i] + timeStepSize * v1_c[i];
    x2[i] = x2[i] + timeStepSize * v2_c[i];
  }
  
  // update velocities for t + 1
  #pragma ivdep 
  for ( int i = 0; i < NumberOfBodies; i++){
    v0[i] = v0[i] + timeStepSize * force0[i] / mass[i];
    v1[i] = v1[i] + timeStepSize * force1[i] / mass[i];
    v2[i] = v2[i] + timeStepSize * force2[i] / mass[i];
  }


  // update merged velocities and positions. Set mass to 0
  if (touched){
    for (int i = 0; i<NumberOfBodies; i++){
      for (int j = 0; j < NumberOfBodies; j++){ 
        if (toBeMerged[i][j]){
          v0[i] = (mass[i] * v0[i] + mass[j] * v0[j])  / (mass[i] + mass[j]) ;
          v1[i] = (mass[i] * v1[i] + mass[j] * v1[j])  / (mass[i] + mass[j]) ;
          v2[i] = (mass[i] * v2[i] + mass[j] * v2[j])  / (mass[i] + mass[j]) ;
          
          x0[i] = ((mass[i] * x0[i])  + (mass[j] * x0[j])) / (mass[i] + mass[j]) ;
          x1[i] = ((mass[i] * x1[i])  + (mass[j] * x1[j])) / (mass[i] + mass[j]) ;
          x2[i] = ((mass[i] * x2[i])  + (mass[j] * x2[j])) / (mass[i] + mass[j]) ;

          mass[i] += mass[j];
          mass[j] = 0.0;
        }
      }
      
    }

    //swap all 0s to end of the array
    for (int anchor = 0, cur = 0; cur < NumberOfBodies; cur++) {
      if (mass[cur] != 0.0) {
        swap(mass[anchor], mass[cur]);
        swap(v0[anchor], v0[cur]);
        swap(v1[anchor], v1[cur]);
        swap(v2[anchor], v2[cur]);
        swap(x0[anchor], x0[cur]);
        swap(x1[anchor], x1[cur]);
        swap(x2[anchor], x2[cur]);
        anchor++; 
      }
    }

    printMass();
    printPos();
    printVel();
  }
  
  
  

  NumberOfBodies -= mergeCount;

  for (int i = 0; i < NumberOfBodies; i++){
    maxV = max( maxV, sqrt( v0[i]*v0[i] + v1[i]*v1[i] + v2[i]*v2[i] ) );
  } 


  t += timeStepSize;

  delete[] force0;
  delete[] force1;
  delete[] force2;
}



/**
 * Main routine.
 *
 * No major changes in assignment. You can add a few initialisation
 * or stuff if you feel the need to do so. But keep in mind that you
 * may not alter what the program plots to the terminal.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  

  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      // printParaviewSnapshot();
      // std::cout << "plot next snapshot"
    	// 	    << ",\t time step=" << timeStepCounter
    	// 	    << ",\t t="         << t
			// 	<< ",\t dt="        << timeStepSize
			// 	<< ",\t v_max="     << maxV
			// 	<< ",\t dx_min="    << minDx
			// 	<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: " << x0[0] << ", " << x1[0] << ", " << x2[0] << std::endl;

  closeParaviewVideoFile();

  return 0;
}