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
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

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

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

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
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
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
    for (int j = 0; j < 3; j++){
      if (j != 2) cout << x[i][j] << ", ";
      else cout << x[i][j];
    }
    cout << ']' << endl;
  }
 
}

void printVel(){
  
  for (int i = 0; i < NumberOfBodies; i++){
    cout << "[";
    for (int j = 0; j < 3; j++){
      if(j!=2) cout << v[i][j] << ", ";
      else cout << v[i][j];
    }
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

  vector<vector<int>> toBeMerged;
  for (int i = 0; i < NumberOfBodies; i++){
    for (int j = i + 1; j < NumberOfBodies; j++){
      const double distance = sqrt(
        (x[i][0]-x[j][0]) * (x[i][0]-x[j][0]) +
        (x[i][1]-x[j][1]) * (x[i][1]-x[j][1]) +
        (x[i][2]-x[j][2]) * (x[i][2]-x[j][2])
      );

      if (distance <= C * (mass[i]*mass[j])) {
        toBeMerged.push_back(vector<int> {i, j});
        touched = true;
      }

      double forceX = (x[j][0]-x[i][0]) * mass[i]*mass[j] / distance / distance / distance ;
      double forceY = (x[j][1]-x[i][1]) * mass[i]*mass[j] / distance / distance / distance ;
      double forceZ = (x[j][2]-x[i][2]) * mass[i]*mass[j] / distance / distance / distance ;
      // x,y,z forces acting on particle i
      force0[i] += forceX;
      force1[i] += forceY;
      force2[i] += forceZ;
      // equal and opposite force
      force0[j] -= forceX;
      force1[j] -= forceY;
      force2[j] -= forceZ;

      minDx = std::min( minDx, distance );

    }
  }
  
  

  // update positions
  for ( int i = 0; i<NumberOfBodies; i++){
    x[i][0] = x[i][0] + timeStepSize * v[i][0];
    x[i][1] = x[i][1] + timeStepSize * v[i][1];
    x[i][2] = x[i][2] + timeStepSize * v[i][2];
  }
  
  // update velocities
  for ( int i = 0; i < NumberOfBodies; i++){
    v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
    v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
    v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];
  }
  
  

  // update merged velocities and positions. Set mass to 0
  for (vector<int> t: toBeMerged){
    v[t[0]][0] = (mass[t[0]] * v[t[0]][0] + mass[t[1]] * v[t[1]][0])  / (mass[t[0]] + mass[t[1]]) ;
    v[t[0]][1] = (mass[t[0]] * v[t[0]][1] + mass[t[1]] * v[t[1]][1])  / (mass[t[0]] + mass[t[1]]) ;
    v[t[0]][2] = (mass[t[0]] * v[t[0]][2] + mass[t[1]] * v[t[1]][2])  / (mass[t[0]] + mass[t[1]]) ;
    
    x[t[0]][0] = ((mass[t[0]] * x[t[0]][0])  + (mass[t[1]] * x[t[1]][0])) / (mass[t[0]] + mass[t[1]]) ;
    x[t[0]][1] = ((mass[t[0]] * x[t[0]][1])  + (mass[t[1]] * x[t[1]][1])) / (mass[t[0]] + mass[t[1]]) ;
    x[t[0]][2] = ((mass[t[0]] * x[t[0]][2])  + (mass[t[1]] * x[t[1]][2])) / (mass[t[0]] + mass[t[1]]) ;

    mass[t[0]] += mass[t[1]];
    mass[t[1]] = 0.0;
  }
  

  //swap all 0s to end of the array
  
  for (int anchor = 0, cur = 0; cur < NumberOfBodies; cur++) {
    if (mass[cur] != 0.0) {
      swap(mass[anchor], mass[cur]);
      swap(v[anchor], v[cur]);
      swap(x[anchor], x[cur]);
      anchor++; 
    }
  }
  
  
  

  NumberOfBodies -= toBeMerged.size();

  
  
  for (int i = 0; i < NumberOfBodies; i++){
    maxV = std::max( maxV, std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] ) );
  }

  if (touched) {
    cout << minDx << endl;
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
      //printParaviewSnapshot();
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
  std::cout << "Position of first remaining object: " << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;

  closeParaviewVideoFile();

  return 0;
}