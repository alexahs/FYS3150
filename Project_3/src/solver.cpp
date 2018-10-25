#include "solver.h"

solver::solver(){

  cout << "Cannot be initialized without arguments"<< endl;
  exit(1);
}

solver::solver(int n, double T, vector<heavenlyBody> planets, int m){
  gridPoints = n;
  simulationTime = T;
  bodies = planets;
  solverType = m;
}

void solver::solve(){
  //precalculated constants
  double h = simulationTime/((double) gridPoints);
  double hHalf = h/2.0;
  double h2Half = h*h/2.0;

  //temporary acceleration values
  double ax, ay, az;

  ofstream outputfile;
  // outputfile.open("data.dat");

  //Timer
  using std::chrono::steady_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  steady_clock::time_point programStart;
  programStart = steady_clock::now();

  if(this->solverType == 0){
    outputfile.open("verlet.dat");
    cout << "Solving with Verlet..." << endl;

    //Velocity Verlet
    for(int i = 0; i < gridPoints; i++){
      for(heavenlyBody &current: bodies){
        current.acceleration[0] = current.acceleration[1] = current.acceleration[2] =
        current.accelerationNew[0] = current.accelerationNew[1] = current.accelerationNew[2] = 0;
        for(heavenlyBody &other: bodies){
          current.acceleration[0] += current.newtonian_acceleration(other, 0);
          current.acceleration[1] += current.newtonian_acceleration(other, 1);
          current.acceleration[2] += current.newtonian_acceleration(other, 2);
          }
        }

        for(heavenlyBody &current: bodies){
        //write positions to file
        outputfile << setw(30) << setprecision(12) << current.position[0];
        outputfile << setw(30) << setprecision(12) << current.position[1];
        outputfile << setw(30) << setprecision(12) << current.position[2];

        current.position[0] = current.position[0] + h*current.velocity[0] + h2Half*current.acceleration[0];
        current.position[1] = current.position[1] + h*current.velocity[1] + h2Half*current.acceleration[1];
        current.position[2] = current.position[2] + h*current.velocity[2] + h2Half*current.acceleration[2];
      }


      for(heavenlyBody &current: bodies){
        for(heavenlyBody &other: bodies)
        {
          current.accelerationNew[0] += current.newtonian_acceleration(other, 0);
          current.accelerationNew[1] += current.newtonian_acceleration(other, 1);
          current.accelerationNew[2] += current.newtonian_acceleration(other, 2);
        }
      }


      for(heavenlyBody &current: bodies){
        current.velocity[0] = current.velocity[0] + hHalf*(current.acceleration[0] + current.accelerationNew[0]);
        current.velocity[1] = current.velocity[1] + hHalf*(current.acceleration[1] + current.accelerationNew[1]);
        current.velocity[2] = current.velocity[2] + hHalf*(current.acceleration[2] + current.accelerationNew[2]);
      }

      outputfile << endl;
    }//end integration loop
  } //end verlet

  else if(this->solverType == 1){
    outputfile.open("euler.dat");
    cout << "Solving with Euler..."<< endl;


    //Euler-Chromer
    for(int i = 0; i < gridPoints; i++){
      for(heavenlyBody &current: bodies){

        //reset acceleration after each loop
        ax = ay = az = 0.0;

        for(heavenlyBody &other: bodies){
          ax += current.newtonian_acceleration(other, 0);
          ay += current.newtonian_acceleration(other, 1);
          az += current.newtonian_acceleration(other, 2);
        }
        current.velocity[0] += h*ax;
        current.velocity[1] += h*ay;
        current.velocity[2] += h*az;

        current.position[0] += h*current.velocity[0];
        current.position[1] += h*current.velocity[1];
        current.position[2] += h*current.velocity[2];

        //write positions to filename
        outputfile << setw(30) << setprecision(12) << current.position[0];
        outputfile << setw(30) << setprecision(12) << current.position[1];
        outputfile << setw(30) << setprecision(12) << current.position[2];


      }//end bodies loop
      outputfile << endl;
    }//end integration loop
  }//end euler

  outputfile.close();
  //Print run time to terminal
  duration<double> programTime = duration_cast<duration<double>>(steady_clock::now() - programStart);
  cout << "Program time: " << programTime.count() << endl;

}
