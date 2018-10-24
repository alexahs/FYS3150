#include "solver.h"


solver::solver(){

  gridPoints = 1000;
  simulationTime = 1.0;
}

solver::solver(int n, double T, vector <heavenlyBody> bodies, int m, string ofname){
  gridPoints = n;
  simulationTime = T;
  allBodies = bodies;
  solverType = m;
  outfileName = ofname;
}

void solver::solve(){
  //precalculated constants
  h = simulationTime/((double) gridPoints);
  hHalf = h/2.0;
  h2Half = h*h/2.0;

  //temporary acceleration values
  double ax, double ay, double az, double axNew, double ayNew, double azNew;
  ofstream outfile;
  outfile.open(outfileName);
  if(this->solverType == 0){
    //Velocity Verlet
    for(int i = 0; i < gridPoints; i++){
      for(heavenlyBody &current: allBodies){

        //reset acceleration after each loop
        ax = ay = az = axNew = ayNew = azNew = 0.0;

        for(heavenlyBody &other: allBodies){
          ax += current.newtonian_acceleration(other, 0);
          ay += current.newtonian_acceleration(other, 1);
          az += current.newtonian_acceleration(other, 2);
        }

        //write positions to file
        outfile << setw(30) << setprecision(12) << current.position[0];
        outfile << setw(30) << setprecision(12) << current.position[1];
        outfile << setw(30) << setprecision(12) << current.position[2];

        current.position[0] = current.position[0] + h*current.velocity[0] + h2Half*accelerationOld[0];
        current.position[1] = current.position[1] + h*current.velocity[1] + h2Half*accelerationOld[1];
        current.position[2] = current.position[2] + h*current.velocity[2] + h2Half*accelerationOld[2];

        for(heavenlyBody &other: allBodies)
        {
          axNew += current.newtonian_acceleration(other, 0);
          ayNew += current.newtonian_acceleration(other, 1);
          azNew += current.newtonian_acceleration(other, 2);
        }

        current.velocity[0] = current.velocity[0] + hHalf*(axNew + ax);
        current.velocity[1] = current.velocity[1] + hHalf*(ayNew + ay);
        current.velocity[2] = current.velocity[2] + hHalf*(azNew + az);



      }//end allBodies loop
      outfile << endl;
    }//end integration loop


  } //end verlet

  if(this->solverType == 1){
    //Euler-Chromer
    for(int i = 0; i < gridPoints; i++){
      for(heavenlyBody &current: allBodies){

        //reset acceleration after each loop
        ax = ay = az = 0.0;

        for(heavenlyBody &other: allBodies){
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
        outfile << setw(30) << setprecision(12) << current.position[0];
        outfile << setw(30) << setprecision(12) << current.position[1];
        outfile << setw(30) << setprecision(12) << current.position[2];


      }//end allBodies loop
      outfile << endl;
    }//end integration loop
  }//end euler
  outfile.close();


}
