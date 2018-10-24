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
  double ax, ay, az, axNew, ayNew, azNew;


  ofstream outfile;
  if(this->solverType == 0){
    outfile.open("verlet.dat");
    cout << "Solving with Verlet..." << endl;
    //Velocity Verlet
    for(int i = 0; i < gridPoints; i++){
      for(heavenlyBody &current: bodies){

        //reset acceleration after each loop
        ax = ay = az = axNew = ayNew = azNew = 0.0;

        for(heavenlyBody &other: bodies){
          ax += current.newtonian_acceleration(other, 0);
          ay += current.newtonian_acceleration(other, 1);
          az += current.newtonian_acceleration(other, 2);
        }

        //write positions to file
        outfile << setw(30) << setprecision(12) << current.position[0];
        outfile << setw(30) << setprecision(12) << current.position[1];
        outfile << setw(30) << setprecision(12) << current.position[2];

        current.position[0] = current.position[0] + h*current.velocity[0] + h2Half*ax;
        current.position[1] = current.position[1] + h*current.velocity[1] + h2Half*ay;
        current.position[2] = current.position[2] + h*current.velocity[2] + h2Half*az;

        for(heavenlyBody &other: bodies)
        {
          axNew += current.newtonian_acceleration(other, 0);
          ayNew += current.newtonian_acceleration(other, 1);
          azNew += current.newtonian_acceleration(other, 2);
        }

        current.velocity[0] = current.velocity[0] + hHalf*(axNew + ax);
        current.velocity[1] = current.velocity[1] + hHalf*(ayNew + ay);
        current.velocity[2] = current.velocity[2] + hHalf*(azNew + az);



      }//end bodies loop
      outfile << endl;
    }//end integration loop


  } //end verlet

  else if(this->solverType == 1){
    outfile.open("euler.dat");
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
        cout << ax << endl;
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


      }//end bodies loop
      outfile << endl;
    }//end integration loop
  }//end euler
  outfile.close();


}
