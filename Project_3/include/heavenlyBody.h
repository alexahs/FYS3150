#ifndef HEAVENLYBODY_H
#define HEAVENLYBODY_H
#include <cmath>
#include <string>


class heavenlyBody{

public:

  //Properties
  double mass;
  double position[3];
  double velocity[3];
  double Epot;
  double Ekin;

  //Initializers
  heavenlyBody();
  heavenlyBody(double m, double x, double y, double z, double vx, double vy, double vz);

  //Methods
  double kinetic_energy();
  double distance(heavenlyBody &otherBody); //calculates the distnce r relative to another planet
  double newtonian_acceleration(heavenlyBody &otherBody, const int dimPos); //calculates the acceleration due to gravity of another planet
  double GR_adjusted_acceleration(heavenlyBody &otherBody, const int dimPos); //takes into account relativistic effects
  double angular_momentum2D(); //calculates angular momentum per mass in two dimensions
  double position_angle2D(); //calculates the angle of the planets position in the 2D plane for a given position


};

#endif
