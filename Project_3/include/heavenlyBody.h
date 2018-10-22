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
  // std::string name;

  //Initializers
  heavenlyBody();
  heavenlyBody(double m, double x, double y, double z, double vx, double vy, double vz);//, std::string newName);

  //Functions
  double distance(heavenlyBody &otherBody);
  double newtonian_acceleration(heavenlyBody &otherBody, const int dimPos);
  double GR_adjusted_acceleration(heavenlyBody &otherBody, const int dimPos);
  double angular_momentum2D();
  double position_angle2D();


};

#endif
