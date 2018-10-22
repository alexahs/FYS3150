#include "heavenlyBody.h"

double G = 4*M_PI*M_PI;
double c = 63241.0;      //AU per yr
heavenlyBody::heavenlyBody()
{
  mass = 1.0;
  position[0] = 1.0;
  position[1] = 0.0;
  position[2] = 0.0;
  velocity[0] = 0.0;
  velocity[1] = 0.0;
  velocity[2] = 0.0;
  Epot = 0.0;
  Ekin = 0.0;
}

heavenlyBody::heavenlyBody(double m, double x, double y, double z, double vx, double vy, double vz)
{
  mass = m;
  position[0] = x;
  position[1] = y;
  position[2] = z;
  velocity[0] = vx;
  velocity[1] = vy;
  velocity[2] = vz;
  Epot = 0.0;
  Ekin = 0.0;
}

double heavenlyBody::kinetic_energy()
{
  vx = this->velocity[0];
  vy = this->velocity[1];
  vz = this->velocity[2];
  return 0.5*this->mass*(vx*vx + vy*vy + vz*vz)
}

double heavenlyBody::distance(heavenlyBody &otherBody)
{
  double x1, x2, y1, y2, z1, z2, X, Y, Z;

  x1 = this->position[0];
  y1 = this->position[1];
  z1 = this->position[2];

  x2 = otherBody.position[0];
  y2 = otherBody.position[1];
  z2 = otherBody.position[2];

  X = x1-x2;
  Y = y1-y2;
  Z = z1-z2;

  return sqrt(X*X + Y*Y + Z*Z);
}
double heavenlyBody::newtonian_acceleration(heavenlyBody &otherBody, const int dim)
{
  double r = this->distance(otherBody);
  double dimPos = this->position[dim] - otherBody.position[dim];
  double otherM = otherBody.mass;
  if(r!=0) return -dimPos*G*otherM/(r*r*r);
  else return 0;
}

double heavenlyBody::GR_adjusted_acceleration(heavenlyBody &otherBody, const int dim)
{
  double l = this->angular_momentum2D();
  double r = this->distance(otherBody);
  if(r!=0) return this->newtonian_acceleration(otherBody, dim)*(1.0 + 3*l*l/(r*r*c*c));
  else return 0;
}


double heavenlyBody::angular_momentum2D()
{
  return this->position[0]*this->velocity[1] - this->position[1]*this->velocity[0];
}

double heavenlyBody::position_angle2D()
{
  return atan(this->position[1]/this->position[0]);
}
