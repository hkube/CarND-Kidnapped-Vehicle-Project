/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <limits>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

ParticleFilter::ParticleFilter()
: num_particles_(100)
, is_initialized_(false)
{}


void ParticleFilter::init(double x, double y, double theta, Sigma3 std) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates
  // of x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

//  std::cout << "ParticleFilter::init(" << x << ", " << y << ", " << theta
//            << ", [" << std.x << ", " << std.y << ", " << std.theta << "])" << std::endl;

  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(x, std.x);
  std::normal_distribution<double> dist_y(y, std.y);
  std::normal_distribution<double> dist_theta(theta, std.theta);
  for (unsigned i = 0; num_particles_ > i; i++)
  {
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;
    particles_.push_back(p);
  }
  is_initialized_ = true;
}

void ParticleFilter::prediction(double delta_t,
                                const ParticleFilter::Sigma3 std_pos,
                                double velocity,
                                double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  std::default_random_engine gen;
  for (unsigned i = 0; num_particles_ > i; i++)
  {
    Particle & p = particles_[i];
    double ofs_x = 0.0;
    double ofs_y = 0.0;
    double ofs_theta = 0.0;
    if (1e-3 < abs(yaw_rate))
    {
      ofs_x = velocity/yaw_rate * (sin(p.theta + yaw_rate*delta_t) - sin(p.theta));
      ofs_y = velocity/yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate*delta_t));
      ofs_theta = yaw_rate * delta_t;
    }
    else
    {
      ofs_x += velocity * delta_t * cos(p.theta);
      ofs_y += velocity * delta_t * sin(p.theta);
      // heading doesn't change!
    }
    std::normal_distribution<double> dist_x(ofs_x, std_pos.x);
    p.x += dist_x(gen);
    std::normal_distribution<double> dist_y(ofs_y, std_pos.y);
    p.y += dist_y(gen);
    std::normal_distribution<double> dist_theta(ofs_theta, std_pos.theta);
    p.theta += dist_theta(gen);
  }
}


void ParticleFilter::dataAssociation(//const std::vector<LandmarkObs> predicted,
                                     const Map &map_landmarks,
                                     std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  for (std::vector<LandmarkObs>::iterator it = observations.begin(); it != observations.end(); it++) {
    double minDist = std::numeric_limits<double>::max();
    for (unsigned lmIdx = 0; map_landmarks.landmark_list.size() > lmIdx; lmIdx++) {
      const auto lm = map_landmarks.landmark_list[lmIdx];
      const double dist = sqrt(pow((lm.x_f - it->x), 2.) + pow((lm.y_f - it->y), 2.));
      if (dist < minDist) {
        minDist = dist;
        it->id = lmIdx;
      }
    }
  }
}


void ParticleFilter::updateWeights(double sensor_range,
                                   const ParticleFilter::Sigma2 std_landmark,
                                   const std::vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

  for (auto & p : particles_) {

    std::vector<LandmarkObs> transObs;

    for (const auto o : observations) {
      LandmarkObs to;
      to.id = 0;
      to.x = o.x * cos(p.theta) - o.y * sin(p.theta) + p.x;
      to.y = o.x * sin(p.theta) + o.y * cos(p.theta) + p.y;
      transObs.push_back(to);
    }

    dataAssociation(map_landmarks, transObs);

    const auto lmListLen = map_landmarks.landmark_list.size();
    const double probDensNormalizer = 1 / (2 * M_PI * std_landmark.x * std_landmark.y);
    double weight = 1.;

    for (const auto o : transObs) {
      const auto lmIdx = o.id;
      if ((0 <= lmIdx) && (lmListLen > lmIdx)) {
        const auto & assLm = map_landmarks.landmark_list[lmIdx];
        const double probDens =
            probDensNormalizer *
            exp(
                -(  pow((o.x - assLm.x_f), 2.) / (2 * pow(std_landmark.x, 2.))
                  + pow((o.y - assLm.y_f), 2.) / (2 * pow(std_landmark.y, 2.)) )
                );
        weight *= probDens;
      }
    }
    p.weight = weight;
  }
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  std::discrete_distribution<unsigned> dist(particles_.begin(), particles_.end());
  std::default_random_engine generator;

  std::vector<Particle> new_particles;
  for (unsigned i = 0; particles_.size() > i; i++) {
    new_particles.push_back( particles_[ dist(generator) ]);
  }

  particles_ = new_particles;
}


ParticleFilter::Particle ParticleFilter::SetAssociations(ParticleFilter::Particle& particle,
                                                         const std::vector<int>& associations,
                                                         const std::vector<double>& sense_x,
                                                         const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world
    // coordinates mapping to associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}


string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
