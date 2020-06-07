#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculates the RMSE.
   */

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.empty() || estimations.size() != ground_truth.size()) {
    std::cout << "The estimations are empty or the ground truth and estimations sizes are not equal" << std::endl;
    return rmse;
  }

  for (unsigned int i = 0; i < estimations.size(); ++i) {
    rmse += (VectorXd) (estimations[i] - ground_truth[i]).array().square();
  }

  rmse /= estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
  /**
   * Calculate a Jacobian here.
   */

  MatrixXd Hj(3, 4);
  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // check division by zero
  if (px == 0 && py == 0) {
    std::cout << "CalculateJacobian() - Error - Division by Zero" << std::endl;
    return Hj;
  }

  // compute the Jacobian matrix
  double distance = px * px + py * py;
  double px_dist = px / sqrt(distance);
  double py_dist = py / sqrt(distance);
  double bottom_left = (vx * py - vy * px) / pow(distance, 3 / 2.);

  Hj << px_dist, py_dist, 0, 0,
      -py / distance, px / distance, 0, 0,
      py * bottom_left, -px * bottom_left, px_dist, py_dist;

  return Hj;

}
