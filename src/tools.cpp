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
   * TODO:
   * Calculate a Jacobian here.
   */
}
