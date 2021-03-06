//#include <iostream>
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * Predict the state
   */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * Updates the state by using Kalman Filter equations
   */
  MatrixXd y = z - H_ * x_;
  MatrixXd H_transposed = H_.transpose();

  MatrixXd S = H_ * P_ * H_transposed + R_;
  MatrixXd K = P_ * H_transposed * S.inverse();
  x_ = x_ + K * y;

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * Updates the state by using Extended Kalman Filter equations
   */
  VectorXd h(3);

  double dist = sqrt(x_(0) * x_(0) + x_(1) * x_(1));

  h << dist, atan2(x_(1), x_(0)), (x_(0) * x_(2) + x_(1) * x_(3)) / dist;

  double angle = z(1);
  double current_angle = h(1);

  if (abs(angle - current_angle) > M_PI / 2) {
    if (angle > current_angle) {
      angle -= M_PI * 2;
    } else {
      angle += M_PI * 2;
    }
  }

  VectorXd newZ(z.size());
  newZ << z(0), angle, z(2);

  MatrixXd y = newZ - h;

//  std::cout << "H_: " << h(1) << "\tz_: " << z(1) << "\tNewZ_: " << newZ(1) << "\tY_: " << y(1) << std::endl;

  MatrixXd H_transposed = H_.transpose();

  MatrixXd S = H_ * P_ * H_transposed + R_;
  MatrixXd K = P_ * H_transposed * S.inverse();
  x_ = x_ + K * y;

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  P_ = (I - K * H_) * P_;
}
