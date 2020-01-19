#include "kalman_filter.h"
#include <iostream>

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
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
  std::cout << "P_ in predict method"<< std::endl;
  std::cout << P_ << std::endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  std::cout << "update method before computation"<< std::endl;
  VectorXd z_pred = H_ * x_;
  std::cout << "z_pred" << std::endl;
  std::cout << z_pred << std::endl;
  VectorXd y = z - z_pred;
  std::cout << "y" << std::endl;
  std::cout << y << std::endl;
  MatrixXd Ht = H_.transpose();
  std::cout << "Ht" << std::endl;
  std::cout << Ht << std::endl;
  std::cout << "P_" << std::endl;
  std::cout << P_ << std::endl;
  std::cout << "H_" << std::endl;
  std::cout << H_ << std::endl;
  std::cout << "R_" << std::endl;
  std::cout << R_ << std::endl;  
  MatrixXd S = H_ * P_ * Ht + R_;
  std::cout << S << std::endl;
  MatrixXd Si = S.inverse();
  std::cout << "Si" << std::endl;
  std::cout << Si << std::endl;
  MatrixXd PHt = P_ * Ht;
  std::cout << "PHt" << std::endl;
  std::cout << PHt << std::endl;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  std::cout << "P_" << std::endl; 
  std::cout << P_ << std::endl; 
  std::cout << "I" << std::endl; 
  std::cout << I << std::endl;
  std::cout << "K" << std::endl;
  std::cout << K << std::endl;
  std::cout << "H_" << std::endl;
  std::cout << H_ << std::endl;
  P_ = (I - K * H_) * P_;
  std::cout << "P_ in update method"<< std::endl;
  std::cout << P_ << std::endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
    float rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  float phi = atan2(x_(1), x_(0));
  float rho_dot;
  if (fabs(rho) < 0.0001) {
    rho_dot = 0;
  } else {
    rho_dot = (x_(0)*x_(2) + x_(1)*x_(3))/rho;
  }
  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
