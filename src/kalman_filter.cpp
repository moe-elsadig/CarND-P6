#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

  // predict the state
  x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

  // update the state by using Kalman Filter equations
  VectorXd y = z - (H_*x_);
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_*Ht;
  MatrixXd S = (H_*PHt) + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = (PHt*Si);

  //new estimate
  x_ = x_ + (K*y);
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = ((I - (K*H_))*P_);
  // Calculate h(x')
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  float hypo_ = sqrt(pow(px,2.)+pow(py,2.));

  VectorXd hx = VectorXd(3);
  if ( px < 0.0001 && py < 0.0001 )
    hx << 0.0, 0.0, 0.0;
  else
    hx << hypo_, atan2(py,px), ((px*vx+py*vy)/hypo_);

  // update the state by using Kalman Filter equations
  VectorXd y = z - hx;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_*Ht;
  MatrixXd S = (H_*PHt) + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = (PHt*Si);

  //new estimate
  x_ = x_ + (K*y);
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = ((I - (K*H_))*P_);
}
