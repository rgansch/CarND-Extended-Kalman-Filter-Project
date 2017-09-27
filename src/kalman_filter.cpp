#include <iostream>
#include "kalman_filter.h"
#include "tools.h"
#include "math.h"

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

void KalmanFilter::Predict(const float delta_T, const float noise_ax, const float noise_ay) {
	/**
	* predict the state
	*/	
	float c_T2 = pow(delta_T,2);
	float c_T3 = pow(delta_T,3)/2;
	float c_T4 = pow(delta_T,4)/4;
	
	F_(0, 2) = delta_T;
	F_(1, 3) = delta_T;
	
	Q_ << c_T4*noise_ax, 0, c_T3*noise_ax, 0,
		  0, c_T4*noise_ay, 0, c_T3*noise_ay,
		  c_T3*noise_ax, 0, c_T2*noise_ax, 0,
		  0, c_T3*noise_ay, 0, c_T2*noise_ay;

	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	/**
	* update the state by using Kalman Filter equations
	*/
  	VectorXd z_pred = H_ * x_;
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

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
	* update the state by using Extended Kalman Filter equations
	*/
	VectorXd x_polar = tools_.Cartesian2Polar(x_);
	// Wrap angle of x into range of z for proper subtraction
	x_polar(1) = tools_.WrapAngle(x_polar(1), z(1)-M_PI, z(1)+M_PI);
	VectorXd y = z - x_polar;
	
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ << x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ << (I - K * H_) * P_;
}
