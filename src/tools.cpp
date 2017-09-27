#include <iostream>
#include "tools.h"
#include "math.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	/**
	* Calculate the RMSE here.
	*/
  	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	/**
	* Calculate a Jacobian here.
	*/
	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if(fabs(c1) < 0.0001){
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}

VectorXd Tools::Cartesian2Polar(const VectorXd& x_cartesian) {
	VectorXd x_polar(3);
	
	float px = x_cartesian(0);
	float py = x_cartesian(1);
	float vx = x_cartesian(2);
	float vy = x_cartesian(3);
	
	// Range
	float rho = sqrt(px*px+py*py);	
	if(fabs(rho) < 0.0001){
		std::cout << "Cartesian2Polar- Error - Division by Zero" << endl;
		return x_polar;
	}	
	
	// Angle
	float phi = atan(py/px);
	if (py > 0) {
		// 1st Quadrant
		if (px > 0){
			phi = phi;
		}
		// 2nd Quadrant
		else {
			phi = M_PI + phi;
		}
	}
	else {
		// 3rd Quadrant
		if (px <0) {
			//phi = -M_PI + phi;
			phi = M_PI + phi;
		}
		// 4th Quadrant
		else {
			//phi = phi;
			phi = 2*M_PI + phi;
		}
	}		
	
	// Range rate
	float drho= (px*vx + py*vy) / rho;
	
	x_polar << rho, phi, drho;
	return x_polar;
}

VectorXd Tools::Polar2Cartesian(const VectorXd& x_polar) {
	VectorXd x_cartesian(4);

	float rho = x_polar(0);
	float phi = x_polar(1);
	float drho = x_polar(2);

	float px = rho * cos(phi);
	float py = rho * sin(phi);
	float vx = drho * cos(phi) * (px>0 ? 1 : -1);
	float vy = drho * sin(phi) * (py>0 ? 1 : -1);

	x_cartesian << px, py, vx, vy;
	return x_cartesian;
}

float Tools::WrapAngle(const float angle, const float min, const float max) {
	if (min > max) {
		std::cout << "WrapAngle- Error - max smaller than min" << endl;
		return angle;
	}	
	
	float wrapped_angle = angle;
	while(wrapped_angle < min) {
		wrapped_angle += (max-min);
	}
	while(wrapped_angle > max) {
		wrapped_angle -= (max-min);
	}
	return wrapped_angle;
}
/*
VectorXd Tools::PolarTransformAngle(const VectorXd& x_polar) {
	VectorXd x_polar_new = x_polar;

	float phi = x_polar(1);
	if (phi < 0) {
		phi = 2*M_PI + phi;
	}
	x_polar_new(1) = phi;
	return x_polar_new;
}*/