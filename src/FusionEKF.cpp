#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include "json.hpp"
#include <fstream>
#include <iostream>

using json = nlohmann::json;
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
	is_initialized_ = false;

	// Load parameters from json file
	std::ostringstream param_buf; 
	std::ifstream param_file("../data/config.json"); 
	param_buf << param_file.rdbuf(); 
	auto param = json::parse(param_buf.str());
	
	previous_timestamp_ = 0;
	
	laser_enabled_ = param["EKF"]["laser_enabled"].get<int>();
	radar_enabled_ = param["EKF"]["radar_enabled"].get<int>();

	// initializing matrices
	VectorXd x(4);
	x << 0, 0, 0, 0;
	
	//measurement covariance matrix - laser	
	R_laser_ = MatrixXd(2, 2);
	R_laser_ << param["EKF"]["noise_laser11"].get<double>(), 0,
			    0, param["EKF"]["noise_laser22"].get<double>();	
			
	//measurement covariance matrix - radar
	R_radar_ = MatrixXd(3, 3);		
	R_radar_ << param["EKF"]["noise_radar11"], 0, 0,
			    0, param["EKF"]["noise_radar22"], 0,
				0, 0, param["EKF"]["noise_radar33"];
							
	H_laser_ = MatrixXd(2, 4);
	H_laser_ << 1, 0, 0, 0,
			    0, 1, 0, 0;
				
	Hj_ = MatrixXd(3, 4);

	// Initialize process noise for process noise covariance matrix (Q)
	noise_ax_ = param["EKF"]["noise_ax"].get<double>();
	noise_ay_ = param["EKF"]["noise_ax"].get<double>();
	MatrixXd Q(4,4);
	Q << 0, 0, 0, 0,
		 0, 0, 0, 0,
		 0, 0, 0, 0,
		 0, 0, 0, 0;
		
	MatrixXd P(4,4);
	P << 1, 0, 0, 0,
	     0, 1, 0, 0,
	     0, 0, 1000, 0,
	     0, 0, 0, 1000;
		  
	MatrixXd F(4,4);
	F << 1, 0, 1, 0,
	     0, 1, 0, 1,
	     0, 0, 1, 0,
	     0, 0, 0, 1;

	// Initialize EKF
	ekf_.Init(x, P, F, H_laser_, R_laser_, Q);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
	// Skip RADAR or LASER sensor type if not configured
	if (radar_enabled_ == 0 && measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		return;
	}
	if (laser_enabled_ == 0 && measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		return;
	}
	
	/*****************************************************************************
	*  Initialization
	****************************************************************************/
	if (!is_initialized_) {
	/**
	  * Initialize the state ekf_.x_ with the first measurement.
	*/
	previous_timestamp_ = measurement_pack.timestamp_;
	
	// first measurement
	cout << "EKF: " << endl;
	
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		/**
		Convert radar from polar to cartesian coordinates and initialize state.
		*/
		VectorXd x_polar = measurement_pack.raw_measurements_;
		VectorXd x_cartesian = tools_.Polar2Cartesian(x_polar);
		ekf_.x_ = x_cartesian;
	}
	else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		/**
		Initialize state.
		*/
		VectorXd x_cartesian = measurement_pack.raw_measurements_;
		ekf_.x_ << x_cartesian(0), x_cartesian(1), 0, 0;
	}
	
	// done initializing, no need to predict or update
	is_initialized_ = true;
	std::cout << ekf_.x_ << std::endl;
	
	return;
	}

	/*****************************************************************************
	*  Prediction
	****************************************************************************/

	/**
	 * Update the state transition matrix F according to the new elapsed time.
	  - Time is measured in seconds.
	 * Update the process noise covariance matrix.
	*/
	// Update of F and Q handled in Predict method of ekf_ class
	float delta_T = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = measurement_pack.timestamp_;

	ekf_.Predict(delta_T, noise_ax_, noise_ay_);

	/*****************************************************************************
	*  Update
	****************************************************************************/

	/**
	 * Use the sensor type to perform the update step.
	 * Update the state and covariance matrices.
	*/

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		// Radar updates
		std::cout << "Radar" << std::endl;
		
		VectorXd x_polar = measurement_pack.raw_measurements_;
		x_polar(1) = tools_.WrapAngle(x_polar(1), 0, 2*M_PI);
		
		VectorXd x_cartesian = tools_.Polar2Cartesian(x_polar);	
		Hj_ = tools_.CalculateJacobian(x_cartesian);

		ekf_.H_ = Hj_;
		ekf_.R_ = R_radar_;
		ekf_.UpdateEKF(x_polar);
	} else {
		// Laser updates
		std::cout << "Laser" << std::endl;
		
		VectorXd x_cartesian = measurement_pack.raw_measurements_;
		
		ekf_.H_ = H_laser_;
		ekf_.R_ = R_laser_;
		ekf_.Update(x_cartesian);	
	}

	// print the output
	cout << "x_ = " << ekf_.x_ << endl;
	cout << "P_ = " << ekf_.P_ << endl;
}
