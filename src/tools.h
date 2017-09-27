#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
	/**
	* Constructor.
	*/
	Tools();

	/**
	* Destructor.
	*/
	virtual ~Tools();

	/**
	* A helper method to calculate RMSE.
	*/
	VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

	/**
	* A helper method to calculate Jacobians.
	*/
	MatrixXd CalculateJacobian(const VectorXd& x_state);

	/**
	* A helper method to convert polar to cartesian coordinates
	*/
	VectorXd Cartesian2Polar(const VectorXd& x_cartesian);
	
	/**
	* A helper method to convert cartesian to polar coordinates
	*/
	VectorXd Polar2Cartesian(const VectorXd& x_polar);
	
	/**
	* A helper method to wrap an angle into the range [min, max], e.g. [0, 2*M_PI]
	*/	
	float WrapAngle(const float angle, const float min, const float max);

};

#endif /* TOOLS_H_ */
