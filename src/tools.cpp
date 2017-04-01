#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  // Calculate the RMSE here.
  VectorXd rmse(4);
	rmse << 0.,0.,0.,0.;

	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		std::cout << "Invalid estimation or ground_truth data" << std::endl;
		return rmse;
	}

	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	rmse = rmse/estimations.size();
	rmse = rmse.array().sqrt();
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  // Calculate a Jacobian here.
  MatrixXd Hj_(3,4);

	double p_x = x_state(0);
	double p_y = x_state(1);
	double v_x = x_state(2);
	double v_y = x_state(3);
  // std::cout << "px,py,vx,vy " << p_x << p_y << v_x << v_y << std::endl;

	if(p_x == 0 || p_y == 0){
	    std::cout << "Division by 0 Error!" << std::endl;
      return Hj_;
	}

  Hj_ << (p_x/sqrt(pow(p_x,2)+pow(p_y,2))),(p_y/sqrt(pow(p_x,2)+pow(p_y,2))),0,0,
    (-p_y/(pow(p_x,2)+pow(p_y,2))),(p_x/(pow(p_x,2)+pow(p_y,2))),0,0,
    (((p_y*(v_x*p_y-v_y*p_x))/(pow((pow(p_x,2)+pow(p_y,2)),(3/2))))),
    (((p_x*(v_y*p_x-v_x*p_y))/(pow((pow(p_x,2)+pow(p_y,2)),(3/2))))),
    (p_x/sqrt(pow(p_x,2)+pow(p_y,2))),
    (p_y/sqrt(pow(p_x,2)+pow(p_y,2)));

  return Hj_;
}
