#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;


  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
    if( 
        (0==estimations.size())
        ||
        (0==ground_truth.size())
        ||
        (estimations.size() != ground_truth.size())
       )
    {
        std::cerr << "Error in inputs of estimation or ground_truth!" << std::endl;
    }
    
  //accumulate squared residuals
  VectorXd residual_sum(4);
  residual_sum << 0,0,0,0;
  for(unsigned int i=0; i < estimations.size(); ++i){
        VectorXd diff  = (estimations[i] - ground_truth[i]);
        VectorXd squared_residual = diff.array()*diff.array();
        residual_sum += squared_residual;
  }

  //calculate the mean
  VectorXd residual_mean = residual_sum/estimations.size();

  //calculate the squared root
  rmse = residual_mean.array().sqrt();


  //return the result
  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3,4);
  //recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  double range_squar = (px*px + py*py);
  //check division by zero
  if(0 != range_squar)
  {
      //TODO: px*px + py*py might be close to zero. What should be done in those cases?
      Hj << px/std::sqrt(range_squar),                    py/std::sqrt(range_squar),                      0,                            0,
            -py/range_squar,                              px/range_squar,                                 0,                            0,
            py*(vx*py-vy*px)/std::pow(range_squar, 1.5),  px*(vy*px-vx*py)/std::pow(range_squar, 1.5),    px/std::sqrt(range_squar),    py/std::sqrt(range_squar);     
              
  }
  else
  {
      std::cerr << "CalculateJacobian() - Error - Division by Zero" << std::endl;
  }
  
  
  //compute the Jacobian matrix

  return Hj;
}


Eigen::VectorXd Tools::ConvCartesianToPolar(const Eigen::VectorXd& x)
{
  double px = x[0];
  double py = x[1];
  double vx = x[2];
  double vy = x[3];

  double distance = px*px + py*py;
  VectorXd ret_val(3);
  if(0 != distance)
  {
    double ro = std::sqrt(distance);
    double phi = std::atan2(py,px);
    double ro_dot = (px*vx+py*vy)/(ro);
    ret_val << ro, phi, ro_dot;
  }
  else
  { 
    double ro_dot = 0.0;
    double phi = std::atan2(py,px);

    //if the speed is along x-axis or y-axis
    if(0==px)
    {
      ro_dot = vy;
    }
    else if (0 == py)
    {
      ro_dot = vx;
    }

    ret_val << 0, phi, ro_dot;
    // std::cerr << "ConvPolarToCartesian() - Error - Division by Zero" << std::endl;
  }

  return ret_val;
}

// make the value fall between -PI and PI.
double Tools::ClampAngleFromNegPiToPi(double x)
{
  //TODO: modulo by 2*PI, although this Udacity input data does not seem to require it
  if (x< -M_PI)
  {
    x+= M_PI;
  }
  else if( x > M_PI)
  {
    x-= M_PI;
  }
  else
  {
    //do nothing if value is between -PI to PI.
  }
  return x;
}