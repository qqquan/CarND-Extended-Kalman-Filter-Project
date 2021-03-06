#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
    * Finish initializing the FusionEKF.
    * Set the process/prediction noise for the process covariance matrix Q
  */
  noise_ax = 7; 
  noise_ay = 8;


  VectorXd x_in(4);
  x_in << 0,0,0,0;

  MatrixXd P_in(4, 4);
  P_in << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;

  MatrixXd F_in(4, 4);
  F_in << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1; 

  MatrixXd Q_in(4, 4);
  Q_in << 0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0; 

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;     

  //init with LASAR data
  ekf_.Init(x_in, P_in, F_in, H_laser_, R_laser_, Q_in);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    double init_x = 0.0;
    double init_y = 0.0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];

      init_x = rho*std::cos(phi);
      init_y = rho*std::sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      init_x = measurement_pack.raw_measurements_[0];
      init_y = measurement_pack.raw_measurements_[1];
    }


    if(( 0==init_x)  && (0 == init_y) )
    {
      //invalid data, discard
      is_initialized_ = false;

    }
    else
    {
      ekf_.x_  = VectorXd(4);
      ekf_.x_ << init_x, init_y, 0 ,0;


      double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
      previous_timestamp_ = measurement_pack.timestamp_;
        
      ekf_.F_(0,2) = dt; 
      ekf_.F_(1,3) = dt; 
      // done initializing, no need to predict or update
      is_initialized_ = true;      
    }

    return;
  }




  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  //1. Modify the F matrix so that the time is integrated
  //compute the time elapsed between the current and previous measurements
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
    
  ekf_.F_(0,2) = dt; 
  ekf_.F_(1,3) = dt; 
  
  //2. Set the process covariance matrix Q
  ekf_.Q_= MatrixXd(4, 4);
    
  ekf_.Q_ <<  std::pow(dt,4)*noise_ax/4,      0,                          std::pow(dt,3)*noise_ax/2,      0,
              0,                              std::pow(dt,4)*noise_ay/4,  0,                              std::pow(dt,3)*noise_ay/2,        
              std::pow(dt,3)*noise_ax/2,      0,                          std::pow(dt,2)*noise_ax,        0,  
              0,                              std::pow(dt,3)*noise_ay/2,  0,                              std::pow(dt,2)*noise_ay;  
  

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }


}

