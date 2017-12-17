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
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  ekf_.x_ = VectorXd(4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  long long cur_timestamp = measurement_pack.timestamp_;

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2];

      // Get four measurements: px, py, vx, vy
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);

      // Save first measurement
      ekf_.x_ << px, py, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      // Get two measurements: px, py
      double px = measurement_pack.raw_measurements_[0];
      double py = measurement_pack.raw_measurements_[1];

      // Laser won't tell use vx, vy, so we use 0 as initial value
      ekf_.x_ << px, py, 0, 0;
    }

    // Set initial value for F. But I think this step is optional since we will reset F when we
    // do ext prediction
    ekf_.F_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 1;

    // State covariance matrix. All of our four state variable are independent of each other
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;

    previous_timestamp_ = cur_timestamp;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // Update state transition matrix
  double delta_t = (cur_timestamp - previous_timestamp_) / 1000000.0;
  ekf_.F_ << 1, 0, delta_t, 0,
             0, 1, 0, delta_t,
             0, 0, 1, 0,
             0, 0, 0, 1;

  // Update process noise covariance matrix (Q is not zero because we are considering ax and ay as
  // process noise
  const double noise_ax = 9;
  const double noise_ay = 9;
  double var_ax = noise_ax;
  double var_ay = noise_ay;
  double delta_t_2 = delta_t * delta_t;
  double delta_t_3 = delta_t * delta_t_2;
  double delta_t_4 = delta_t * delta_t_3;
  ekf_.Q_ << delta_t_4 / 4 * var_ax, 0, delta_t_3 / 2 * var_ax, 0,
             0, delta_t_4 / 4 * var_ay, 0, delta_t_3 / 2 * var_ay,
             delta_t_3 / 2 * var_ax, 0, delta_t_2 * var_ax, 0,
             0, delta_t_3 / 2 * var_ay, 0, delta_t_2 * var_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Update measurement
    VectorXd measurement(3);
    measurement << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1],
                   measurement_pack.raw_measurements_[2];

    // Update measurement noise covariance matrix
    ekf_.R_ = R_radar_;

    // Update measurement matrix
    ekf_.H_ = MatrixXd(3, 4);
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

    ekf_.UpdateEKF(measurement);
  } else { // Laser
    // Update measurement
    VectorXd measurement(2);
    measurement << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];

    // Update measurement covariance matrix
    ekf_.R_ = R_laser_;

    // Update measurement covariance matrix
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ << 1, 0, 0, 0,
               0, 1, 0, 0;

    ekf_.Update(measurement);
  }

  // Update timestamp
  previous_timestamp_ = cur_timestamp;

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
