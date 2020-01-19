#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;
  cout << "starting initialize Matrises" << endl;
  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;
  cout << "R_laser" << endl;
  cout << R_laser_ << endl;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  cout << "R_radar" << endl;
  cout << R_radar_ << endl;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  cout << "H_laser" << endl;
  cout << H_laser_ << endl;

  Hj_ << 1, 1, 0, 0,
         1, 1, 0, 0,
         1, 1, 1, 1; 
  cout << "Hj_" << endl;
  cout << Hj_ << endl;
  
  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;
  cout << "ekf_.F_" << endl;
  cout << ekf_.F_ << endl;

  //state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
  cout << "edf_.P_" << endl;
  cout << ekf_.P_ << endl;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "first measurement" << endl;
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    // ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      float ro     = measurement_pack.raw_measurements_(0);
      float phi    = measurement_pack.raw_measurements_(1);
      float ro_dot = measurement_pack.raw_measurements_(2);
      ekf_.x_(0) = ro     * cos(phi);
      ekf_.x_(1) = ro     * sin(phi);      
      ekf_.x_(2) = ro_dot * cos(phi);
      ekf_.x_(3) = ro_dot * sin(phi);
      cout << "radar measurement" << endl;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.

      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
      ekf_.x_(2) = 0;
      ekf_.x_(3) = 0;
      cout << "lidar measurement" << endl;
    }
    
    cout << "ekf_.x_" << endl;  
    cout << ekf_.x_ << endl;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  cout << "prediction step" << endl;
    //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt   * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  cout << "ekf_.F_" << endl;
  cout << ekf_.F_ << endl;
  //set the acceleration noise components
  float noise_ax = 9;
  float noise_ay = 9;

  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
             0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
             dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;
  cout << "ekf_.Q_" << endl;
  cout << ekf_.Q_ << endl;
  ekf_.Predict();


  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */
  cout << "update step" << endl;
  cout << "ekf_.x_ after predict step" << endl;
  cout << ekf_.x_ << endl;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    Tools tools;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    cout << "H_ in radar update" << endl;
    cout << ekf_.H_ << endl;
    ekf_.R_ = R_radar_;
    cout << "R_ in radar update" << endl;
    cout << ekf_.R_ << endl;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    cout << "H_ in laser update" << endl;
    cout << ekf_.H_ << endl;
    cout << "R_ in laser update" << endl;
    cout << ekf_.R_ << endl;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
