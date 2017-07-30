#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* if this is true, sigma points are also used for lidar update
  bool use_sigma_points_for_lidar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Number of sigma points
  int n_sig_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* measurement covariance matrices
  MatrixXd R_laser_;
  MatrixXd R_radar_;
  ///* transformation matrix
  MatrixXd H_laser_;

  ///* NIS measure
  double nis_laser_;
  double nis_radar_;

  int nis_laser_95_cnt_;
  int nis_radar_95_cnt_;

  int laser_cnt_;
  int radar_cnt_;

  ///* chi squared critical values at 0.95 level for 2 and 3 degrees of freedom
  static const double CHI2_2DF_;
  static const double CHI2_3DF_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   *  Normalize angle to be between -pi and pi
   */
  void NormalizeAngle(double &angle);

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Predicts sigma points, the state, and the state covariance matrix.
   * @param {double} delta_t the change in time (in seconds) between the last
   * measurement and this one.
   */
  MatrixXd Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a laser measurement and sigma points
   * @param meas_package The measurement at k+1
   * @param Xsig_pred Sigma points
   */
  void UpdateLidar(MeasurementPackage meas_package, const MatrixXd &Xsig_pred);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   * @param Xsig_pred Sigma points
   */
  void UpdateRadar(MeasurementPackage meas_package, const MatrixXd &Xsig_pred);

private:
  /**
   * Sigma points update subroutine
   * @param meas_package The measurement at k+1
   * @param n_z Dimension of measurement
   * @param Zsig Sigma points in measurement space
   * @param Xsig_pred Sigma points in state space
   */
  void UpdateWithSigmaPoints(MeasurementPackage meas_package, const int &n_z, const MatrixXd &Zsig, const MatrixXd &Xsig_pred);
};

#endif /* UKF_H */
