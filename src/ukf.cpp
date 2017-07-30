#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // if this is true, sigma points are also used for lidar update
  use_sigma_points_for_lidar_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.75;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = n_x_ + 2;

  // Number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(n_sig_);
  // Weights of sigma points
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  // all other weights
  weights_(1) = 0.5 / (lambda_ + n_aug_);
  for (int i = 2; i < n_sig_; i++) {
    weights_(i) = weights_(1);
  }

  // time when the state is true, in us
  time_us_ = 0;

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0,                       std_laspy_ * std_laspy_;


  // measurement covariance matrix - laser
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_ * std_radr_, 0,                         0,
              0,                     std_radphi_ * std_radphi_, 0,
              0,                     0,                         std_radrd_ * std_radrd_;

  // measurement matrix
  H_laser_ = MatrixXd(2, 5);
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;

  // NIS metric initialization
  nis_laser_ = 0;
  nis_radar_ = 0;

  nis_laser_95_cnt_ = 0;
  nis_radar_95_cnt_ = 0;

  laser_cnt_ = 0;
  radar_cnt_ = 0;
}

///* chi squared critical values at 0.95 level for 2 and 3 degrees of freedom
const double UKF::CHI2_2DF_ = 5.991;
const double UKF::CHI2_3DF_ = 7.815;

UKF::~UKF() {}

/**
 *  Normalize angle to be between -pi and pi
 */
void UKF::NormalizeAngle(double &angle) {
  while (angle > M_PI) {
    angle -= 2. * M_PI;
  }
  while (angle < -M_PI) {
    angle += 2. * M_PI;
  }
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    cout << "Kalman Filter Initialization " << endl;
    // first measurement
    // NOTE: covariance matrices are initialized in the constructor
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // convert radar data from polar to cartesian coordinates and initialize state with zero velocity
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      x_ << rho * cos(phi), rho * sin(phi), 0, 0, 0;
    } else {  // LASER
      // set the state with the initial location and zero velocity
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    // set initial timestamp
    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  if ((use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
      || (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)) {
    // compute the time elapsed between the current and previous measurements
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;  // dt - expressed in seconds
    time_us_ = meas_package.timestamp_;

    // predict
    MatrixXd Xsig_pred = Prediction(dt);

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar measurement update
      UpdateRadar(meas_package, Xsig_pred);
      radar_cnt_++;
    } else {
      // Laser measurement update
      if (use_sigma_points_for_lidar_) {
        UpdateLidar(meas_package, Xsig_pred);
      } else {
        UpdateLidar(meas_package);
      }
      laser_cnt_++;
    }

    // print the output

    if (laser_cnt_ > 0) {
      double nis_percent_laser = double(nis_laser_95_cnt_) / double(laser_cnt_);
      cout << "\nNIS under 95% Laser = " << nis_percent_laser << endl;
    }
    if (radar_cnt_ > 0) {
      double nis_percent_radar = double(nis_radar_95_cnt_) / double(radar_cnt_);
      cout << "NIS under 95% Radar = " << nis_percent_radar << endl;
    }
    double nis_percent_all = (double(nis_laser_95_cnt_) + double(nis_radar_95_cnt_))  / (double(laser_cnt_) + double(radar_cnt_));
    cout << "NIS under 95% All = " << nis_percent_all << endl;

    cout << "\nx_ = " << x_ << endl;
    cout << "P_ = " << P_ << endl;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
MatrixXd UKF::Prediction(double delta_t) {
  // AUGMENT STATE AND GENERATE SIGMA POINTS

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug.tail(2) << 0, 0;
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  double scale = sqrt(lambda_ + n_aug_);

  Xsig_aug.col(0) = x_aug;
  for (unsigned int i = 0; i < n_aug_; i++) {
      Xsig_aug.col(i + 1)          = x_aug + scale * L.col(i);
      Xsig_aug.col(n_aug_ + i + 1) = x_aug - scale * L.col(i);
  }

  // PREDICTION OF SIGMA POINTS

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, n_sig_);

  VectorXd x_new = VectorXd(n_x_);
  //predict sigma points
  for (unsigned int i = 0; i < Xsig_aug.cols(); i++) {
    x_new = Xsig_aug.col(i).head(n_x_);
    double px = x_new(0);
    double py = x_new(1);
    double v = x_new(2);
    double phi = x_new(3);
    double phi_dot = x_new(4);
    double noise_a = Xsig_aug(5, i);
    double noise_phi = Xsig_aug(6, i);

    double sin_phi = sin(phi);
    double cos_phi = cos(phi);
    double delta_t2 = delta_t * delta_t;

    // transition
    if (fabs(phi_dot) < 0.001) {
      x_new(0) += v * cos_phi * delta_t;
      x_new(1) += v * sin_phi * delta_t;
    } else {
      x_new(0) += v / phi_dot * (sin(phi + phi_dot * delta_t) - sin_phi);
      x_new(1) += v / phi_dot * (-cos(phi + phi_dot * delta_t) + cos_phi);
    }
    x_new(3) += phi_dot * delta_t;

    // process noise
    x_new(0) += 0.5 * delta_t2 * cos_phi * noise_a;
    x_new(1) += 0.5 * delta_t2 * sin_phi * noise_a;
    x_new(2) += delta_t * noise_a;

    x_new(3) += 0.5 * delta_t2 * noise_phi;
    x_new(4) += delta_t * noise_phi;

    // write vector
    Xsig_pred.col(i) = x_new;
  }

  // PREDICT MEAN AND COVARIANCE

  // predict state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    x_ += weights_(i) * Xsig_pred.col(i);
  }

  //predict state covariance matrix
  P_.fill(0.0);
  for (unsigned int i = 0; i < n_sig_; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;

    //angle normalization
    NormalizeAngle(x_diff(3));

    P_ += weights_(i) * (x_diff * x_diff.transpose());
  }

  return Xsig_pred;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  // measurement vector
  VectorXd z = meas_package.raw_measurements_;

  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;

  nis_laser_ = y.transpose() * Si * y;
  if (nis_laser_ < CHI2_2DF_) {
    nis_laser_95_cnt_++;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement and sigma points
 * @param meas_package The measurement at k+1
 * @param Xsig_pred Sigma points
 */
void UKF::UpdateLidar(MeasurementPackage meas_package, const MatrixXd &Xsig_pred) {
  // TRANSFORM SIGMA POINTS INTO MEASUREMENT SPACE

  // measurement dimension
  int n_z = 2;
  MatrixXd Zsig = Xsig_pred.topLeftCorner(n_z, n_sig_);

  // UPDATE STEP USING SIGMA POINTS
  UpdateWithSigmaPoints(meas_package, n_z, Zsig, Xsig_pred);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement
 * @param meas_package The measurement at k+1
 * @param Xsig_pred Sigma points
 */
void UKF::UpdateRadar(MeasurementPackage meas_package, const MatrixXd &Xsig_pred) {
  // TRANSFORM SIGMA POINTS INTO MEASUREMENT SPACE

  // measurement dimension
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  for (unsigned int i = 0; i < n_sig_; i++) {
    double px  = Xsig_pred(0,i);
    double py  = Xsig_pred(1,i);
    double v   = Xsig_pred(2,i);
    double phi = Xsig_pred(3,i);

    Zsig(0, i) = sqrt(px * px + py * py);
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px * cos(phi) * v + py * sin(phi) * v) / Zsig(0, i);
  }

  // UPDATE STEP USING SIGMA POINTS
  UpdateWithSigmaPoints(meas_package, n_z, Zsig, Xsig_pred);
}

/**
 * Sigma points update subroutine
 * @param meas_package The measurement at k+1
 * @param n_z Dimension of measurement
 * @param Zsig Sigma points in measurement space
 * @param Xsig_pred Sigma points in state space
 */
void UKF::UpdateWithSigmaPoints(MeasurementPackage meas_package, const int &n_z, const MatrixXd &Zsig, const MatrixXd &Xsig_pred) {
  // measurement vector
  VectorXd z = meas_package.raw_measurements_;

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points
    z_pred += weights_(i) * Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points
    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      //angle normalization
      NormalizeAngle(z_diff(1));
    }

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    S += R_radar_;
  } else {
    S += R_laser_;
  }

  // UPDATE STATE AND COVARIANCE

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (unsigned int i = 0; i < n_sig_; i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      //angle normalization
      NormalizeAngle(z_diff(1));
    }

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    //angle normalization
    NormalizeAngle(x_diff(3));

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd Si = S.inverse();
  MatrixXd K = Tc * Si;

  //residual
  VectorXd z_diff = z - z_pred;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    //angle normalization
    NormalizeAngle(z_diff(1));
  }

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    nis_radar_ = z_diff.transpose() * Si * z_diff;
    if (nis_radar_ < CHI2_3DF_) {
      nis_radar_95_cnt_++;
    }
  } else {  // LASER
    nis_laser_ = z_diff.transpose() * Si * z_diff;
    if (nis_laser_ < CHI2_2DF_) {
      nis_laser_95_cnt_++;
    }
  }
}
