#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// state dimension
	n_x_ = 5;

	// initial state vector
	x_ = VectorXd(n_x_);

	// initial covariance matrix
	P_ = MatrixXd(n_x_, n_x_);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 2.0;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.8;

	// vector nu
	nu_ = VectorXd(2);
	nu_ << std_a_,
			std_yawdd_;

	// augmented state dimension
	n_aug_ = n_x_ + nu_.size();

	// number of sigma points
	n_sig_ = 2.0 * n_aug_ + 1;

	// predicted sigma points matrix
	Xsig_pred_ = MatrixXd(n_x_, n_sig_);

	// sigma point spread factor
	lambda_ = 3 - n_x_;

	// weight vector
	weights_ = VectorXd(n_sig_);

	// NIS
	NIS_radar_ = 0.0;
	NIS_laser_ = 0.0;

	//DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
	//DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

	/**
	 TODO:

	 Complete the initialization. See ukf.h for other member properties.

	 Hint: one or more values initialized above might be wildly off...
	 */
}

UKF::~UKF() {
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	/**
	 * TODO:
	 Complete this function! Make sure you switch between lidar and radar
	 measurements.
	 */
	if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
			|| (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)) {
		if (!is_initialized_) {
			Initialize_Process(meas_package);
		}

		/*****************************************************************************
		 *  Prediction
		 ****************************************************************************/
		float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
		time_us_ = meas_package.timestamp_;

		Prediction(dt);

		/*****************************************************************************
		 *  Update
		 ****************************************************************************/
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			// Radar updates
			UpdateRadar(meas_package);
		} else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			// Laser updates
			UpdateLidar(meas_package);
		}
	}
}

/**
 * Initializes the state and the state covariance matrix using current sensor type.
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::Initialize_Process(MeasurementPackage meas_package) {
	// init state
	x_.fill(0.0);

	// init covariance matrix
	P_ = MatrixXd::Identity(n_x_, n_x_);

	// init timestamp
	time_us_ = meas_package.timestamp_;

	if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
		/**
		 Initialize state from lidar measurements
		 */
		x_(0) = meas_package.raw_measurements_(0);
		x_(1) = meas_package.raw_measurements_(1);

		P_(0, 0) = std_laspx_ * std_laspx_;
		P_(1, 1) = std_laspy_ * std_laspy_;
	} else if (meas_package.sensor_type_ == MeasurementPackage::RADAR
			&& use_radar_) {
		/**
		 Convert radar from polar to cartesian coordinates and initialize state.
		 */
		float rho = meas_package.raw_measurements_(0);
		float phi = meas_package.raw_measurements_(1);
		x_(0) = rho * cos(phi);
		x_(1) = rho * sin(phi);

		//calculate variance for position measurement
		// assumption: tangential std_dev can be approximated with formula {1}
		double norm = sqrt(x_(0) * x_(0) + x_(1) * x_(1));
		double std_radxy = norm * sin(std_radphi_); // formula {1}
		if (std_radxy < std_radr_) {
			std_radxy = std_radr_;
		}
		std::cout << "std_radxy: " << std_radxy << std::endl;
		P_(0, 0) = std_radxy * std_radxy;
		P_(1, 1) = std_radxy * std_radxy;
	}

	// done initializing, no need to predict or update
	//std::cout << "..initialized!" << std::endl;
	is_initialized_ = true;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	/**
	 TODO:

	 Complete this function! Estimate the object's location. Modify the state
	 vector, x_. Predict sigma points, the state, and the state covariance matrix.
	 */

	//create sigma point matrix
	MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

	//calculate square root of P
	MatrixXd A = P_.llt().matrixL();

	//set first column of sigma point matrix
	Xsig.col(0) = x_;

	//set remaining sigma points
	for (int i = 0; i < n_x_; i++) {
		Xsig.col(i + 1) = x_ + sqrt(lambda_ + n_x_) * A.col(i);
		Xsig.col(i + 1 + n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
	}

	/*****************************************************************************
	 *  Augment Sigma Points
	 ****************************************************************************/

	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

	// sigma point spread factor for augmented
	lambda_ = 3 - n_aug_;

	//create augmented mean state
	x_aug.fill(0.0);
	x_aug.head(n_x_) = x_;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	for (int i = 0; i < nu_.size(); i++) {
		P_aug(n_x_ + i, n_x_ + i) = nu_(i) * nu_(i);
	}

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i < n_aug_; i++) {
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
	}

	/*****************************************************************************

	 *  Predict Sigma Points

	 ****************************************************************************/

	for (int i = 0; i < n_sig_; i++) {
		//extract values for better readability
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		//predicted state values
		VectorXd dx(n_x_);
		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			dx.head(2) << v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw)),
					v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
		} else {
			dx.head(2) << v * delta_t * cos(yaw),
					v * delta_t * sin(yaw);
		}
		// write common values
		dx.tail(n_x_ - 2) << 0,
				yawd * delta_t,
				0;

		//calculate noise vector
		VectorXd sig_noise(n_x_);
		sig_noise << 0.5 * nu_a * delta_t * delta_t * cos(yaw),
				0.5 * nu_a * delta_t * delta_t * sin(yaw),
				nu_a * delta_t,
				0.5 * nu_yawdd * delta_t * delta_t,
				nu_yawdd * delta_t;

		//write predicted sigma point into right column (x + dx + noise)
		Xsig_pred_.col(i) = Xsig_aug.col(i).head(n_x_) + dx + sig_noise;
	}

	/*****************************************************************************

	 *  Convert Predicted Sigma Points to Mean/Covariance

	 ****************************************************************************/

	// set weights
	double weight_0 = lambda_ / (lambda_ + n_aug_);
	weights_(0) = weight_0;
	for (int i = 1; i < n_sig_; i++) {  //2n+1 weights
		double weight = 0.5 / (n_aug_ + lambda_);
		weights_(i) = weight;
	}

	//predicted state mean
	x_.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points
		x_ += weights_(i) * Xsig_pred_.col(i);
	}

	//predicted state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3) > M_PI)
			x_diff(3) -= 2. * M_PI;
		while (x_diff(3) < -M_PI)
			x_diff(3) += 2. * M_PI;
		P_ += weights_(i) * x_diff * x_diff.transpose();
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	/**
	 TODO:

	 Complete this function! Use lidar data to update the belief about the object's
	 position. Modify the state vector, x_, and covariance, P_.

	 You'll also need to calculate the lidar NIS.
	 */

	//extract measurement as VectorXd
	VectorXd z = meas_package.raw_measurements_;

	//set measurement dimension, lidar can measure p_x and p_y
	int n_z = 2;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, n_sig_);

	//transform sigma points into measurement space
	for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points
		// extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);

		// measurement model
		Zsig(0, i) = p_x;
		Zsig(1, i) = p_y;
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		z_pred += weights_(i) * Zsig.col(i);
	}

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		S += weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_laspx_ * std_laspx_, 0,
			0, std_laspy_ * std_laspy_;

	S += R;

	/*****************************************************************************
	 *  UKF Update for Lidar
	 ****************************************************************************/
	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		Tc += weights_(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z_diff = z - z_pred;

	//calculate NIS
	NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
	std::cout << "NIS laser: " << NIS_laser_ << std::endl;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	/**
	 TODO:

	 Complete this function! Use radar data to update the belief about the object's
	 position. Modify the state vector, x_, and covariance, P_.

	 You'll also need to calculate the radar NIS.
	 */

	//extract measurement as VectorXd
	VectorXd z = meas_package.raw_measurements_;

	//set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 3;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, n_sig_);

	//transform sigma points into measurement space
	for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points
		// extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw) * v;
		double v2 = sin(yaw) * v;

		// measurement model
		Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);													//rho
		Zsig(1, i) = atan2(p_y, p_x);             												//phi
		Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y); //rho_dot
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		z_pred += weights_(i) * Zsig.col(i);
	}

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1) > M_PI)
			z_diff(1) -= 2. * M_PI;
		while (z_diff(1) < -M_PI)
			z_diff(1) += 2. * M_PI;

		S += weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_radr_ * std_radr_, 0, 0,
			0, std_radphi_ * std_radphi_, 0,
			0, 0, std_radrd_ * std_radrd_;

	S += R;

	/*****************************************************************************
	 *  UKF Update for Radar
	 ****************************************************************************/
	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1) > M_PI)
			z_diff(1) -= 2. * M_PI;
		while (z_diff(1) < -M_PI)
			z_diff(1) += 2. * M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		//angle normalization
		while (x_diff(3) > M_PI)
			x_diff(3) -= 2. * M_PI;
		while (x_diff(3) < -M_PI)
			x_diff(3) += 2. * M_PI;

		Tc += weights_(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z_diff = z - z_pred;

	//angle normalization
	while (z_diff(1) > M_PI)
		z_diff(1) -= 2. * M_PI;
	while (z_diff(1) < -M_PI)
		z_diff(1) += 2. * M_PI;

	//calculate NIS
	NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
	std::cout << "NIS radar: " << NIS_radar_ << std::endl;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();
}
