#include "ukf.h"
#include "tools.h"
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
    
    is_initialized_ = false;
    
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;
    
    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;
    
    // initial state vector
    x_ = VectorXd(5);
    
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 0.4;
    
    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.6;
    
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
    
    //set state dimension
    n_x_ = 5;
    
    //set augmented dimension
    n_aug_ = n_x_ + 2;
    
    //define spreading parameter
    lambda_ = 3 - n_aug_;
    
    //initial covariance matrix
    P_ = MatrixXd(5, 5);
    
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;
    
    //create sigma point matrix
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    
    //weights
    weights_ = VectorXd(2*n_aug_+1);
    long double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
        long double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }
    
    previous_timestamp_ = 0;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	
    if (!is_initialized_) {
        
        long double px = 0;
        long double py = 0;
        long double v = 0;
        
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            
            //initialize rho, phi and calculate px,py
            long double rho = meas_package.raw_measurements_[0];
            long double phi = meas_package.raw_measurements_[1];
            long double rho_dot = meas_package.raw_measurements_[2];
            
            px = rho * cos(phi);
            py = rho * sin(phi);
            v = rho_dot;
            
            if(fabs(px) < 0.0001){
                px = 0.0001;
                
            }
            if(fabs(py) < 0.0001){
                py = 0.0001;
                
            }
            
            x_ << px, py, v, 0, 0;
            previous_timestamp_ = meas_package.timestamp_;
            is_initialized_ = true;
            
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            
            //initialize px,py from laser data
            px = meas_package.raw_measurements_[0];
            py = meas_package.raw_measurements_[1];
            
            if(fabs(px) < 0.0001){
                px = 0.0001;
                
            }
            if(fabs(py) < 0.0001){
                py = 0.0001;
                
            }
            
            x_ << px, py, 0, 0, 0;
            previous_timestamp_ = meas_package.timestamp_;
            is_initialized_ = true;
            
        }
        return;
        
    }
    
    //checks if laser and/or radar data is used and if the measurement is laser or radar
    if (!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
        return;
    
    if (!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
        return;
    
    long double current_timestamp = meas_package.timestamp_;
    long double deltat = (current_timestamp - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = meas_package.timestamp_;
    
    Prediction(deltat);
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        UpdateRadar(meas_package);
    }
    else {
        // Laser updates
        UpdateLidar(meas_package);
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */

void UKF::Prediction(double delta_t) {

  	//create augmented mean vector
  	VectorXd x_aug_ = VectorXd(7);

  	//create augmented state covariance
  	MatrixXd P_aug_ = MatrixXd(7, 7);

  	//create sigma point matrix
  	MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  	//create augmented mean state
  	x_aug_.head(5) = x_;
  	x_aug_(5) = 0;
  	x_aug_(6) = 0;

  	//create augmented covariance matrix
  	P_aug_.fill(0.0);
  	P_aug_.topLeftCorner(5,5) = P_;
  	P_aug_(5,5) = std_a_*std_a_;
  	P_aug_(6,6) = std_yawdd_*std_yawdd_;

  	//create square root matrix
  	MatrixXd L = P_aug_.llt().matrixL();

  	//create augmented sigma points
  	Xsig_aug_.col(0)  = x_aug_;
  	for (int i = 0; i< n_aug_; i++)
  	{
    	Xsig_aug_.col(i+1)        = x_aug_ + sqrt(lambda_+n_aug_) * L.col(i);
    	Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L.col(i);
  	}

  	//predict sigma points
  	for (int i = 0; i< 2*n_aug_+1; i++)
  	{
    	//extract values for better readability
    	long double p_x = Xsig_aug_(0,i);
    	long double p_y = Xsig_aug_(1,i);
    	long double v = Xsig_aug_(2,i);
    	long double yaw = Xsig_aug_(3,i);
    	long double yawd = Xsig_aug_(4,i);
    	long double nu_a = Xsig_aug_(5,i);
    	long double nu_yawdd = Xsig_aug_(6,i);

    	//predicted state values
    	long double px_p, py_p;

    	//avoid division by zero
    	if (fabs(yawd) > 0.001) {
        	px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        	py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    	}
    	else {
        	px_p = p_x + v*delta_t*cos(yaw);
        	py_p = p_y + v*delta_t*sin(yaw);
    	}

    	long double v_p = v;
    	long double yaw_p = yaw + yawd*delta_t;
    	long double yawd_p = yawd;

    	//add noise
    	px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    	py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    	v_p = v_p + nu_a*delta_t;

    	yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    	yawd_p = yawd_p + nu_yawdd*delta_t;

    	//write predicted sigma point into right column
    	Xsig_pred_(0,i) = px_p;
    	Xsig_pred_(1,i) = py_p;
    	Xsig_pred_(2,i) = v_p;
    	Xsig_pred_(3,i) = yaw_p;
    	Xsig_pred_(4,i) = yawd_p;
  	}

  	//create vector for predicted state
  	VectorXd x = VectorXd(n_x_);

  	//create covariance matrix for prediction
  	MatrixXd P = MatrixXd(n_x_, n_x_);

  	//predicted state mean
  	x.fill(0.0);
  	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    	x = x+ weights_(i) * Xsig_pred_.col(i);
  	}

  	//predicted state covariance matrix
  	P.fill(0.0);
  	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    	// state difference
    	VectorXd x_diff = Xsig_pred_.col(i) - x;

    	//angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    	P = P + weights_(i) * x_diff * x_diff.transpose() ;
  	}

  	x_ = x;
  	P_ = P;

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::UpdateLidar(MeasurementPackage meas_package) {

	//set measurement dimension, lidar can measure x, y
  	int n_z_laser = 2;

	//create matrix for sigma points in measurement space
  	MatrixXd Zsig_ = MatrixXd(n_z_laser, 2 * n_aug_ + 1);

  	//transform sigma points into measurement space
  	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    	// extract values for better readibility
    	long double p_x = Xsig_pred_(0,i);
    	long double p_y = Xsig_pred_(1,i);

    	// measurement model
    	Zsig_(0,i) = p_x;
    	Zsig_(1,i) = p_y;
    	
  	}

  	//mean predicted measurement
  	VectorXd z_pred = VectorXd(n_z_laser);
  	z_pred.fill(0.0);
  	for (int i=0; i < 2*n_aug_+1; i++) {
    	z_pred = z_pred + weights_(i) * Zsig_.col(i);
  	}

  	//measurement covariance matrix S
  	MatrixXd S = MatrixXd(n_z_laser,n_z_laser);
  	S.fill(0.0);
  	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    	//residual
    	VectorXd z_diff = Zsig_.col(i) - z_pred;

    	//angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    	S = S + weights_(i) * z_diff * z_diff.transpose();
  	}

  	//add measurement noise covariance matrix
  	MatrixXd R_laser = MatrixXd(n_z_laser,n_z_laser);
  	
    R_laser << std_laspx_*std_laspx_, 0,
          	   0, std_laspy_*std_laspy_;
          
  	S = S + R_laser;

  	//create matrix for cross correlation Tc
  	MatrixXd Tc = MatrixXd(n_x_, n_z_laser);

  	//measurement matrix
    VectorXd z = meas_package.raw_measurements_;

  	//calculate cross correlation matrix
  	Tc.fill(0.0);
  	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    	//residual
    	VectorXd z_diff = Zsig_.col(i) - z_pred;
    	//angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    	// state difference
    	VectorXd x_diff = Xsig_pred_.col(i) - x_;
    	//angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    	Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  	}

  	//Kalman gain K;
  	MatrixXd K = Tc * S.inverse();

  	//residual
  	VectorXd z_diff = z - z_pred;

  	//angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  	//update state mean and covariance matrix
  	x_ = x_ + K * z_diff;
  	P_ = P_ - K*S*K.transpose();

    // calculate NIS Laser
  	NIS_laser_ = z_diff.transpose()*S.inverse()*z_diff;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::UpdateRadar(MeasurementPackage meas_package) {
 
	//set measurement dimension, radar can measure r, phi, and r_dot
  	int n_z_radar = 3;

	//create matrix for sigma points in measurement space
  	MatrixXd Zsig_ = MatrixXd(n_z_radar, 2 * n_aug_ + 1);

  	//transform sigma points into measurement space
  	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    	// extract values for better readibility
    	long double p_x = Xsig_pred_(0,i);
    	long double p_y = Xsig_pred_(1,i);
    	long double v  = Xsig_pred_(2,i);
    	long double yaw = Xsig_pred_(3,i);

    	long double v1 = cos(yaw)*v;
    	long double v2 = sin(yaw)*v;

    	// measurement model
    	Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    	Zsig_(1,i) = atan2(p_y,p_x);                                 //phi
    	Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  	}

  	//mean predicted measurement
  	VectorXd z_pred = VectorXd(n_z_radar);
  	z_pred.fill(0.0);
  	for (int i=0; i < 2*n_aug_+1; i++) {
    	z_pred = z_pred + weights_(i) * Zsig_.col(i);
  	}

  	//measurement covariance matrix S
  	MatrixXd S = MatrixXd(n_z_radar,n_z_radar);
  	S.fill(0.0);
  	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    	//residual
    	VectorXd z_diff = Zsig_.col(i) - z_pred;
    	
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    	S = S + weights_(i) * z_diff * z_diff.transpose();
  	}

  	//add measurement noise covariance matrix
  	MatrixXd R_radar = MatrixXd(n_z_radar,n_z_radar);
  	
    R_radar << std_radr_*std_radr_, 0, 0,
               0, std_radphi_*std_radphi_, 0,
               0, 0,std_radrd_*std_radrd_;
  	
    S = S + R_radar;

  	//create matrix for cross correlation Tc
  	MatrixXd Tc = MatrixXd(n_x_, n_z_radar);

  	//measurement matrix
    VectorXd z = meas_package.raw_measurements_;

  	//calculate cross correlation matrix
  	Tc.fill(0.0);
  	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    	//residual
    	VectorXd z_diff = Zsig_.col(i) - z_pred;
    	//angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    	// state difference
    	VectorXd x_diff = Xsig_pred_.col(i) - x_;
    	//angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    	Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  	}

  	//Kalman gain K;
  	MatrixXd K = Tc * S.inverse();

  	//residual
  	VectorXd z_diff = z - z_pred;

  	//angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  	//update state mean and covariance matrix
  	x_ = x_ + K * z_diff;
  	P_ = P_ - K*S*K.transpose();

    // calculate NIS Radar
  	NIS_radar_ = z_diff.transpose()*S.inverse()*z_diff;

}
