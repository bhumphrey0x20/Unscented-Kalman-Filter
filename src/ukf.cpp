#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

int counter = 0; 

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  ///* initially set to false, set to true in first call of ProcessMeasurement
	is_initialized_ = false; 

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.2; //1.788; //30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.788; //3.1975; //0.25; //30;
  
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
  
	time_us_ = 0;
	n_x_ = 5; 
	n_aug_ = 7;

	lambda_ = 3 - n_aug_;

	weights_ = VectorXd(2 * n_aug_ + 1); 
	weights_(0) = lambda_/(n_aug_ + lambda_); 
	for(int i = 1; i< (2*n_aug_+1); i++){
		weights_(i) = 0.5/(n_aug_ + lambda_); 
	}


	//create matrix with predicted sigma points as columns
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	// flag to write NIS values to .dat file
	WRITE_NIS = false; 
	NIS = 0.0; 

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	//static int cnt = 0;  // cnt input measurements, for debugging

	if( !is_initialized_){
	// First measurement: init x_ with measurement

	

		if(meas_package.sensor_type_ == MeasurementPackage::LASER){
			double px = meas_package.raw_measurements_(0);
			double py = meas_package.raw_measurements_(1);
			double v = sqrt(px*px + py*py);

			x_ << px, py, v, 0, 0;
			
			/*** Initial P_ values taken from:
					https://github.com/mvirgo/Unscented-Kalman-Filter/blob/master/src/ukf.cpp  ***/

			P_ << std_radr_ * std_radr_, 0, 0, 0, 0,
						0, std_radr_ * std_radr_, 0, 0, 0,
						0, 0, 1, 0, 0,
						0, 0, 0, std_radphi_, 0,
						0, 0, 0, 0, std_radphi_;

			is_initialized_ = true;

		}else{
			VectorXd z_in = meas_package.raw_measurements_; //.head(3);
			update_x_polar_to_Cart(meas_package.raw_measurements_); //z_in);


			/*** Initial P_ values taken from:
					https://github.com/mvirgo/Unscented-Kalman-Filter/blob/master/src/ukf.cpp  ***/
		
      P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
            0, std_laspy_*std_laspy_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

	is_initialized_ = true;
		}
	 
	}else{
	// measurement is not first, continue with prediction and update
		
		double delta_t = (double)(meas_package.timestamp_ - time_us_)/1000000; 

		Prediction(delta_t); 	

		if(meas_package.sensor_type_ == MeasurementPackage::LASER){
			// Predict Lidar measurements and update state vector
			UpdateLidar(meas_package); 
			
		}else{
			// Predict radar measurements and update state vector
			UpdateRadar(meas_package);
					
		}
	}
	// set time_us_ as previous timestamp
	time_us_ = meas_package.timestamp_; 
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	
	// Generate Sigma Points
	lambda_ = 3-n_x_;
	
	MatrixXd Xsig_ = MatrixXd(n_x_, 2*n_x_+1); 
//	MatrixXd P_temp = lambda_ * P_; 
//	MatrixXd A = P_temp.llt().matrixL();
	MatrixXd A = P_.llt().matrixL();

	Xsig_.col(0) = x_;
	for( int col = 0; col < n_x_; col++){
		Xsig_.col(col+1) = x_ + sqrt(lambda_ + n_x_)*A.col(col);
		Xsig_.col(col+1+n_x_) = x_ - sqrt(lambda_ + n_x_)*A.col(col);
	}

	//*** Generate Augmented Sigma Points***
	//change lambda to match aug sig point matrix dims
	lambda_ = 3 - n_aug_;


	//create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
	x_aug.fill(0.0); 

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_); 
	P_aug.fill(0.0); 

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	Xsig_aug.fill(0.0); 

	x_aug.head(n_x_) = x_; 
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(5,5) = (std_a_ * std_a_); 
	P_aug(6,6) = (std_yawdd_ * std_yawdd_); 
		
  //create square root matrix
//	P_aug = (lambda_ + n_aug_) * P_aug;

	A = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for(int col = 0; col < n_aug_; col++){
//		Xsig_aug.col(col + 1) = x_aug + A.col(col);
//		Xsig_aug.col(col +1 + n_aug_) = x_aug - A.col(col);
		Xsig_aug.col(col + 1) = x_aug + sqrt(lambda_+n_aug_)*A.col(col);
		Xsig_aug.col(col +1 + n_aug_) = x_aug - sqrt(lambda_+n_aug_)*A.col(col);
	
	}


	//*** Predict Sigma Points ***

	//create matrix with predicted sigma points as columns
	// and create temp matrices to hold kinematic diff eq'n

	VectorXd x_temp = VectorXd(n_aug_); 
	VectorXd x = VectorXd(n_x_); 
	VectorXd x_int = VectorXd(n_x_); 
	VectorXd x_acc = VectorXd(n_x_); 

	Xsig_pred_.fill(0.0); 

	double coeff = 0; 
	double px= 0;
	double py = 0;
	double v = 0;
	double psi = 0; 
	double psi_dot = 0;
	double dt_2 = (delta_t * delta_t); 
	double nu_a = 0; 
	double nu_psi_dbl_dot = 0; 

	for(int i= 0; i< Xsig_aug.cols(); i++){
		x_temp 	= Xsig_aug.col(i); 
		px 			= x_temp[0];
		py 			= x_temp[1];
		v  			= x_temp[2];
		psi 		= x_temp[3];
		psi_dot = x_temp[4];
		nu_a 		= x_temp[5]; 
		nu_psi_dbl_dot = x_temp[6];

		//set x_int
		x << px, py, v, psi, psi_dot; 

		// Check for divide-by-zero
		if( psi_dot == 0){
			x_int[0] = v * cos(psi) * delta_t; 
			x_int[1] = v * sin(psi) * delta_t; 
			x_int[2] = 0; 
			x_int[3] = 0; //psi_dot * delta_t; 
			x_int[4] = 0; 
		}else{
			coeff = v/psi_dot; 
			x_int[0] = coeff * (sin(psi + (psi_dot * delta_t)) - sin(psi)); 
			x_int[1] = coeff * (-1 * cos(psi + (psi_dot * delta_t)) + cos(psi)); 
			x_int[2] = 0; 
			x_int[3] = psi_dot * delta_t;
			x_int[4] = 0; 
		}

	// calc acc component of X_k+1
			x_acc[0] = 0.5 * dt_2 * cos(psi) * nu_a; 
			x_acc[1] = 0.5 * dt_2 * sin(psi) * nu_a; 
			x_acc[2] = delta_t * nu_a; 	
			x_acc[3] = 0.5 * dt_2 * nu_psi_dbl_dot;
			x_acc[4] = delta_t * nu_psi_dbl_dot; 

		Xsig_pred_.col(i) = x + x_int + x_acc; 
	}//end for


	//*** Predict Mean and Covariance ***

	// Create Vector for predicted mean and Matr for predicted covariance
	x_mean = VectorXd(n_x_);
	x_mean.fill(0.0);
	P_pred = MatrixXd(n_x_, n_x_);
	P_pred.fill(0.0); 


  //predict state mean
	for(int col = 0; col < (2*n_aug_+1); col++){
		x_mean = x_mean + ( weights_(col) * Xsig_pred_.col(col) );  
	}

  //predict state covariance matrix

	for (int i = 0; i < (2*n_aug_+1); i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_mean;

    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_pred = P_pred + weights_(i) * x_diff * x_diff.transpose() ;
	}
	x_= x_mean;
	P_ = P_pred; 

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

	// *** Predict Lidar Measurement ***

	VectorXd z_diff; 

	int n_z_lidar = 2; 
	
	MatrixXd R_lidar = MatrixXd(2,2);
	R_lidar << std_radr_* std_radr_, 0,
        			0, std_radr_* std_radr_;

	//define spreading parameter
  lambda_ = 3 - n_aug_;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig_meas_ = MatrixXd(n_z_lidar, 2 * n_aug_ + 1);
	Zsig_meas_.fill(0.0); 

  //mean predicted measurement
  VectorXd z_mean_meas_ = VectorXd(n_z_lidar);
  z_mean_meas_.fill(0.0);

  //measurement covariance matrix S
  S = MatrixXd(n_z_lidar,n_z_lidar);
	S.fill(0.0);

  //transform sigma points into Lidar measurement space
	VectorXd z_temp  = Xsig_pred_.col(0); 
	Zsig_meas_ = Xsig_pred_.topRows(n_z_lidar); 


  //calculate mean predicted measurement
	for(int col = 0; col < (2*n_aug_ + 1); col++){
		z_mean_meas_ = z_mean_meas_ + ( weights_(col) * Zsig_meas_.col(col) );  
	}

	//Calculate Covar matrix S_
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    z_diff = Zsig_meas_.col(i) - z_mean_meas_;

    S = S + weights_(i) * z_diff * z_diff.transpose() ;
	}

	S = S + R_lidar; 

	//*** Update State  ***

	VectorXd z_in = meas_package.raw_measurements_.head(n_z_lidar);

//create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_lidar);
	Tc.fill(0.0); 

  //calculate cross correlation matrix
	VectorXd x_pred_diff = VectorXd(n_x_); 
	VectorXd z_pred_diff = VectorXd(n_z_lidar);

	for(int col = 0; col < (2*n_aug_ +1); col++){
		x_pred_diff = Xsig_pred_.col(col) - x_;
		while(x_pred_diff(3) > M_PI) x_pred_diff(3) -= (2.0 * M_PI);
		while(x_pred_diff(3) < -M_PI) x_pred_diff(3) += (2.0 * M_PI);

		z_pred_diff = Zsig_meas_.col(col) - z_mean_meas_;
//		while(z_pred_diff(1) > M_PI) z_pred_diff(1) -= (2.0 * M_PI);
//		while(z_pred_diff(1) < M_PI) z_pred_diff(1) += (2.0 * M_PI);

		Tc = Tc + (weights_(col) * x_pred_diff * z_pred_diff.transpose()); 	
	}

  //calculate Kalman gain K;
	MatrixXd Si = S.inverse();
	MatrixXd K = Tc * Si; 

  //update state mean and covariance matrix
	MatrixXd error = z_in - z_mean_meas_;
	x_ = x_ + K*error; 

	P_ = P_ - (K * S * K.transpose());

		VectorXd nis_temp = error.transpose() * Si * error; 

	NIS = nis_temp(0); 

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
 
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

	VectorXd z_in = meas_package.raw_measurements_.head(3); 

  //define spreading parameter
	lambda_ = 3 - n_aug_;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	Zsig.fill(0.0);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  //measurement covariance matrix S
  S = MatrixXd(n_z,n_z);
	S.fill(0.0);

	// Radar Measurement Noise R
	MatrixXd R = MatrixXd(n_z, n_z); 
	R << (std_radr_ * std_radr_), 0, 0,
								0, (std_radphi_ * std_radphi_), 0,
								0, 0, (std_radrd_ * std_radrd_); 


  //transform sigma points into measurement space
	//	counter++;

  //calculate mean predicted measurement
	for(int col = 0; col < (2*n_aug_+1); col++){
		VectorXd Xsig_vect = Xsig_pred_.col(col);

		double px = Xsig_vect(0);
		double py = Xsig_vect(1);
		double v 	= Xsig_vect(2);
		double psi = Xsig_vect(3);
		double rho = sqrt((px* px) + (py * py));
		double phi = atan2(py, px);
		double rho_dot = ( px * cos(psi) * v) + (py * sin(psi) * v);
		rho_dot = rho_dot/rho;
		Zsig.col(col) << rho, phi, rho_dot;
		z_pred = z_pred + ( weights_(col) * Zsig.col(col) );  

//		z_pred = z_pred + ( weights_(col) * Zsig.col(col) );  
	}


  //calculate innovation covariance matrix S	
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  

    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
//cout << "z_diff" << z_diff << endl << endl;
    //angle normalization

    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose() ;
	}

	S = S+R; 


	//*** Update State ***


  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);

  //calculate cross correlation matrix
	VectorXd X_temp = VectorXd(n_x_); 
	VectorXd Z_temp = VectorXd(n_z);
	  
	for(int col = 0; col < (2*n_aug_ +1); col++){
		X_temp = Xsig_pred_.col(col) - x_mean;
		while(X_temp(3) > M_PI) X_temp(3) -= 2.0 * M_PI;
		while(X_temp(3) < -M_PI) X_temp(3) += 2.0* M_PI;


		Z_temp = Zsig.col(col) - z_pred;
	/*** Normalizing Z_temp(1) values taken from:
			https://github.com/mvirgo/Unscented-Kalman-Filter/blob/master/src/ukf.cpp  ***/

		while(Z_temp(1) > M_PI) Z_temp(1)-= 2. * M_PI;
		while(Z_temp(1) <-M_PI) Z_temp(1)+= 2. * M_PI;

		Tc = Tc + (weights_(col) * X_temp * Z_temp.transpose()); 	
 	}

  //calculate Kalman gain K;
	MatrixXd Si = S.inverse();
	MatrixXd K = Tc * Si; 

  //update state mean and covariance matrix
	MatrixXd error = z_in-z_pred;

	x_ = x_ + K*error; 

	P_ = P_ - (K * S * K.transpose());

	VectorXd nis_temp = error.transpose() * Si * error; 

	NIS = nis_temp(0); 
}

void UKF::update_x_polar_to_Cart(const Eigen::VectorXd z_in){

		float px, py; 
		float rho = z_in[0];
		float phi = z_in[1];
		float rho_dot = z_in[2];
		float cos_p = cos(phi); 
		float sin_p = sin(phi); 

		px = rho * cos_p; 
		py = rho * sin_p; 

		x_ << px, py, rho_dot, rho_dot * cos_p, rho_dot * sin_p; 
}


