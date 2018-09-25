#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

	// Variables: 
	VectorXd rmse(4); 
	VectorXd temp(4);

	rmse.fill(0.0);
	temp.fill(0.0);


	/* Check validity of input:
		* estimation/ground truth size > 0
		* estimation == ground truth  */
	int est_size = estimations.size();
	int truth_size = ground_truth.size();

	if(est_size == 0){
		cout << "RMSE Error: esitmation size == 0!" << endl;
	}
	if(est_size != truth_size){
		cout << "RMSE Error: Estimation Size and Ground Truth Size Not Equal Length!" << endl;
	}

	// Accumulator
	for(int i=0; i< est_size; i++){
		temp = estimations[i] - ground_truth[i];
		temp = temp.array() * temp.array();
		rmse += temp;  
	}	

	// Calc mean then RMSE
	rmse = rmse/est_size; 
	rmse = rmse.array().sqrt();

	return rmse; 
}
