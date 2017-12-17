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
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() != ground_truth.size() || ground_truth.empty()) {
    cout << "Invalid ground truth data" << endl;
  }

  for (auto i = 0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  rmse /= estimations.size();

  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // Predefine some terms to make the code look shorter
  float t1 = px * px + py * py;
  float t2 = sqrt(t1);
  float t3 = (t1 * t2);

  // Check division by zero
  if (fabs(t1) < 0.0001) {
    cout << "Divided by zero" << endl;
    return Hj;
  }

  Hj << (px / t2), (py / t2), 0, 0,
        -(py / t1), (px / t1), 0, 0,
        py * (vx * py - vy * px) / t3, px * (vy * px - vx * py) / t3, px / t2, py / t2;

  return Hj;
}

VectorXd Tools::PloarToCatesian(const VectorXd& polar) {

  return Eigen::VectorXd();
}
