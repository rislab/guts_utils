#include <guts_utils/GUTSUtils.h>

#define LAMBDA 0.01 // from Equation 4.5 of Nikhil Bakshi's MSR thesis
#define NOISE 0.008 // Noise value

namespace gu=guts_utils;
using Point=gu::Grid::Point;

gu::GUTS::GUTS()
{
  EMitr_ = 1;
  a = 0.1;
  b = 1;
  eps_ = 1e-08;
}

gu::GUTS::~GUTS()
{
}

gu::Matrix gu::GUTS::multivariateNormal(const gu::Vector& mean,
					const gu::Matrix& cov)
{
  StandardNormal standard_normal;
  gu::Vector x = standard_normal.sample(mean.rows());
  Eigen::BDCSVD<gu::Matrix> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
  gu::Matrix U = svd.matrixU();
  gu::Matrix S = svd.singularValues();

  // See https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
  // for information about how to sample from a multivariate Gaussian distribution
  return U * S.cwiseSqrt().asDiagonal() * x + mean;
}


void gu::GUTS::posterior(const gu::Matrix& X, const gu::Matrix& Y)
{
  T a = 0.1;
  unsigned int b = 1;
  unsigned int n = X.cols();

  Vector beta_hat = Vector::Zero(n, 1);
  Vector gamma = Vector::Constant(n, 1, 1e8);
  Matrix Sig_beta = Matrix::Zero(n, n);
  Matrix ginv_xtsx = Matrix::Zero(n, n);
  Matrix Gamma = Matrix::Zero(n, n);
  Matrix inv_y_one;

  std::vector<bool> idx(n);
  for (unsigned int i = 0; i < n; ++i)
    idx[i] = false;

  for (unsigned int i = 0; i < EMitr_; ++i)
  {
    Gamma = gamma.asDiagonal();
    inv_y_one = Y.col(1).cwiseInverse().asDiagonal();

    Matrix Gamma_diag_inv = Gamma.diagonal().cwiseInverse().asDiagonal();

    // Equation (3), X^T \Sigma X + \Gamma^{-1}
    ginv_xtsx = X.transpose() * inv_y_one * X + Gamma_diag_inv;

    // Equation (3), V == Sig_beta
    Sig_beta = ginv_xtsx.diagonal().cwiseInverse().asDiagonal();

    Vector beta_hat_tmp  = (Sig_beta * X.transpose())*(inv_y_one * Y.col(0));

    // Equation (4)
    for (unsigned int j = 0; j < n; ++j)
    {
      gamma(j,0) = (2*b + beta_hat_tmp(j,0)*beta_hat_tmp(j,0) +
		    Sig_beta(j,j)) / (1 + 2 * a);

      if (idx[j] == false)
	beta_hat[j, 0] = beta_hat_tmp[j];

      if (gamma(j, 0) < eps_)
      {
	gamma(j, 0) = eps_;
	beta_hat(j,0) = 0.0;
	idx[j] = true;
      }
    }
  }

  Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(Sig_beta);
  if (eigensolver.info() == Eigen::NumericalIssue)
  {
    std::cerr << "WARNING! Sig_beta is not PSD" << std::endl;
    makeCovariancePSD(Sig_beta);
    eigensolver = Eigen::SelfAdjointEigenSolver<Matrix>(Sig_beta);
  }
  T min_eig = eigensolver.eigenvalues().minCoeff();

  // This is a modification made in Nikhil Bakshi's MSR thesis because
  // beta_tilde is the posterior, which is a probability.
  // Probabilities should be between zero and one; however, note that
  // the values may be larger than 1.0.
  Vector beta_tilde = multivariateNormal(beta_hat, Sig_beta);
  for (int i = 0; i < beta_tilde.rows(); ++i)
    if (beta_tilde(i) < 0.0)
      beta_tilde(i) = 0.0;

  ///// Debugging
  beta_tilde_ = beta_tilde;
  Gamma_ = Gamma;
  Sig_beta_ = Sig_beta;
  beta_hat_ = beta_hat;
  ginv_xtsx_ = ginv_xtsx;
  inv_y_one_ = Y.col(1).cwiseInverse().asDiagonal();
  gamma_ = gamma;
  min_eig_ = min_eig;
  X_ = X;
  Y_ = Y;
  ///// End debugging
}


void gu::GUTS::makeCovariancePSD(gu::Matrix& C)
{
  C = 0.5 * ( C + C.transpose() );
}


Point toGridPoint(const gu::Vector2& p)
{
  return p.cast<float>();
}


bool gu::GUTS::createDirectionalSensor(const gu::Vector2& curr_pos, const gu::Vector2& next_pos,
				       const Grid& map, gu::Matrix& x, gu::Vector& noise_var)
{
  Point cp = toGridPoint(curr_pos);
  Point np = toGridPoint(next_pos);

  if (map.inGrid(cp) && map.inGrid(np))
  {
    std::pair<bool, std::vector<Cell>> ret = map.raytrace(cp, np);
    bool success = ret.first;
    std::vector<Cell> cells = ret.second;

    if (success == false)
    {
      std::cerr << "GUTS Failed to raytrace" << std::endl;
      return false;
    }

    x = Matrix::Zero(cells.size(), map.data.size());
    for (int i = 0; i < cells.size(); ++i)
    {
      unsigned int idx = map.cellToIndex(cells[i]);
      x(i, idx) = 1;
    }

    noise_var = Vector::Constant(cells.size(), 1, NOISE);

    return true;
  }
  else
  {
    if (!map.inGrid(cp))
      std::cerr << "current pose not in map" << std::endl;

    if (!map.inGrid(np))
      std::cerr << "next pose not in map" << std::endl;
  }

  return false;
}

// This is Equation 4.5 of Nikhil Bakshi's MS thesis
gu::T gu::GUTS::loss(const gu::Vector2& curr_pos, const gu::Vector2& next_pos,
		 const Grid& map)
{
  Matrix x;
  Vector noise_var;

  Matrix ginv_xtsx = ginv_xtsx_;
  Matrix X_sparseT = X_.transpose();
  Matrix y_inv_one_sparse = inv_y_one_;
  Matrix XTyinv_sparse = X_sparseT * y_inv_one_sparse;
  Matrix Y = Y_;

  if (createDirectionalSensor(curr_pos, next_pos, map, x, noise_var))
  {
    Matrix ginv_xtsx_t = ginv_xtsx + x.transpose() * noise_var.cwiseInverse().asDiagonal() * x;
    Matrix Sig = ginv_xtsx_t.diagonal().cwiseInverse().asDiagonal();

    Matrix hstack(XTyinv_sparse.rows(), XTyinv_sparse.cols()+x.rows());
    hstack << XTyinv_sparse, x.transpose() * noise_var.cwiseInverse().asDiagonal();
    Matrix b = Sig * hstack;

    // Portion of b associated with the observation we already have
    // Columns 0...x.shape[0]-1
    // Example: If x.shape = (3,4) and b.shape = (4,5),
    // we return columns 0,1 of b
    Matrix b1 = b.leftCols(b.cols()-x.rows());

    // Portion of b associated with candidate x
    // Columns (b.shape[1]-x.shape[0])...b.shape[1]
    // Example: If x.shape = (3,4) and b.shape = (4,5), we return columns
    // 2,3,4 of b
    Matrix b2 = b.rightCols(x.rows());


    Matrix xTb = x * beta_tilde_;

    Vector beta_hat_t = b1 * Y.col(0) + b2 * xTb;
    Vector beta_tilde = beta_tilde_;

    gu::T beta_hat_t_max = beta_hat_t.maxCoeff();
    gu::T beta_tilde_max = beta_tilde.maxCoeff();
    gu::T coeff_all_same = 0.0;
    for (int i = 0; i < beta_tilde.rows(); ++i)
    {
      bool est_t = std::round(beta_hat_t(i)) > beta_hat_t_max / 2;
      bool real_t = std::round(beta_tilde(i)) > beta_tilde_max / 2;
      if (est_t != real_t)
      {
	coeff_all_same = 1.0;
      }
    }

    return (beta_hat_t - beta_tilde).norm() + LAMBDA * coeff_all_same;
  }

  // Should never hit this case
  std::cerr << "Something went wrong in the loss function" << std::endl;
  return 0;
}


gu::T gu::GUTS::score(const gu::Vector2& current_pose, const gu::Vector2& next_pose,
		      const Grid& map)
{
  return loss(current_pose, next_pose, map);
}
