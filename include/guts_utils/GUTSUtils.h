#pragma once

/*
Copyright 2023 Wennie Tabib

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
“AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#include <random>
#include <cstdlib>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Dense>

#include <guts_utils/Grid.h>

// References:
// [1] "Multi Agent Active Search using Realistic Depth-Aware Noise
// Model" DOI: 10.1109/ICRA48506.2021.9561598

namespace guts_utils
{

  using T = double;
  using Matrix = Eigen::Matrix<T, -1, -1>;
  using Vector = Eigen::Matrix<T, -1, 1>;
  using Vector2 = Eigen::Matrix<T, 2, 1>;
  using Grid = guts_utils::Grid;
  using Point = guts_utils::Grid::Point;
  using Cell = guts_utils::Cell;

  class StandardNormalGenerator
  {
  public:
    StandardNormalGenerator() {}
    ~StandardNormalGenerator() {}
    virtual Vector sample(unsigned int nsamples) = 0;
  };

  class StandardNormal : public StandardNormalGenerator
  {
  public:
    StandardNormal() {}
    ~StandardNormal() {}

    Vector sample(unsigned int nsamples)
    {
      // Change the seed every time the sample function is called
      // Ensures random numbers are different every time
      std::random_device rd;
      std::default_random_engine generator;
      generator_.seed( rd() );

      Vector samples = Vector::Zero(nsamples,1);
      for (unsigned int i = 0; i < nsamples; ++i)
      {
	samples(i) = distribution_(generator_);
      }
      return samples;
    }

  private:
    std::default_random_engine generator_;
    std::normal_distribution<double> distribution_{0.0,1.0};
  };

  struct GUTS
  {

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Ptr = std::shared_ptr<GUTS>;
    using ConstPtr = std::shared_ptr<const GUTS>;

    GUTS();
    ~GUTS();

    // M is the number of gridded environment cells
    // Input:
    // mean: (M x 1)
    // cov:  (M x M)
    //
    // Output: (M x 1)
    Matrix multivariateNormal(const Vector& mean, const Matrix& cov);

    // This function calculates Equation (2)--(4) of of [1]
    //
    // Note: in this implementation, Q is 1.
    //
    // B is the gridded environment with dimension M1 x M2
    // beta is the flattened version of B, meaning it is M x 1, where M = M1*M2
    // X is the sensing matrix for time 1...t-1 available to the agent at time step t (comes from D)
    //   it has dimension (num_observations x M). While the paper discusses a Q variables, which
    //   represents the FOV grid points available, the implementation here assumes Q = 1.
    // Y combines all the yt from [1] and includes the noise vector.
    //   yt combines the y variable in [1], which is Q x 1 with the noise vector, n,
    //   which is Q x 1 (i.e., Y is Q x 2).
    //
    // posterior(beta | data) = N(mu, V)
    // V = Sig_beta
    // mu = beta_hat_tmp
    void posterior(const Matrix& X, const Matrix& Y);

    // Uses C = (A+A.T) * 0.5 according to https://stackoverflow.com/a/41692936
    void makeCovariancePSD(Matrix& C);

    // Creates observations from current pose to goal pose
    //
    // Assumption: the inclusion zone is convex.
    //
    // createDirectionalSensor creates observations raytraces a line
    // from curr_pos to next_pos
    bool createDirectionalSensor(const Vector2& curr_pos, const Vector2& next_pos,
				 const Grid& map, Matrix& x, Vector& noise_var);

    // This is Equation 4.5 of Nikhil Bakshi's MS thesis
    T loss(const Vector2& curr_pos, const Vector2& next_pos, const Grid& map);


    // It is assumed that at the start of each planning round, the user
    // will call the posterior function.
    //
    // Then the user will call the score function on one or more candidate
    // grid locations.
    T score(const Vector2& current_pose, const Vector2& next_pose,
	    const Grid& map);

    ///// debugging
    Vector beta_hat_;
    Vector gamma_;
    Vector beta_tilde_;
    T min_eig_;
    Matrix Sig_beta_;
    Matrix ginv_xtsx_;
    Matrix Gamma_;
    Matrix inv_y_one_;
    Matrix X_;
    Matrix Y_;
    ///// End debugging

    unsigned int EMitr_;
    unsigned int b;

    T a;
    T eps_;

  };
}
