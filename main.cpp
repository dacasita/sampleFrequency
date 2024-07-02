#define _USE_MATH_DEFINES
#include <cmath>
#include "../matplotlibcpp.h"
#include "myfunctions.h"
#include <vector>
#include <complex>
#include <unordered_map>
#include <iostream>

namespace plt = matplotlibcpp;
#include <opencv2/opencv.hpp>

using namespace cv;

int main() 
{
  int k = 11;
  int N = 1 + 2 * k;

  // step 1
  std::unordered_map<double, std::complex<double>> Hd;
  std::vector<double> frequencies;
  std::vector<double> magnitudes;
  for(double m = 0; m < N; m++)
  {
    double f = m / N;
    Hd[f] = std::complex<double>(0, 1) * 2.0 * M_PI * f;

    frequencies.push_back(f);
    magnitudes.push_back(std::abs(Hd[f])); 
  }

  plt::plot(frequencies, magnitudes);
  plt::title("Frequency vs Magnitude");
  plt::xlabel("Frequency");
  plt::ylabel("Magnitude");
  plt::show();

  // step 2
  std::vector<std::complex<double>> r1_n(N, std::complex<double>(0,0));
  for(int n = 0; n < N; n++)
  {
    for(double m = 0; m < N; m++)
    {
      double angle = 2 * M_PI * m * n / N;
      r1_n[n] += Hd[m / N] * std::exp(std::complex<double>(0, angle));
    }
    r1_n[n] /= N;
  }

  // step 2 to step 4
  std::vector<std::complex<double>> h = makeR1intoH(N, r1_n);

    std::vector<double> abs_h;
    for(auto h_i : h)
    {
      abs_h.push_back(std::abs(h_i));
    }

    plt::plot(abs_h);
    plt::title("abs of h after sample in step 4");
    plt::show();

    std::unordered_map<double, std::complex<double>> R_F;
    std::vector<double> F_sample;
    std::vector<double> R_magnitudes;
    double df = 0.01;
    for(double F = 0; F <= 1; F += df)
    {
      R_F[F] = 0;
      for(int n = -N / 2; n <= N / 2; n++)
      {
        double angle = -2.0 * M_PI * F * n;
        // step 4, h[n]=r[n-8]
        R_F[F] += h[n + N/2] * std::exp(std::complex<double>(0, angle));
        //std::cout << F << " : " << std::abs(h[n + N/2] ) << std::endl;
      }
      F_sample.push_back(F);
      R_magnitudes.push_back(std::abs(R_F[F])); 
      //std::cout << std::abs(R_F[F]) << std::endl;
    }

    std::vector<double> abs_R1;
    for(double F = 0; F <= 1; F += df)
    {
      R_F[F] = 0;
      for(int n = 0; n < N; n++)
      {
        double angle = -2.0 * M_PI * F * n;
        R_F[F] += r1_n[n] * std::exp(std::complex<double>(0, angle));
      }
      abs_R1.push_back(std::abs(R_F[F])); 
    }
    plt::plot(F_sample, R_magnitudes, {{"label", "R = fft(r[n])"}});
    plt::plot(F_sample, abs_R1, {{"label", "R1 = fft(r1[n])"}});
    plt::title("R(F)");
    plt::xlabel("F");
    plt::ylabel("R(F)");
    plt::show();

  return 0;
}

//0~N-1
