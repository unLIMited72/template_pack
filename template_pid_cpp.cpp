#include <iostream>
#include <memory>
#include <string>
#include <functional>
#include <chrono>
#include <cmath>
#include <optional>

enum class DSatMode { none, cut, blend };

struct GainParam { double kp {1.0}, ki {0.0}, kd {0.0}; };

struct SatControlInput { double max {INFINITY}, min {-INFINITY}; };

struct VarForIControl
{
  double integral {0.0}, limit {INFINITY};
  double err_overshoot_std {0.0}, ratio_with_p {0.5}, leak_ratio {0.0};
  std::optional<double> prev_err;
  bool aw_enable {true};
};

struct VarForDControl
{
    double tau {0.0}, limit {INFINITY}, deadband {0.0};
    double prev_signal {0.0}, prev_filter {0.0}, sat_gain {0.3};
    bool on_measure_sw {true};
    DSatMode sat_mode {DSatMode::cut};
};

class TemplatePID
{
private:
  GainParam gains_;
  SatControlInput u_sat_;
  VarForIControl i_;
  VarForDControl d_;
public:
  TemplatePID(
    double kp, double ki, double kd, double u_max, double u_min,
    double i_limit, bool i_aw_enable, double i_err_overshoot_std, 
    double i_ratio_with_p, double i_leak_ratio, double d_tau, 
    double d_limit, double d_deadband, bool d_on_measure_sw, 
    DSatMode d_sat_mode, double d_sat_gain)
  : gains_{kp, ki, kd}, u_sat_{u_max, u_min},
    i_{0.0, i_limit, 0.0, i_err_overshoot_std, i_ratio_with_p, i_leak_ratio, i_aw_enable},
    d_{d_tau, d_limit, d_deadband, 0.0, 0.0, d_sat_gain, d_on_measure_sw, d_sat_mode} {}

  void update_gains(double, double, double);
  void reset();
};


void TemplatePID::update_gains(double kp, double ki, double kd){
  this->gains_.kp = kp;
  this->gains_.ki = ki;
  this->gains_.kd = kd;
}

void TemplatePID::reset(){
  this->i_.integral = 0.0;
  this->
}