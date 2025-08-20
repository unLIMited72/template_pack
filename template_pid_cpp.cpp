#include <optional>
#include <cmath>
#include <limits>
#include <algorithm>

enum class DSatMode { none, cut, blend }; 

struct GainParam { double kp {1.0}, ki {0.0}, kd {0.0}; }; 

struct SatControlInput { 
  double max {std::numeric_limits<double>::infinity()};
  double min {-std::numeric_limits<double>::infinity()}; 
}; 

struct VarForIControl { 
  double limit {std::numeric_limits<double>::infinity()}; 
  double err_overshoot_std {0.0};
  double ratio_with_p {0.5}; 
  double leak_ratio {0.0}; 
  bool aw_enable {true};
  double integral {0.0};
  std::optional<double> prev_err;
}; 

struct VarForDControl { 
  double tau {0.0};
  double limit {std::numeric_limits<double>::infinity()};
  double deadband {0.0};
  double sat_gain {0.3}; 
  bool on_measure {true}; 
  DSatMode sat_mode {DSatMode::cut};
  double prev_filter {0.0};
  std::optional<double> prev_signal;
}; 

class TemplatePID { 
private: 
  GainParam gains_; 
  SatControlInput u_sat_; 
  VarForIControl i_; 
  VarForDControl d_; 
public: 
  TemplatePID(
    double kp, double ki, double kd, double u_max, double u_min,
    double i_limit, double i_err_overshoot_std, double i_ratio_with_p,
    double i_leak_ratio, bool i_aw_enable, double d_tau, double d_limit,
    double d_deadband, double sat_gain, bool on_measure, 
    DSatMode sat_mode)
    : gains_{kp, ki, kd}, u_sat_{u_max, u_min},
      i_{i_limit, i_err_overshoot_std, i_ratio_with_p, i_leak_ratio, i_aw_enable},
      d_{d_tau, d_limit, d_deadband, sat_gain, on_measure, sat_mode} {}
      
  void update_gains(double, double, double);
  void reset();

  double discrete_LPF(double, double);
  void integral(double, double);
  double derivative(double, double, double, double, double);
  double loop(double, double, double);
};

void TemplatePID::update_gains(double kp, double ki, double kd) {
  this->gains_.kp = kp;
  this->gains_.ki = ki;
  this->gains_.kd = kd;
}

void TemplatePID::reset() {
  this->i_.integral = 0.0;
  this->i_.prev_err.reset();
  this->d_.prev_signal.reset();
  this->d_.prev_filter = 0.0;
}

double TemplatePID::discrete_LPF(double u, double Ts) {
  if (this->d_.tau > 0.0) {
    double alpha = Ts / (this->d_.tau + Ts);
    double filter = ((1 - alpha) * this->d_.prev_filter) + (alpha * u);
    return filter;
  } else {
    return u;
  }
}

void TemplatePID::integral(double err, double dt) {
  double future_integral = this->i_.integral + (err * dt);
  bool skip_flag = false;

  double future_p = this->gains_.kp * err;
  double future_i = this->gains_.ki * future_integral;
  double future_d = (this->gains_.kd != 0.0) ? this->gains_.kd * this->d_.prev_filter : 0.0;
  double future_u = future_p + future_i + future_d;

  bool aw_hinder_flag = 
    (this->i_.prev_err.has_value()) && 
    (((this->i_.prev_err.value()) * err) < 0.0) &&
    (std::abs(err) > this->i_.err_overshoot_std) &&
    ((future_p * future_i) < 0.0) &&
    (std::abs(future_i) > (this->i_.ratio_with_p * std::abs(future_p)));

  if (this->i_.aw_enable) {
    if(((future_u > this->u_sat_.max)) && (err > 0) || ((future_u < this->u_sat_.min) && (err < 0)))
      skip_flag = true;
    if (aw_hinder_flag) {
      if (this->i_.leak_ratio > 0.0) future_integral = (1.0 - this->i_.leak_ratio) * this->i_.integral;
      else skip_flag = true;
    }
  }

  if (!skip_flag) this->i_.integral = future_integral;

  this->i_.integral = std::clamp(this->i_.integral, -this->i_.limit, this->i_.limit);
}

double TemplatePID::derivative(double ref, double feedback, double dt, double f_p, double f_i) {
  if (this->gains_.kd == 0.0) return 0.0;

  double signal = this->d_.on_measure ? feedback : (ref - feedback);
  if (!this->d_.prev_signal.has_value()) {
    this->d_.prev_signal = signal;
    this->d_.prev_filter = 0.0;
    return 0.0;
  }

  double delta = signal - this->d_.prev_signal.value();
  double raw_deriv = (std::abs(delta) < this->d_.deadband) ? 0.0 : (delta / dt);
  raw_deriv = std::clamp(raw_deriv, -this->d_.limit, this->d_.limit);

  double filter_deriv = this->discrete_LPF(raw_deriv, dt);
  filter_deriv = std::clamp(filter_deriv, -this->d_.limit, this->d_.limit);

  if (this->d_.sat_mode != DSatMode::none) {
    double f_d = this->gains_.kd * filter_deriv;
    double f_u = f_p + f_i + f_d;
    if (((f_u > this->u_sat_.max) && (f_d > 0.0)) || ((f_u < this->u_sat_.min) && (f_d < 0.0))) {
      if (this->d_.sat_mode == DSatMode::cut) filter_deriv = 0.0;
      else if (this->d_.sat_mode == DSatMode::blend) filter_deriv = filter_deriv * this->d_.sat_gain;
    }
  }
  this->d_.prev_signal = signal;
  this->d_.prev_filter = filter_deriv;
  return filter_deriv;
}

double TemplatePID::loop(double ref, double feedback, double dt) {
  double err = ref - feedback;
  this->integral(err, dt);

  double control_p = this->gains_.kp * err;
  double control_i = this->gains_.ki * this->i_.integral;
  double d_term = this->derivative(ref, feedback, dt, control_p, control_i);
  double control_d = this->gains_.kd * d_term;

  double u_unsat = control_p + control_i + control_d;
  double u_t = std::clamp(u_unsat, this->u_sat_.min, this->u_sat_.max);

  this->i_.prev_err = err;

  return u_t;
}