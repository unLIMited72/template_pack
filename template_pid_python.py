class TemplatePID:
  def __init__(
      self, kp: float = 1.0, ki: float = 0.0, kd: float = 0.0,
      u_max: float = float('inf'), u_min: float = -float('inf'),
      i_limit: float = float('inf'), i_aw_enable: bool = True,
      i_err_overshoot_std: float = 0.0, i_ratio_with_p: float = 0.5,
      i_leak_ratio: float = 0.0, d_tau: float = 0.0, d_limit: float = float('inf'), 
      d_deadband: float = 0.0, d_on_measure_sw: bool = True, 
      d_sat_mode: str = "cut", d_sat_gain: float = 0.3):
# Gain Parameters
    self.kp = kp
    self.ki = ki
    self.kd = kd
# Saturation Control Input
    self.u_max = u_max
    self.u_min = u_min
# For I Controller
    self.i_integral = 0.0
    self.i_limit = i_limit
    self.i_prev_err = None
    self.i_aw_enable = i_aw_enable
    self.i_err_overshoot_std = i_err_overshoot_std
    self.i_ratio_with_p = i_ratio_with_p
    self.i_leak_ratio = i_leak_ratio
# For D Controller
    self.d_tau = d_tau
    self.d_limit = d_limit
    self.d_deadband = d_deadband
    self.d_on_measure_sw = d_on_measure_sw
    self.d_prev_signal = None
    self.d_prev_filter = 0.0
    self.d_sat_mode = d_sat_mode # none | cut | blend
    self.d_sat_gain = d_sat_gain
  
  def _update_gains(self, kp: float, ki: float, kd: float):
    self.kp = kp
    self.ki = ki
    self.kd = kd

  def _reset(self):
    self.i_integral = 0.0
    self.i_prev_err = None
    self.d_prev_signal = None
    self.d_prev_filter = 0.0

  def _clamp(self, x: float, x_min: float, x_max: float) -> float:
    if x > x_max : return x_max
    if x < x_min : return x_min
    return x
  
  def _discrete_LPF(self, u: float, Ts: float) -> float:
    if self.d_tau > 0.0:
      alpha = Ts / (self.d_tau + Ts)
      filter = ((1 - alpha) * self.d_prev_filter) + (alpha * u)
    else:
      filter = u
    return filter
  
  def _integral(self, err: float, dt: float):
    future_integral = self.i_integral + (err * dt)
    
    skip_flag = False
    future_p = self.kp * err
    future_i = self.ki * future_integral
    future_d = self.kd * self.d_prev_filter if self.kd != 0.0 else 0.0
    future_u = future_p + future_i + future_d
    overshoot_flag = (
      (self.i_prev_err is not None) and (err * self.i_prev_err < 0.0) and
      (abs(err) > self.i_err_overshoot_std))
    opposes_flag = (future_p * future_i) < 0.0
    ratio_flag = abs(future_i) > (self.i_ratio_with_p * abs(future_p))
    i_aw_hinder_flag = overshoot_flag and opposes_flag and ratio_flag
    if self.i_aw_enable:
      if ((future_u > self.u_max) and (err > 0)) or ((future_u < self.u_min) and (err < 0)):
        skip_flag = True
      if i_aw_hinder_flag:
        if self.i_leak_ratio > 0.0:
          future_integral = (1.0 - self.i_leak_ratio) * self.i_integral
        else:
          skip_flag = True
    if not skip_flag: self.i_integral = future_integral

    self.i_integral = self._clamp(self.i_integral, -self.i_limit, self.i_limit)

  def _derivative(self, ref: float, feedback: float, dt: float, future_p: float, future_i: float) -> float:
    if self.kd == 0.0: return 0.0

    signal = feedback if self.d_on_measure_sw else (ref - feedback)
    if self.d_prev_signal is None:
      self.d_prev_signal = signal
      self.d_prev_filter = 0.0
      return 0.0
    
    delta = signal - self.d_prev_signal
    raw_deriv = 0.0 if (abs(delta) < self.d_deadband) else (delta / dt)
    raw_deriv = self._clamp(raw_deriv, -self.d_limit, self.d_limit)

    filter_deriv = self._discrete_LPF(raw_deriv, dt)
    filter_deriv = self._clamp(filter_deriv, -self.d_limit, self.d_limit)

    if self.d_sat_mode != "none":
      future_d = self.kd * filter_deriv
      future_u = future_p + future_i + future_d
      saturation_flag = (
        (future_u > self.u_max) and (future_d > 0.0) or 
        (future_u < self.u_min) and (future_d < 0.0))
      if saturation_flag:
        if self.d_sat_mode == "cut":
          filter_deriv = 0.0
        elif self.d_sat_mode == "blend":
          filter_deriv = self.d_sat_gain * filter_deriv
        else:
          pass

    self.d_prev_signal = signal
    self.d_prev_filter = filter_deriv
    return filter_deriv
  
  def loop(self, ref: float, feedback: float, dt: float) -> float:
    err = ref - feedback

    self._integral(err = err, dt = dt)

    control_input_p = self.kp * err
    control_input_i = self.ki * self.i_integral

    d_term = self._derivative(
      ref = ref, feedback = feedback, dt = dt, 
      future_p = control_input_p,  future_i = control_input_i)
    control_input_d = self.kd * d_term

    u_unsat = control_input_p + control_input_i + control_input_d
    u_t = self._clamp(u_unsat, self.u_min, self.u_max)

    self.i_prev_err = err

    return u_t