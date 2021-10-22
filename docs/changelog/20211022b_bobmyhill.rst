* A new BoundaryLayerPerturbation class has been implemented that allows the user
  to create a thermal perturbation to a planetary layer according to the model proposed
  by :cite:`Richter1981`. This perturbation is exponential, taking the form:
  T = a*exp((r - r1)/(r0 - r1)*c) + b*exp((r - r0)/(r1 - r0)*c). The user defines
  a, b and c using three related parameters: rayleigh_number (Ra = c^4),
  temperature_change (the total difference in temperature perturbation
  across the layer, (a + b)*exp(c)), and boundary_layer_ratio (a/b).
  Once initialised, a BoundaryLayerPerturbation object can be modified
  using the function set_model_thermal_gradients,
  which sets the thermal gradient (dT/dr) at the bottom and top of the layer).

  *Bob Myhill, 2021/10/22*
