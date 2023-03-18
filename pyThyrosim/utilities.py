import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate._ivp.common import norm
from scipy._lib._util import check_random_state
from scipy.sparse import issparse
from scipy.sparse.linalg import norm as sparse_norm
from ThyrosimIM import initialize, find_patient_ic, thyrosim

"""
    simulate(h, w, sex, ...)

Simulate a person of known height, weight, and gender for 30 days (default).

If `warmup = True`, will first run the model for 30 days, assuming healthy
thyroid function, to get approximate initial condition. 

Note: The last 5 parameters are optional scaling parameters, but in our final model,
only `scale_plasma_ode` and `scale_clearance_by_gender` are set to `True`. Thus only 
these 2 are True by default. 
"""
def simulate(
    h, # units meters
    w, # units kg
    sex:bool, # True = male, False = female
    days = 30, 
    dial=[1.0, 0.88, 1.0, 0.88], 
    T4dose=0.0, # mcgs
    T3dose=0.0, # mcgs
    dosing_interval=24.0, #hours
    warmup:bool = True,
    fitting_index = [],
    parameters = [],
    fixed_parameters= tuple([]),
    scale_plasma_ode=True,
    scale_slow_ode=False,
    scale_fast_ode=False,
    scale_allometric_exponent = False,
    scale_clearance_by_gender = True,
    ):
    def add_dose(integrator):
        integrator.u[9] += integrator.p[54]
        integrator.u[11] += integrator.p[55]
        return integrator
    
    cbk = PeriodicCallback(add_dose, dosing_interval)

    # initialize thyrosim parameters
    ic, p = initialize(dial, True, h, w, sex, 
        fitting_index=fitting_index, p_being_optimized=parameters,
        fixed_parameters=fixed_parameters,
        scale_plasma_ode=scale_plasma_ode, scale_slow_ode=scale_slow_ode,
        scale_fast_ode=scale_fast_ode,
        scale_allometric_exponent=scale_allometric_exponent,
        scale_clearance_by_gender=scale_clearance_by_gender)
    p[fitting_index] = parameters

    # run simulation for 30 days to get approximate steady state conditions
    # this assumes healthy patient without dose
    ic = find_patient_ic(ic, p, 30) if warmup else ic

    # setup daily dosing and fitting parameters 
    p[54] = T4dose / 777.0 # daily dose
    p[55] = T3dose / 651.0 # daily dose
    p[56:60] = dial #set dial

    # solve and return ode solution
    tspan = (0.0, 24*days)
    prob = solve_ivp(thyrosim,tspan,ic,args=(p))
    sol = [sol.x, prob.y]

    return sol




class PeriodicCallback: # This is super hacky and probably won't work, but let's thyrosim run properly in absenece of T3/T4 dosing
    def __init__(self, callback, period, args=(), kwargs={}, first_call=None):
        self.callback = callback
        self.period = period
        self.args = args
        self.kwargs = kwargs
        self.first_call = first_call
        self.t_last = None

    def __call__(self, t, y):
        if self.t_last is None:
            if self.first_call is not None:
                self.callback(self.first_call, y, *self.args, **self.kwargs)
            self.t_last = t
            return None
        if (t - self.t_last) >= self.period:
            n_periods = int((t - self.t_last) // self.period)
            t_event = self.t_last + n_periods * self.period
            self.callback(t_event, y, *self.args, **self.kwargs)
            self.t_last = t_event
        return None

