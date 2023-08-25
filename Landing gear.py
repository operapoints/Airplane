import scipy.optimize as opt
import math
import numpy as np
from module1 import landing_gear

lg = landing_gear()

epsilon = 1e-4
#t w L theta c
bounds=[(epsilon/2,1e-2), #cap thickness
        (epsilon,10e-2), #width
        (epsilon,10e-1), #length
        (epsilon,math.pi/2), #angle with ground
        (epsilon,1e-1), #core thickness
        (1,10) #taper ratio
        ]

res = opt.shgo(
    lg.objective,
    bounds,
    n=10000,
    constraints=[
        {'type':'ineq','fun':lg.con1_max_stress},
        {'type':'ineq','fun':lg.con2_max_deflection},
        {'type':'ineq','fun':lg.con3_max_z_accel},
        {'type':'ineq','fun':lg.con4_buckling}
        ],
    sampling_method = 'sobol'
    )


def print_dict(d):
    for key,value in d.items():
        #value = float('%.3g' % value)
        print(f"{key}: {value}")
    print('\n')


print(res)
#print_dict(solar_plane.design_parameters)

soln_list = res.xl[:20] if len(res.xl >= 20) else res.xl

for solution in soln_list:
    lg.evaluate_design(solution)
    print_dict(lg.design_parameters)

    print("Constraint values are the safety margin, negative values mean constraint is violated")
    print_dict(lg.constraints)
