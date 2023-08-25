import scipy.optimize as opt
import math
import numpy as np
from module1 import Design

solar_plane = Design()

epsilon = 1e-3
bounds=[(epsilon,0.8),(epsilon,4),(0.6,1.2)]

res = opt.shgo(
    solar_plane.objective,
    bounds,
    n=250000,
    constraints=[
        {'type':'ineq','fun':solar_plane.con1},
        {'type':'ineq','fun':solar_plane.con2},
        {'type':'ineq','fun':solar_plane.con3},
        {'type':'ineq','fun':solar_plane.con4}
        ],
    sampling_method = 'sobol'
    )


def print_dict(d):
    for key,value in d.items():
        print(f"{key}: {value}")
    print('\n')


#print(res)
#print_dict(solar_plane.design_parameters)
soln_list = res.xl[:20] if len(res.xl >= 20) else res.xl
for solution in soln_list:
    solar_plane.evaluate_design(solution)
    print_dict(solar_plane.design_parameters)

    print("Constraint values are the safety margin, negative values mean constraint is violated")
    print_dict(solar_plane.constraints)


