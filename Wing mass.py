from module1 import Design

wm = Design()

wm.evaluate_design([0.2325,0.9,1])
def print_dict(d):
    for key,value in d.items():
        print(f"{key}: {value}")
    print('\n')
print_dict(wm.design_parameters)
