from ChemCompute import *
from ChemCompute.Kinetic import KineticalCalculator
import numpy as np
import time


# rxn1 = Reaction.from_string_simple_syntax("A > B" , [1.1 , 1] , enthalpy = 100, entropy = 100, K = 100 , T=500)
# rxn2 = Reaction.from_string_simple_syntax("B > C" , [1 , 5] , enthalpy = 100, entropy = 100, K = 100)
# rxn3 = Reaction.from_string_simple_syntax("C > 2D " , [1 , 0])


rxn1 = Reaction.from_string_simple_syntax("g + a > 2a + g" , kf = 0.1 , kb = 0)
rxn2 = Reaction.from_string_simple_syntax("a + b > 2b" , kf = 0.1 , kb = 0)
rxn3 = Reaction.from_string_simple_syntax("b > d" , kf = 0.1 , kb = 0)

env = Enviroment(rxn1 , rxn2 , rxn3)
env.concentrations = [3,2,1,0]
kc = KineticalCalculator(1e-3)
kc.fit(env)



t1 = time.time()
concentrations = kc.calculate(time = 100 , plot = "save")
print(concentrations)
t2 = time.time()
print(t2-t1)