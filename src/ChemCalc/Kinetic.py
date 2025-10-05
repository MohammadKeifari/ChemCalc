from _general import Enviroment
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import random
from itertools import count
import threading
import time

class KineticalCalculator:
    """
    Simulates chemical reaction kinetics within an Enviroment instance.

    This class numerically integrates concentration changes of compounds over time
    using the rate constants and stoichiometric relationships defined in an `Enviroment` object.

    Attributes:
        accuracy (float): Time step for numerical integration (default: 1e-3).
        fitted (bool): Indicates whether the calculator has been linked to an `Enviroment` instance.
        enviroment (Enviroment): The fitted reaction environment (after calling `fit`).
        rate_constants (list[list[float]]): List of forward and backward rate constants for each reaction.
        reactions_by_index (list[list[list[int]]]): Index mapping of reactants and products per reaction.
        stoichiometric_coefficient_by_reaction (list[list[list[float]]]): Stoichiometric coefficients per reaction.
        rate_dependency_by_reaction (list[list[list[float]]]): Reaction rate dependencies on each compound.
        number_of_reactions (int): Number of reactions in the environment.
        concentrations (list[float]): Current concentration values for each compound in the environment.
    """
    def __init__(self , accuracy = 1e-3):
        """
        Initialize the kinetic calculator with a specified numerical accuracy.

        Args:
            accuracy (float, optional): Time step (Δt) for concentration updates.
                Smaller values yield higher accuracy but slower computation.
                Default is 1e-3.
        """
        self.accuracy = accuracy
        self.fitted = False
    def fit(self , enviroment):
        """
        Link the calculator to an existing `Enviroment` instance.

        Args:
            enviroment (Enviroment): The reaction environment containing all reactions and compounds.

        Raises:
            ValueError: If `enviroment` is not an instance of `Enviroment`.
        """
        if type(enviroment) == Enviroment :
            self.enviroment = enviroment
        else:
            raise ValueError("The input should be an instance of Enviroment class")
        self.rate_constants = enviroment.rate_constants
        self.reactions_by_index = enviroment.reaction_by_index 
        self.stoichiometric_coefficient_by_reaction = enviroment.stoichiometric_coefficient_by_reaction
        self.rate_dependency_by_reaction = enviroment.rate_dependency_by_reaction
        self.number_of_reactions = len(enviroment)
        self.concentrations = []
        for compound in enviroment.compounds_concentration :
            self.concentrations.append(compound["concentration"])
        self.fitted = True
    def calculate(self  , time , checkpoint_time = [] , plot = False  ,concentration_below_zero = "SetToZero"):
        """
        Numerically integrate reaction kinetics over a given time period.

        Supports plotting concentration vs. time and extracting concentration values at checkpoints.

        Args:
            time (float): Total simulation time.
            checkpoint_time (list[float], optional): Specific times to record concentration snapshots.
            plot (bool, optional): Whether to visualize the reaction progress in real-time.
            concentration_below_zero (str, optional): Behavior when a concentration would fall below zero.
                - "SetToZero": Clamp to zero (default)
                - "DoNotChange": Keep previous value
                - "GoNegetive": Allow negative concentrations

        Returns:
            list[list]: List of checkpoint data as `[time, concentrations]` pairs.

        Raises:
            NameError: If `fit` was not called before calculation.
            ValueError: For invalid `concentration_below_zero` parameter.
        """

        if not self.fitted :
            raise NameError("You should fir the model to an enviromt object before calculation")
        
        plt.figure()
        colors = []
        if plot:
            for i in self.enviroment.compounds :
                colors.append(((random.randint(0, 95)/100) , (random.randint(0, 95)/100) , (random.randint(0, 95)/100)))
                plt.xlabel("time")
                plt.ylabel("concentration")
                
        
        compounds_lenght = len(self.concentrations)
        checkpoints = [["time" , self.enviroment.compounds_unicode_formula]]
        def calculate_concentration_change():
            """Compute instantaneous change in concentrations for all compounds."""
            concentration_change = [0] * compounds_lenght
            for rxn_index in range(self.number_of_reactions) :
                rf = self.accuracy * self.rate_constants[rxn_index][0]
                counter = 0
                for compound in self.reactions_by_index[rxn_index][0]:
                    if not(self.concentrations[compound] == 0 and self.rate_dependency_by_reaction[rxn_index][0][counter] < 0):
                        rf = rf * (self.concentrations[compound] ** self.rate_dependency_by_reaction[rxn_index][0][counter])
                    else: 
                        rf = 0
                    counter += 1
                counter = 0
                rb = self.accuracy * self.rate_constants[rxn_index][1]
                for compound in self.reactions_by_index[rxn_index][1]:
                    if not(self.concentrations[compound] == 0 and self.rate_dependency_by_reaction[rxn_index][1][counter] < 0):
                        rb = rb * (self.concentrations[compound] ** self.rate_dependency_by_reaction[rxn_index][1][counter])
                    else:
                        rb = 0
                    counter += 1
                counter = 0
                for reactant in self.reactions_by_index[rxn_index][0]:
                    concentration_change[reactant] += (rb - rf) * self.stoichiometric_coefficient_by_reaction[rxn_index][0][counter] 
                    counter += 1
                counter = 0
                for product in self.reactions_by_index[rxn_index][1]:
                    concentration_change[product] += (rf - rb) * self.stoichiometric_coefficient_by_reaction[rxn_index][1][counter]
                    counter += 1
            return concentration_change
        t = 0
        for i in range(int(time/self.accuracy+1)) :
            delta_concentration = calculate_concentration_change()
            for j in range(len(delta_concentration)):
                concentration = self.concentrations[j] + delta_concentration[j]
                if concentration > 0 :
                    self.concentrations[j] = concentration
                elif concentration_below_zero == "DoNotChange":
                    pass
                elif concentration_below_zero == "GoNegetive":
                    self.concentrations[j] = concentration
                elif concentration_below_zero == "SetToZero":
                    self.concentrations[j] = 0
                else:
                    raise ValueError("concentration_below_zero must be one of: "
                                     "SetToZero, DoNotChange, GoNegetive")
            if plot :
                for k in range(len(self.concentrations)):
                    plt.plot([t , t-self.accuracy],[self.concentrations[k] , self.concentrations[k]-delta_concentration[k]] , color = colors[k])
            for checkpoint_t in checkpoint_time:
                
                if t <= checkpoint_t < t + self.accuracy:
                    checkpoints.append([checkpoint_t , self.concentrations.copy()])
                    
            t += self.accuracy
        if plot :
            for k in range(len(self.concentrations)):
                plt.plot([0 , 0],[0 , 0] , color = colors[k], label = self.enviroment.compounds[k].unicode_formula)
            plt.legend()
            plt.show(block = False)
            
            print("Type 'exit' to close the plot:")
            exited = False
            while not exited:
                cmd = input().strip().lower()
                if cmd == "exit":
                    plt.close()
                    break
                else:
                    print("Invalid input.")
        checkpoints.append([time , self.concentrations])
        
        
        return checkpoints
    
    def calculate_responsively(self, checkpoint_time = [] ,  animation_update_interval = 0.1 , plot = False , concentration_below_zero = "SetToZero"):
        """
        Continuously update and visualize reaction progress in real-time.

        This method is similar to `calculate()`, but uses an interactive animation
        (`matplotlib.animation.FuncAnimation`) for responsive plotting.

        Args:
            checkpoint_time (list[float], optional): Times to record intermediate concentrations.
            animation_update_interval (float, optional): Update interval for live animation (seconds).
            plot (bool, optional): Whether to visualize concentrations dynamically.
            concentration_below_zero (str, optional): How to handle negative concentrations.
                Same as in `calculate()`.

        Returns:
            list[list]: List of checkpoint data `[time, concentrations]`.
        """
        if not self.fitted :
            raise NameError("You must fit the model to an Enviroment before calculation.")
        
        plt.figure()
            
        colors = []
        if plot:
            for i in self.enviroment.compounds :
                colors.append(((random.randint(0, 95)/100) , (random.randint(0, 95)/100) , (random.randint(0, 95)/100)))
                plt.xlabel("time")
                plt.ylabel("concentration")
                
        
        compounds_lenght = len(self.concentrations)
        checkpoints = [["time" , self.enviroment.compounds_unicode_formula]]
        def calculate_concentration_change():
            """Compute instantaneous concentration change per reaction step."""
            concentration_change = [0] * compounds_lenght
            for rxn_index in range(self.number_of_reactions) :
                rf = self.accuracy * self.rate_constants[rxn_index][0]
                counter = 0
                for compound in self.reactions_by_index[rxn_index][0]:
                    if not(self.concentrations[compound] == 0 and self.rate_dependency_by_reaction[rxn_index][0][counter] < 0):
                        rf = rf * (self.concentrations[compound] ** self.rate_dependency_by_reaction[rxn_index][0][counter])
                    else: 
                        rf = 0
                    counter += 1
                counter = 0
                rb = self.accuracy * self.rate_constants[rxn_index][1]
                for compound in self.reactions_by_index[rxn_index][1]:
                    if not(self.concentrations[compound] == 0 and self.rate_dependency_by_reaction[rxn_index][1][counter] < 0):
                        rb = rb * (self.concentrations[compound] ** self.rate_dependency_by_reaction[rxn_index][1][counter])
                    else:
                        rb = 0
                    counter += 1
                counter = 0
                for reactant in self.reactions_by_index[rxn_index][0]:
                    concentration_change[reactant] += (rb - rf) * self.stoichiometric_coefficient_by_reaction[rxn_index][0][counter] 
                    counter += 1
                counter = 0
                for product in self.reactions_by_index[rxn_index][1]:
                    concentration_change[product] += (rf - rb) * self.stoichiometric_coefficient_by_reaction[rxn_index][1][counter]
                    counter += 1
            return concentration_change
        time = count()
        def animate(i):
            """Animation update loop for real-time kinetics visualization."""
            t = self.accuracy * next(time)
            delta_concentration = calculate_concentration_change()
            for j in range(len(delta_concentration)):
                concentration = self.concentrations[j] + delta_concentration[j]
                if concentration > 0 :
                    self.concentrations[j] = concentration
                elif concentration_below_zero == "DoNotChange":
                    pass
                elif concentration_below_zero == "GoNegetive":
                    self.concentrations[j] = concentration
                elif concentration_below_zero == "SetToZero":
                    self.concentrations[j] = 0
                else:
                    raise ValueError("concentration_below_zero should have one of the following values: SetToZero\\DoNotChange\\GoNegetive")
            if plot :
                for k in range(len(self.concentrations)):
                    plt.plot([t , t-self.accuracy],[self.concentrations[k] , self.concentrations[k]-delta_concentration[k]] , color = colors[k])
            for checkpoint_t in checkpoint_time:
                
                if t <= checkpoint_t < t + self.accuracy:
                    checkpoints.append([checkpoint_t , self.concentrations.copy()])


        ani = FuncAnimation(plt.gcf() , animate , interval = animation_update_interval , cache_frame_data=False)
        if plot :
            for k in range(len(self.concentrations)):
                plt.plot([0 , 0],[0 , 0] , color = colors[k], label = self.enviroment.compounds[k].unicode_formula)

            plt.legend()
            plt.show(block = False)
            print("Type 'exit' to close / 'stop' to pause / 'resume' to continue:")
            while True:
                cmd = input().strip().lower()
                if cmd == "exit":
                    ani.pause()
                    plt.close()
                    break
                elif cmd == "stop":
                    ani.pause()
                elif cmd == "resume":
                    ani.resume()
                else:
                    print("invalid_input")
        checkpoints.append(["end" , self.concentrations])
        return checkpoints