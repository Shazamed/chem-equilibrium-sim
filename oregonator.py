import numpy as np
import matplotlib.pyplot as plt
import pickle
from timeit import default_timer as timer


# constants
TIME_STEP = 1e-6
MAX_TIME = 90

def int_divide_round(time, precision):
    product = float(10**precision)
    return int(time * product + 0.5)/product

class System:
    def __init__(self, state_dict):
        self.time_elapsed = state_dict["time"]

        self.add_reactions()
        self.add_species()

        # initialize matrices
        self.conc = np.zeros(len(self.species_list))
        self.del_conc = np.zeros(len(self.species_list))

        self.initialize_concentrations(state_dict)
        self.generate_exponential_matrix()
        self.generate_rate_constants_matrix()


    def add_reactions(self):
        with open("oregonator-reactions.txt", 'r') as f: 
            reactions_read = f.read().splitlines()

        self.reactions_list = []
        self.species_list = []
        self.rate_constant_list = []
        for reaction_line in reactions_read: # iterate over different reaction lines to get k and the reagents and products
            rate_constant = float(reaction_line.split(" ")[1]) 
            reaction_line = reaction_line.split(" ")[0]
            self.rate_constant_list.append(rate_constant)
            self.reactions_list.append([reaction_side.split(",") for reaction_side in reaction_line.split(">")]) # save reactions into regeants and products

    def add_species(self):
        for reaction in self.reactions_list:
            for side in reaction:
                for species in side:
                    self.species_list.append(species)
        self.species_list = list(dict.fromkeys(sorted(self.species_list))) # get unique species

    def generate_exponential_matrix(self):
        self.powers = np.zeros([len(self.species_list),len(self.reactions_list)])
        self.coeff = np.zeros([len(self.species_list),len(self.reactions_list)])
        for reaction_idx, reaction in enumerate(self.reactions_list):
            for reagent in reaction[0]:
                species_idx = self.species_list.index(reagent)
                self.powers[species_idx, reaction_idx] += 1
                self.coeff[species_idx, reaction_idx] -= 1
            for product in reaction[1]:
                species_idx = self.species_list.index(product)
                self.coeff[species_idx, reaction_idx] += 1

    def initialize_concentrations(self, state_dict):
        for species_idx, species in enumerate(self.species_list): # initialize concentrations based off of the same index in the np array as the list
            self.conc[species_idx] = state_dict[species]
    
    def generate_rate_constants_matrix(self):
        self.k = np.array(self.rate_constant_list)

    def calc_rate(self):
        self.rate = np.prod(np.array([self.conc]*len(self.reactions_list)).T**self.powers, axis=0)*self.k # matrix before collapsing: rows = concentration, columns = reaction

    def calc_del_conc(self):
        self.del_conc = np.sum(np.array([self.rate]*len(self.species_list))*self.coeff,axis=1)*TIME_STEP

    def calc_conc(self):
        self.conc += self.del_conc

    def step(self):
        self.calc_rate()
        self.calc_del_conc()
        self.calc_conc()

        self.time_elapsed += TIME_STEP
        self.time_elapsed = int_divide_round(self.time_elapsed,6)

        
    def save(self):
        state_dict = {}
        for species_idx, species in enumerate(self.species_list):
            state_dict[species] = self.conc[species_idx]
        state_dict["time"] = self.time_elapsed
        return state_dict
    

def load_initial_state():
    with open('oregonator-concentrations.txt', 'r') as f:
        concentration_lines = f.read().splitlines()
        state_dict = {}
        for line in concentration_lines:
            state_dict[line.split(" ")[0]] = float(eval(line.split(" ")[1]))
    state_dict["time"] = 0
    return state_dict

# arrays for plotting
time_list = []
x_list = []
y_list = []
z_list = []
state_list = []

try:
    state_list = pickle.load(open("oregonator-states.pickle", "rb"))
    last_state = state_list[-1]
    state = System(last_state)
    for state_dict in state_list:    
        time_list.append(state_dict["time"])
        x_list.append(state_dict["X"])
        y_list.append(state_dict["Y"])
        z_list.append(state_dict["Z"])
    print("Data loaded")
except (OSError) as e:
    print("file does not exist")
    state = System(load_initial_state())


while state.time_elapsed < MAX_TIME:
    state.step()
    if state.time_elapsed == int_divide_round(state.time_elapsed,2):
        new_state = state.save()
        state_list.append(new_state)
        time_list.append(new_state["time"])
        x_list.append(new_state["X"])
        y_list.append(new_state["Y"])
        z_list.append(new_state["Z"])
        print(f"time elapsed={state.time_elapsed}s")

print("Iterations done")

with open("oregonator-states.pickle", 'wb') as f:
    pickle.dump(state_list, f)

fig = plt.figure()
ax = fig.add_subplot()
ax.set_yscale('log')
ax.set_title("Concentrations of relevant species in the Oregonator over time")
ax.plot(time_list, x_list,label="X")
ax.plot(time_list, y_list,label="Y")
ax.plot(time_list, z_list,label="Z")
ax.set_xlim(0,MAX_TIME)
ax.set_ylabel("Concentration / M")
ax.set_xlabel("t / s")
ax.legend()
plt.show()