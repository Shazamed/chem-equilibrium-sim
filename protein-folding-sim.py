import numpy as np
import matplotlib.pyplot as plt
import pickle

# constants
TIME_STEP = 1e-6
MAX_TIME = 5

urea_conc_range = np.arange(3,6.5,0.2) # specify range of urea concentrations


def int_divide_round(time, precision):
    product = float(10**precision)
    return int(time * product + 0.5)/product

class System:
    def __init__(self, state_dict):
        self.time_elapsed = state_dict["time"]
        self.urea_conc = state_dict["urea"]
        self.add_reactions()
        self.add_species()

        # initialize matrices
        self.conc = np.zeros(len(self.species_list))
        self.del_conc = np.zeros(len(self.species_list))

        self.initialize_concentrations(state_dict)
        self.generate_exponential_matrix()
        self.generate_rate_constants_matrix()


    def add_reactions(self):
        with open("protein-reactions.txt", 'r') as f: 
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

    def initialize_urea_scaling(self):
        urea_scale = np.array([
            np.e**(-1.68*self.urea_conc),
            np.e**(0.95*self.urea_conc),
            np.e**(-1.72*self.urea_conc),
            np.e**(1.20*self.urea_conc),
        ])
        return urea_scale
    
    def generate_rate_constants_matrix(self):
        self.k = np.array(self.rate_constant_list) * np.array(self.initialize_urea_scaling())

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
        return {"time": self.time_elapsed, 
                "D": self.conc[0], 
                "I": self.conc[1], 
                "N": self.conc[2], 
                "urea": self.urea_conc
                }
    

def load_initial_state(urea_conc):
    with open('protein-concentrations.txt', 'r') as f:
        concentration_lines = f.read().splitlines()
        state_dict = {}
        for line in concentration_lines:
            state_dict[line.split(" ")[0]] = float(eval(line.split(" ")[1]))
    state_dict["time"] = 0
    state_dict["urea"] = urea_conc
    state_list = [state_dict]
    return state_list


urea_list = []
d_list = []
i_list = []
n_list = []
eq_list = []


try:
    eq_list = pickle.load(open("protein-eq.pickle", "rb"))
    last_state = eq_list[-1]
    state = System(last_state)
    for state_dict in eq_list:    
        urea_list.append(state_dict["urea"])
        d_list.append(state_dict["D"])
        i_list.append(state_dict["I"])
        n_list.append(state_dict["N"])
    print("Data loaded")
except (OSError) as e:
    print("file does not exist")



for urea in urea_conc_range:
    state = System(load_initial_state(urea)[0])
    while state.time_elapsed < MAX_TIME:
        state.step()
        # if int_divide_round(state.time_elapsed,6) == int_divide_round(state.time_elapsed,2):
        #     new_state = state.save()
        #     state_list.append(new_state)
        #     time_list.append(new_state["time"])
        #     d_list.append(new_state["D"])
        #     i_list.append(new_state["I"])
        #     n_list.append(new_state["N"])
        #     print(f"time elapsed={state.time_elapsed}"
    new_state = state.save()
    eq_list.append(new_state)
    urea_list.append(urea)
    d_list.append(new_state["D"])
    i_list.append(new_state["I"])
    n_list.append(new_state["N"])
    print(urea)
        

print("Iterations done")

with open("protein-eq.pickle", 'wb') as f:
    pickle.dump(eq_list, f)

fig = plt.figure()
ax = fig.add_subplot()
# ax.set_yscale('log')
ax.set_title("Equilibrium fraction of protein states against urea concentration")
ax.scatter(urea_list, d_list,label="D", marker='x')
ax.scatter(urea_list, i_list,label="I", marker='+')
ax.scatter(urea_list, n_list,label="N", marker='.')
ax.legend()

plt.show()

