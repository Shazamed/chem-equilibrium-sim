import numpy as np
import matplotlib.pyplot as plt
import pickle
from system_class import System

# constants
TIME_STEP = 1e-6
MAX_TIME = 15
RELATIVE_TOLERANCE = 1e-3 # for comparing if equilibrium has been reached
REACTION_PATH = "protein_inputs/protein-reactions.txt"

# specify range of urea concentrations 
urea_conc_range = np.concatenate([np.arange(0,8.5,0.5),np.arange(3,6.5,0.2)])

urea_conc_range = np.unique(urea_conc_range.round(decimals=6)) # remove floating point errors and then remove repeated values

def int_divide_round(time, precision):
    '''Faster rounding to the decimal place indicated by precision'''
    product = float(10**precision)
    return int(time * product + 0.5)/product

def load_initial_state(urea_conc, prev_conc_dict=None):
    '''Read input file for initial concentrations and return the initial state dictionary eg.:
    {
    "X": 0.01,
    "Y": 0.5,
    "time" = 0.0
    }
    '''
    state_dict = {}
    if prev_conc_dict == None:
        with open('protein_inputs/protein-concentrations.txt', 'r') as f:
            concentration_lines = f.read().splitlines()
            for line in concentration_lines:
                state_dict[line.split(" ")[0]] = float(eval(line.split(" ")[1]))
    else:
        state_dict = prev_conc_dict
    state_dict["time"] = 0
    state_dict["urea"] = urea_conc
    return state_dict

class Cell(System): # inherit from System and modify methods to adjust for urea concentrations
    def __init__(self, state_dict, time_step, reaction_path):
        self.urea_conc = state_dict["urea"]
        super().__init__(state_dict, time_step, reaction_path)
    
    def initialize_urea_scaling(self):
        '''Generate numpy array for scaling the rate constant'''
        urea_scaling = np.array([
            np.e**(-1.68*self.urea_conc),
            np.e**(0.95*self.urea_conc),
            np.e**(-1.72*self.urea_conc),
            np.e**(1.20*self.urea_conc),
        ])
        return urea_scaling

    def generate_rate_constants_matrix(self):
        '''Generate the rate constant k'''
        self.k = np.array(self.rate_constant_list) * np.array(self.initialize_urea_scaling()) # multiply by urea exponential scaling
    
    def save(self):
        '''save last state as a dictionary
        {
        "X": 0.01,
        "Y": 0.5,
        "time" = 3.9,
        "urea" = 3.0
        }
        '''
        state_dict = {}
        for species_idx, species in enumerate(self.species_list):
            state_dict[species] = self.conc[species_idx]
        state_dict["time"] = self.time_elapsed
        state_dict['urea'] = self.urea_conc
        return state_dict


# lists for plotting
urea_list = []
d_list = []
i_list = []
n_list = []
eq_list = []

try: # try to load the file if file exist    
    eq_list = pickle.load(open("protein-eq.pickle", "rb"))
    last_state = eq_list[-1]
    state = Cell(last_state, TIME_STEP, REACTION_PATH)  # generate initial conditions from last state saved
    for state_dict in eq_list:    # add saved states to list for plotting
        urea_list.append(state_dict["urea"])
        d_list.append(state_dict["D"])
        i_list.append(state_dict["I"])
        n_list.append(state_dict["N"])
    print("Data loaded")
except (OSError) as e:
    print("file does not exist")

saved_state = None

for urea in urea_conc_range: #  loop over all concentrations in urea range
    if urea in urea_list: # if concentration already exist in saved list, skip calculating
        continue
    print(f"Urea concentration: {urea}M")
    state = Cell(load_initial_state(urea, saved_state), TIME_STEP, REACTION_PATH) # generate state from initial concentrations of the last concentration state or from input for the first run
    while state.time_elapsed < MAX_TIME: # set max time if equilibration takes too long 
        state.step()
        if state.time_elapsed == int_divide_round(state.time_elapsed,1):
            print(f"time elapsed={state.time_elapsed}s")
            if state.time_elapsed >= 0.2:
                if np.allclose(prev_conc, state.conc, atol=0, rtol=RELATIVE_TOLERANCE): # check if concentrations equilibrated within tolerance
                    break
            prev_conc = state.conc.copy() # get previous concentration for comparison


    # get the last state as a dictionary and save to lists for plotting
    saved_state = state.save()
    eq_list.append(saved_state.copy())
    urea_list.append(urea)
    d_list.append(saved_state["D"])
    i_list.append(saved_state["I"])
    n_list.append(saved_state["N"])
    
print("Iterations done")

with open("protein-eq.pickle", 'wb') as f:
    pickle.dump(eq_list, f)

total_conc = np.array(d_list) + np.array(i_list) + np.array(n_list)

# plotting
fig = plt.figure()
ax = fig.add_subplot()
ax.set_title("Equilibrium fraction of protein states against urea concentration")
ax.scatter(urea_list, np.array(d_list)/total_conc,label="D", marker='x')
ax.scatter(urea_list, np.array(i_list)/total_conc,label="I", marker='+')
ax.scatter(urea_list, np.array(n_list)/total_conc,label="N", marker='.')
ax.set_ylabel("Fraction of Species")
ax.set_xlabel("[Urea]/M")
ax.legend()

plt.show()