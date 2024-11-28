import numpy as np
import matplotlib.pyplot as plt
import pickle
from system_class import System
from timeit import default_timer as timer


# constants
TIME_STEP = 1e-6
MAX_TIME = 120
REACTION_PATH = "oregonator_inputs/oregonator-reactions.txt"

def int_divide_round(time, precision):
    '''Faster rounding to the decimal place indicated by precision'''
    product = float(10**precision)
    return int(time * product + 0.5)/product

class Cell(System):
    def __init__(self, state_dict, time_step, reaction_path):
        super().__init__(state_dict, time_step, reaction_path)


def load_initial_state():
    '''Read input file for initial concentrations and return the initial state dictionary eg.:
    {
    "X": 0.01,
    "Y": 0.5,
    "time" = 0.0
    }
    '''
    with open('oregonator_inputs/oregonator-concentrations.txt', 'r') as f:
        concentration_lines = f.read().splitlines()
        state_dict = {}
        for line in concentration_lines:
            state_dict[line.split(" ")[0]] = float(eval(line.split(" ")[1]))
    state_dict["time"] = 0
    return state_dict

# lists for plotting
time_list = []
x_list = []
y_list = []
z_list = []
state_list = []

try: # try to load the file if file exist
    state_list = pickle.load(open("oregonator-states.pickle", "rb"))
    last_state = state_list[-1]
    state = Cell(last_state, TIME_STEP, REACTION_PATH) # generate initial conditions from last state saved
    for state_dict in state_list:    # add saved states to list for plotting
        time_list.append(state_dict["time"])
        x_list.append(state_dict["X"])
        y_list.append(state_dict["Y"])
        z_list.append(state_dict["Z"])
    print("Data loaded")
except (OSError) as e:
    print("file does not exist")
    state = Cell(load_initial_state()) # if file does not exist generate from the concentration file


while state.time_elapsed < MAX_TIME:
    '''Loop to run simulation until MAX_TIME'''
    state.step()
    if state.time_elapsed == int_divide_round(state.time_elapsed,2): # save values every 0.01s of simulation time
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

# plotting
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