import numpy as np

def int_divide_round(time, precision):
    product = float(10**precision)
    return int(time * product + 0.5)/product

class System:
    def __init__(self, state_dict, time_step, reaction_path):
        self.time_elapsed = state_dict["time"]
        self.time_step = time_step
        self.reaction_path = reaction_path
        self.add_reactions()
        self.add_species()

        # initialize matrices
        self.conc = np.zeros(len(self.species_list))
        self.del_conc = np.zeros(len(self.species_list))

        self.initialize_concentrations(state_dict)
        self.generate_exponential_matrix()
        self.generate_rate_constants_matrix()


    def add_reactions(self):
        '''Generate list of reactions and rate constants from input'''
        with open(f"{self.reaction_path}", 'r') as f: 
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
                species_idx = self.species_list.index(reagent) # get index of reagent
                self.powers[species_idx, reaction_idx] += 1 # add one to the power term 
                self.coeff[species_idx, reaction_idx] -= 1
            for product in reaction[1]:
                species_idx = self.species_list.index(product) # get index of product
                self.coeff[species_idx, reaction_idx] += 1

    def initialize_concentrations(self, state_dict):
        for species_idx, species in enumerate(self.species_list): # initialize concentrations based off of the same index in the np array as the list
            self.conc[species_idx] = state_dict[species]
    
    def generate_rate_constants_matrix(self):
        self.k = np.array(self.rate_constant_list)

    def calc_rate(self):
        '''generates the rate from matrix operations
        np.array([self.conc]*len(self.reactions_list)).T gives an array where
        rows are the concentrations and columns are the reactions'''
        self.rate = np.prod(np.array([self.conc]*len(self.reactions_list)).T**self.powers, axis=0)*self.k # matrix before collapsing: rows = concentration, columns = reaction

    def calc_del_conc(self):
        self.del_conc = np.sum(np.array([self.rate]*len(self.species_list))*self.coeff,axis=1)*self.time_step

    def calc_conc(self):
        self.conc += self.del_conc

    def step(self):
        self.calc_rate()
        self.calc_del_conc()
        self.calc_conc()

        self.time_elapsed += self.time_step
        self.time_elapsed = int_divide_round(self.time_elapsed,6) # round the time elapsed to 6 d.p.

        
    def save(self):
        state_dict = {}
        for species_idx, species in enumerate(self.species_list):
            state_dict[species] = self.conc[species_idx]
        state_dict["time"] = self.time_elapsed
        return state_dict
