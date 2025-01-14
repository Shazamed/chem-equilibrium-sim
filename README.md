# chem-equilibrium-sim
Calculate equilibrium concentration and changes in concentration over time of a 1D chemical system.
## Description
The script is able to read concentration and reaction data from a set of .txt files and calculates the rate equation for each of the reactions in the system. 

In this repository, there are two different scenarios saved which utilises the System class. The first is protein-folding-sim.py which calculates the equilibrium concentration of the system at different urea concentrations. The second is the oregonator.py, which shows how the concentration of a system evolves over time in a oscillatory reaction.
## Usage
Navigate to the directory where the repository is cloned.

Simply run either protein-folding-sim.py or oregonator.py with python:
```
python protein-folding-sim.py
```
```
python oregonator.py
```