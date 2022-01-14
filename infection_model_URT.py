import numpy as np
import time
from infection_model_cython_URT import infection_model
import pdb
import os
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as mticker
import matplotlib.ticker as mtick

class Simulation:
    def __init__(self, Dv, PCL_THICKNESS, MUCUS_THICKNESS, CELL_DIAM, CIRCUM, \
                infection_prob, NUM_VIRUS, v_adv,infectible_perc, MEMORY_CUTOFF, T, \
                virions_shed_per_day, num_days, latency_time, virus_at_interface, gen, restrict_to_gen, lung):
        self.Dv = Dv
        self.PCL_THICKNESS = PCL_THICKNESS
        self.MUCUS_THICKNESS = MUCUS_THICKNESS
        self.CELL_DIAM = CELL_DIAM
        self.CIRCUM = CIRCUM
        self.infection_prob = infection_prob
        self.NUM_VIRUS = NUM_VIRUS
        self.v_adv = v_adv
        self.infectible_perc = infectible_perc
        self.MEMORY_CUTOFF = MEMORY_CUTOFF
        self.T = T
        self.virions_shed_per_day = virions_shed_per_day
        self.num_days = num_days
        self.latency_time = latency_time
        self.virus_at_interface = virus_at_interface
        self.gen = gen
        self.restrict_to_gen = restrict_to_gen
        self.lung = lung
        
def run_sim(sim):
    a = time.time()

    pos_x, pos_y, pos_z, arrival_times = np.array([]), np.array([]), np.array([]), np.array([])
    
    infected_dict = {0: {0: 1}}
    total_time = 0
    
    # simulate in blocks of size MEMORY_CUTOFF to save storage space
    for i in range((sim.NUM_VIRUS - 1)//sim.MEMORY_CUTOFF + 1):
        num = ((sim.NUM_VIRUS - 1) % sim.MEMORY_CUTOFF + 1) if (i == (sim.NUM_VIRUS - 1)//sim.MEMORY_CUTOFF) else sim.MEMORY_CUTOFF
        pos_x_, pos_y_, pos_z_, arrival_times_, infected_dict = infection_model(sim, infected_dict, num)
        pos_x = np.append(pos_x, pos_x_)
        pos_y = np.append(pos_y, pos_y_)
        pos_z = np.append(pos_z, pos_z_)
        arrival_times = np.append(arrival_times, arrival_times_)

    assert len(pos_x) == sim.NUM_VIRUS    
    print("Num virus simulated = ", len(pos_x))
    print("Number of virus arrived = ", len(arrival_times))    
    if len(arrival_times) > 0:
        max_arrival_time, mean_arrival_time = np.max(arrival_times), np.mean(arrival_times)
        print("Max arrival time = ", max_arrival_time)
        print("Mean arrival time = ", mean_arrival_time)
    else: 
        max_arrival_time, mean_arrival_time = np.nan, np.nan
 
    return pos_x, pos_y, pos_z, arrival_times, infected_dict

def save_results(filename, pos_y, pos_z, arrival_times):
    with open(filename, 'wb') as f:
        np.save(f, np.array(pos_y))
        np.save(f, np.array(pos_z))
        np.save(f, np.array(arrival_times))

    
if __name__ == "__main__":
    Dv = 1.27 # um^2/s 
    PCL_THICKNESS = 7 # um
    MUCUS_THICKNESS = 50 # um
    CELL_DIAM = 4 # um
    CIRCUM = 150000 # diameter of nasal passage = 5 cm, so circumference = 50000\pi um
    infection_prob = 0.3
    T = 100000 # not a big deal most of the time, just choose something large enough so that they all infect
    MEMORY_CUTOFF = 10000
    NUM_VIRUS = 100000
    v_adv = 0 # um/s
    num_cells = 1001 # should be odd so that virions are shed from middle of first cell
    infectible_perc = 0.5
    virions_shed_per_day = 2000
    num_days = 2
    latency_time = 0.5
    
    sim = Simulation(Dv, PCL_THICKNESS, MUCUS_THICKNESS, CELL_DIAM, CIRCUM, \
                infection_prob, NUM_VIRUS, v_adv, infectible_perc, T, MEMORY_CUTOFF, num_cells, \
                virions_per_day, num_days, latency_time)
    a = time.time()
    pos_x, pos_y, pos_z, arrival_times, infected_cells = run_sim(sim)
    print(a - time.time())
    save_results('save_arr_v' + str(int(v_adv)) + '.npy', pos_y, pos_z, arrival_times)
            
    pdb.set_trace()