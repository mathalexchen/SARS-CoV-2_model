import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
import matplotlib.ticker as mtick
import numpy as np
import os
from copy import copy
import pdb
import pandas as pd

plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2

# Usage: summarize.print_num_virus(self.cell_li)
def print_num_virus(cell_li):
    total_virions = sum([cell.num_virus for cell in cell_li])
    print("Total virus =", total_virions)
    return total_virions

# Usage: summarize.num_cells_infected_by_time(self.cell_li, t)        
def num_cells_infected_by_time(cell_li, t):
    count = 0
    for cell in cell_li:
        if cell.t < t:
            count += 1
    return count
                        
# Usage: summarize.find_infectible_perc(self.inf_model.infected_dict)
def find_infectible_perc(infected_dict):
    count1 = 0
    count2 = 0
    for x in infected_dict:
        for y in infected_dict[x]:
            if infected_dict[x][y] == 1:
                count1 += 1
            else:
                count2 += 1
    print("Num infected cells =", count1)
    print("Num uninfected cells =", count2)

def output_to_csv(virus_prod_times, fname):
    df = pd.DataFrame(virus_prod_times)
    df.to_csv(fname)
    
# Usage: summarize.virus_production_times(self.cell_li, self.sim_folder, self.sim, "viral_load_adv")
def virus_production_times(cell_li, sim_folder, sim, name = "viral_load_adv"):
    virus_prod_times = np.zeros(int(sim.num_days*1440) + 1)
    virus_prod_times = virus_prod_times[720:] # cut off first 0.5 days
    virus_per_min = sim.virions_shed_per_day/(24*60)
    for cell in cell_li:
        start = int(cell.t/60 + sim.latency_time*24*60) - 720 # cell.t is in s, sim.latency_time in days
        finish = int(sim.num_days*24*60) - 720 + 1
        virus_prod_times[start:finish] = virus_prod_times[start:finish] + virus_per_min*np.arange(1,finish-start+1)

    fig, ax = plt.subplots()
    ax.plot(np.arange(finish)/1440, virus_prod_times)
    ax.ticklabel_format(axis = "y", scilimits = (0, 0), style = "sci")
    plt.xticks(np.arange(0, sim.num_days, 0.5))

    plt.xlim((0, sim.num_days - 0.5))
    plt.ylim((0, np.max(virus_prod_times)))
    plt.title("Total viral load over time")
    plt.xlabel("Time (days)")
    plt.ylabel("Number of virions")
    plt.savefig(os.path.join(sim_folder, "viral_load_adv" + str(sim.v_adv) + ".png"))
    output_to_csv(virus_prod_times, os.path.join(sim_folder, "total_viral_load_over_time.csv"))

# Usage: summarize.free_virions_over_time(self.free_virion_dist, self.sim_folder, self.sim, "free_virions")
def free_virions_over_time(free_virion_dist, sim_folder, sim, name = "free_virions"):
    plt.clf()
    free_virion_dist = free_virion_dist[720:] # cut off first 0.5 days
    free_virion_dist = free_virion_dist[:int(len(free_virion_dist)//(1440/2)*1440/2) + 1]
    plt.plot(np.arange(len(free_virion_dist))/1440, free_virion_dist)
    plt.xlabel("Time (days)")
    plt.ylabel("Number of free virions")
    plt.title("Free virions over time")
    plt.xticks(np.arange(0, sim.num_days - 0.5 + 0.1, 0.5))
    plt.ticklabel_format(axis = "y", style = "sci", scilimits = (-3, 3))
    plt.tight_layout()
    plt.savefig(os.path.join(sim_folder, name + ".png"))

def flux_out(exit_times, sim_folder, name = "flux"):
    np.savetxt(os.path.join(sim_folder, name + ".txt"), exit_times/(3600*24), fmt = "%10.5f")
    
def get_gen_len(sim):
    df = pd.read_csv("C:/Users/charg/Desktop/Jan 2021/lung_parameters.csv", sep = "\t")
    length = float(df[df["Generation"] == sim.gen]["Length (cm)"])*10000
    max_x = length//4
    return max_x

# Usage: summarize.infection_map_adv(self.infected_coord_dict, self.sim_folder, self.sim, "default", "infection_map_adv")
def infection_map_adv(infected_coord_dict, sim_folder, sim, mode, name = "infection_map_adv", plot_only_gen = True):
    min_x, max_x, min_y, max_y = np.inf, -np.inf, np.inf, -np.inf
    for pos in infected_coord_dict:
        min_x, max_x, min_y, max_y = min(min_x, pos[0]), max(max_x, pos[0]), min(min_y, pos[1]), max(max_y, pos[1])

    if plot_only_gen:
        max_x = get_gen_len(sim)
        
    min_x = int(min_x)
    max_x = int(max_x) + 1
    # For y (horizontal in plot), start the first infected cell in the center.
    max_y = int(abs(np.maximum(-min_y, max_y + 1)))
    min_y = -max_y
    
    map_x = max_x - min_x + 1
    map_y = max_y - min_y + 1
    if mode == "default":
        map_img = np.zeros((max_x - min_x + 1, max_y - min_y + 1))
        zero_loc = (-min_x, -min_y)
    elif mode == "compressed":
        map_img = np.zeros((int((max_x - min_x)/sim.v_adv) + 1, max_y - min_y + 1))
        zero_loc = (int(-min_x/sim.v_adv), -min_y)
    else: # log
        map_img = np.zeros((max_y - min_y + 1, max_y - min_y + 1))
        scale = (max_y - min_y)/np.log10(max_x - min_x + 1)
        
    for pos in infected_coord_dict:
        infected_time = infected_coord_dict[pos]/(24*3600)
        climit = round(2*sim.num_days)/2
        if infected_time <= climit:
            if 1:
            #try:
                if mode == "default":
                    cur_pos = (int(pos[0] - min_x), int(pos[1] - min_y))
                elif mode == "compressed":
                    cur_pos = (int((pos[0] - min_x)/sim.v_adv), int(pos[1] - min_y))
                else:
                    cur_pos = (int(np.log10(pos[0] - min_x + 1)*scale), int(pos[1] - min_y))

                if cur_pos[0] < map_img.shape[0]:
                    cur_val = map_img[cur_pos]
                    if cur_val > 0:
                        map_img[cur_pos] = np.minimum(infected_time - 0.5, cur_val)
                    else:
                        map_img[cur_pos] = infected_time - 0.5

    map_img = map_img/(map_img > 0)
    fig, ax = plt.subplots()
    newcmp = copy(plt.get_cmap('coolwarm'))
    newcmp.set_bad('black')
    
    if sim.gen != "nasal":
        pos = ax.imshow(map_img, cmap = newcmp, interpolation = "none", origin = "lower", vmin = 0, vmax = sim.num_days - 0.5)
    else:
        pos = ax.imshow(map_img, cmap = newcmp, interpolation = "none", vmin = 0, vmax = sim.num_days - 0.5)

    cbar = fig.colorbar(pos)
    cbar.set_label("Infection time (days)")
    cbar.set_ticks(np.arange(0, sim.num_days - 0.5, 0.25))
    ax.axis("off")
    plt.savefig(os.path.join(sim_folder, "infection_map_adv" + str(sim.v_adv) + "_" + mode + ".png"))
    np.savetxt(os.path.join(sim_folder, "infection_map_adv" + str(sim.v_adv) + "_" + mode + ".txt"), map_img, fmt = "%10.5f")
    return (map_x, map_y), map_img.shape, zero_loc

# Usage: summarize.write_summary_data(self.cell_li, total_virions, map, self.sim_folder, self.sim, "summary_data")    
def write_summary_data(cell_li, total_virions, map, map_img_dim, zero_loc, sim_folder, sim, name = "summary_data"):
    with open(os.path.join(sim_folder, name + ".txt"), "w") as f:
        f.writelines("Generation: " + str(sim.gen) + "\n")
        f.writelines("Advection: " + str(sim.v_adv) + "\n")
        f.writelines("Cell dimensions: " + str(map[0]) + " x " + str(map[1]) + " (# cells)" + "\n")
        f.writelines("Map dimensions: " + str(map_img_dim[0]) + " x " + str(map_img_dim[1]) + "\n")
        f.writelines("Origin location: " + str(zero_loc[0]) + " x " + str(zero_loc[1]) + "\n")
        f.writelines("Total virus: " + str(total_virions) + "\n")
        f.writelines("Total cells infected: " + str(len(cell_li)) + "\n")
        for i in range(1, int(2*sim.num_days) + 1):
            f.writelines("Total cells infected by time t = " + str((i - 1)/2) + " days: " + \
                        str(num_cells_infected_by_time(cell_li, i*12*3600)) + "\n")
        