import numpy as np
import time
try:
    from infection_model_cython_URT import infection_model
    from infection_model_URT import *
except:
    print("You may need to run setup on setup_infection_model_URT. Use 64-bit python in order to not run into memory issues.")
import pdb
import bisect
import itertools
import pandas as pd

# internal modules
import save_data
import summarize_data_and_plot as summarize
from copy_sim_dir import *
from time_profiler import Timer


class InfectionModelSimulation:
    def __init__(self, sim):
        a = time.time()
        _, self.pos_y, self.pos_z, self.arrival_times, self.infected_dict = run_sim(sim)
        print("Time =", time.time()-a)

class InfectedCell:
    def __init__(self, sim, pos, t, inf_model):
        self.t = t
        self.pos = pos
        self.CELL_DIAM = sim.CELL_DIAM 
        self.length = sim.lung.length[sim.gen]
        self.v_adv = sim.v_adv
        
        self.get_virus_prod_times(sim, inf_model)
        self.get_sample(inf_model)
                
    def get_virus_prod_times(self, sim, inf_model):
        virion_prod_time = 24*3600/sim.virions_shed_per_day
        self.num_virus = max(0, int(round((sim.num_days*24*3600 - self.t - sim.latency_time*24*3600)/virion_prod_time)))
        self.prod_times = np.round(np.arange(self.t + sim.latency_time*24*3600, sim.num_days*24*3600, virion_prod_time),2)
        self.prod_times = self.prod_times[:self.num_virus]
        
    def get_sample(self, inf_model):
        choice = np.random.choice(range(len(inf_model.pos_y)), self.num_virus, replace = False)
        self.pos_y = inf_model.pos_y[choice]
        self.pos_z = inf_model.pos_z[choice]
        self.arrival_times = inf_model.arrival_times[choice] + self.prod_times
        
        idx = np.where((self.pos[0] + self.pos_y)*self.CELL_DIAM <= self.length)
        exit_idx = list(set(range(self.num_virus)).difference(idx[0]))
        self.exit_times = self.prod_times + (self.length - (self.pos[0]*self.CELL_DIAM))/self.v_adv
        self.exit_times = self.exit_times[exit_idx]
        self.prod_times_exit = self.prod_times[exit_idx]            

        self.pos_y = self.pos_y[idx]
        self.pos_z = self.pos_z[idx]
        self.arrival_times = self.arrival_times[idx]
        self.prod_times = self.prod_times[idx]
        self.flux_out = self.num_virus - len(idx[0])
        
        del choice

    # clears up pos_y, pos_z, prod_times, arrival_times, the largest storage
    def reduce_storage_space(self):
        del self.pos_y
        del self.pos_z
        del self.prod_times
        del self.arrival_times
            
class DistanceRanking:
    def __init__(self, dim):
        li = list(itertools.product(range(-dim, dim + 1), range(-dim, dim + 1)))
        sorted_li = sorted(li, key = lambda x: x[0]**2 + x[1]**2)
        self.dist_li = [[]]
        cur_dist = 0
        cur_element = 0
        for element in sorted_li:
            if element[0]**2 + element[1]**2 != cur_dist:
                self.dist_li.append([])
                cur_element += 1
                cur_dist = element[0]**2 + element[1]**2
            self.dist_li[cur_element].append(element)

        self.dist_li = self.dist_li[1:] # cut out (0, 0)
        self.NUM_ROWS = 1000
        self.row_num = np.zeros(len(self.dist_li)) # how many rows of each distance rank have been used up
        self.dist_group = 0
        self.count = 0 # how many entries of the given dist_li[i] used up
        self.distances = []
        
        for i in range(len(self.dist_li)):
            arr = np.tile(np.arange(len(self.dist_li[i])), (self.NUM_ROWS, 1))
            self.distances.append(arr)
    
    def get_next_element(self):
        if self.count == len(self.dist_li[self.dist_group]): # go to next dist group
            self.row_num[self.dist_group] += 1
            self.dist_group += 1
            self.count = 0
            
        if self.row_num[self.dist_group] % self.NUM_ROWS == 0:
            arr = self.permute_each_row(self.distances[self.dist_group])

        self.element = self.dist_li[self.dist_group][self.count]
        self.count += 1
        return self.element

    def update_infected_cell_found(self):
        self.row_num[self.dist_group] += 1
        self.count = 0
        self.dist_group = 0

    def permute_each_row(self, s):
        return s[np.arange(len(s))[:,None], np.random.randn(*s.shape).argsort(axis=1)]
    
class FullModel:
    def __init__(self, sim, sim_folder):
        self.sim = sim
        self.sim_folder = sim_folder
        
        self.inf_model = InfectionModelSimulation(sim)
        self.dist_rank = DistanceRanking(3)
        
        self.infection_times = [0]
        self.infected_coord_dict = {(0,0): 0}
        self.infected_time_dict = {0: [(0,0)]}
        self.t = 0
        self.cell_li = []
        self.free_virion_dist = np.zeros(int(sim.num_days*1440))
        self.randu_count = 0
                
        num = 0
        a = time.time()
        while num < len(self.infection_times):
            try:
                coords = self.infected_time_dict[self.infection_times[num]]
            except:
                pdb.set_trace()

            for i in range(len(coords)):
                if not sim.restrict_to_gen or coords[i][0]*sim.CELL_DIAM <= self.sim.lung.length[self.sim.gen]:
                    cell = InfectedCell(sim, coords[i], self.infection_times[num], self.inf_model)
                    self.update_infected_dict(cell) 
                    self.update_free_virion_dist(cell)
                    cell.reduce_storage_space()
                    self.cell_li.append(cell)
                    if len(self.cell_li) % 1000 == 0:
                        print("At cell: ", len(self.cell_li), cell.t/3600/24)
            num += len(coords)
        print("Time =", time.time() - a)   
        self.calculate_flux_out()
        self.summarize_and_plot()

    def update_free_virion_dist(self, cell):
        for i in range(len(cell.arrival_times)):
            start = int(cell.prod_times[i]//60)
            finish = int(cell.arrival_times[i]//60)
            self.free_virion_dist[start:finish] = self.free_virion_dist[start:finish] + 1
            
        for i in range(len(cell.prod_times_exit)):
            start = int(cell.prod_times_exit[i]//60)
            finish = int(cell.exit_times[i]//60)
            self.free_virion_dist[start:finish] = self.free_virion_dist[start:finish] + 1

    def calculate_flux_out(self):
        self.flux_out = 0
        self.exit_times = np.array([])
        with open(os.path.join(self.sim_folder, "flux.txt"), "ab") as f:
            for cell in self.cell_li:
                self.flux_out += cell.flux_out
                np.savetxt(f, cell.exit_times/(3600*24) - 0.5, fmt = "%10.5f")
            
    def update_infected_dict(self, cell):
        for i in range(len(cell.pos_y)):
            candidate_y = cell.pos[0] + cell.pos_y[i]
            candidate_z = cell.pos[1] + cell.pos_z[i]
            candidate_y, candidate_z = self.adjust_pos((cell.pos[0] + cell.pos_y[i], cell.pos[1] + cell.pos_z[i]))
            if (candidate_y, candidate_z) in self.infected_coord_dict:
                if self.infected_coord_dict[(candidate_y, candidate_z)] > cell.arrival_times[i]:
                    if len(self.infected_time_dict[self.infected_coord_dict[(candidate_y, candidate_z)]]) == 1:
                        self.infected_time_dict.pop(self.infected_coord_dict[(candidate_y, candidate_z)])
                    else:
                        self.infected_time_dict[self.infected_coord_dict[(candidate_y, candidate_z)]].remove((candidate_y, candidate_z))
                    del self.infection_times[bisect.bisect_left(self.infection_times, self.infected_coord_dict[(candidate_y, candidate_z)])]
                    self.add_time_to_dict(candidate_y, candidate_z, cell.arrival_times[i])
            else:
                self.add_time_to_dict(candidate_y, candidate_z, cell.arrival_times[i])
            
    def add_time_to_dict(self, pos_y, pos_z, arrival_time):
        bisect.insort(self.infection_times, arrival_time)
        self.infected_coord_dict[(pos_y, pos_z)] = arrival_time
        if arrival_time in self.infected_time_dict:
            self.infected_time_dict[arrival_time].append((pos_y, pos_z))
        else:
            self.infected_time_dict[arrival_time] = [(pos_y, pos_z)]

    def adjust_pos(self, pos):
        new_pos = pos
        while True:
            self.simulate_cell(new_pos)
            if self.inf_model.infected_dict[new_pos[0]][new_pos[1]]:
                self.dist_rank.update_infected_cell_found()
                return new_pos
            new_pos = np.array(pos) + self.dist_rank.get_next_element()
    
    def simulate_cell(self, pos):
        if pos[0] in self.inf_model.infected_dict and pos[1] in self.inf_model.infected_dict[pos[0]]:
            return 
            
        if pos[0] not in self.inf_model.infected_dict:
            self.inf_model.infected_dict[pos[0]] = {}

        if self.randu_count % 10000 == 0:
            self.randu = np.random.rand(10000)
            
        if self.randu[self.randu_count % 10000] < self.sim.infectible_perc:
            self.inf_model.infected_dict[pos[0]][pos[1]] = 1
        else:
            self.inf_model.infected_dict[pos[0]][pos[1]] = 0
        self.randu_count += 1

    def summarize_and_plot(self):
        a = time.time()
        total_virions = summarize.print_num_virus(self.cell_li)
        try:
            save_data.save_obj_components(self, self.sim_folder, "data")
        except:
            print("Error saving data!")
            pdb.set_trace()

        summarize.free_virions_over_time(self.free_virion_dist, self.sim_folder, self.sim, "free_virions")

        if self.sim.v_adv != 0:
            try:
                if self.sim.gen in ["10", "15"]:
                    map, map_img_dim, zero_loc = summarize.infection_map_adv(self.infected_coord_dict, self.sim_folder, self.sim, "default", "infection_map_adv")
                else:
                    map, map_img_dim, zero_loc = summarize.infection_map_adv(self.infected_coord_dict, self.sim_folder, self.sim, "compressed", "infection_map_adv")
            except:
                print("Error creating infection map!")
                pdb.set_trace()
        else:
            map, map_img_dim, zero_loc = summarize.infection_map_adv(self.infected_coord_dict, self.sim_folder, self.sim, "default", "infection_map_adv")

        try:
            summarize.write_summary_data(self.cell_li, total_virions, map, map_img_dim, zero_loc, self.sim_folder, self.sim, "summary_data")
        except:
            print("Error summarizing data!")
            pdb.set_trace()

        try:
            summarize.virus_production_times(self.cell_li, self.sim_folder, self.sim, "viral_load_adv")
        except:
            print("Error computing viral load!")
            pdb.set_trace()
                    
        print("Time =", time.time() - a)    

class LungParam:
    def __init__(self, param_csv):
        self.param_df = pd.read_csv(param_csv, sep = "\t")
        self.param_df = self.param_df.set_index("Generation")
        self.length = (self.param_df["Length (cm)"]*10000).to_dict()

def main():
    Dv = 1.27 # um^2/s 
    PCL_THICKNESS = 7 # um
    CELL_DIAM = 4 # um
    CIRCUM = 150000 # diameter of nasal passage = 5 cm, so circumference = 50000\pi um
    infection_prob = 0.3
    NUM_VIRUS = 10000
    infectible_perc = 0.5
    virions_per_day = 2000
    num_days = 2.02
    latency_time = 0.5
    virus_at_interface = 0
    restrict_to_gen = False
    
    MEMORY_CUTOFF = 10000 # simulate in smaller blocks
    MAX_SIM_TIME_kernel = 100000 # usually infects very quickly, so a large number can be chosen
    
    sim_dir = "Simulations"
    ext = ["py", "pyx"]
    
    gens_to_sim = ["nasal","0","5","10","15"]
    
    param_csv = "lung_parameters.csv"
    param_df = pd.read_csv(param_csv, sep = "\t")
    param_df = param_df[param_df["Generation"].isin(gens_to_sim)]
    lung = LungParam(param_csv)

    for par in range(len(param_df)):
        v_adv = param_df.iloc[par]["Advection (mm/min)"]*1000/60 # mm/min -> um/s
        MUCUS_THICKNESS = 0.5*(param_df.iloc[par]["ASL Height (lower, um)"] + param_df.iloc[par]["ASL Height (upper, um)"])
        gen = param_df.iloc[par]["Generation"]
        
        sim_num = write_sim_folder(sim_dir, ext)
        sim_folder = os.path.join(sim_dir, sim_num)
        
        sim = Simulation(Dv, PCL_THICKNESS, MUCUS_THICKNESS, CELL_DIAM, CIRCUM, \
                    infection_prob, NUM_VIRUS, v_adv, infectible_perc, MEMORY_CUTOFF, MAX_SIM_TIME_kernel, \
                    virions_per_day, num_days, latency_time, virus_at_interface, gen, restrict_to_gen, lung)
        
        FullModel(sim, sim_folder)

if __name__ == "__main__":
    main()
