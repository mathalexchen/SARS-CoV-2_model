import matplotlib.pyplot as plt
from save_data import *
import os

plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2

gens = ["Nasal", "Gen0", "Gen5", "Gen10", "Gen15"]
dir = "Simulations"
free_virion_data = []
for gen in gens:
    sim_folder = os.path.join(dir, gen)
    data_folder = os.path.join(sim_folder, "data")

    infected_coord_dict = load_obj(data_folder, name = "infected_coord_dict")
    sim = load_obj(data_folder, name = "sim")
    cell_li = load_obj(data_folder, name = "cell_li")
    free_virion_dist = load_obj(data_folder, name = "free_virion_dist")
    free_virion_data.append(free_virion_dist[720:2*1440 + 1])
    a = time.time()
    total_virions = summarize.print_num_virus(cell_li)
    summarize.free_virions_over_time(free_virion_dist, sim_folder, sim, "free_virions")

free_virion_data = pd.DataFrame(free_virion_data).T
free_virion_data.columns = ["Nasal", "Gen 0", "Gen 5", "Gen 10", "Gen 15"]
free_virion_data.to_csv("free_virion_data.csv")
pdb.set_trace()
