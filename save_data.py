import os
import pickle
from kernel_model_URT import *

def save_obj(data, dir, name = "data"):
    os.makedirs(os.path.join(dir, name))
    with open(os.path.join(dir, name + '.pkl'), 'wb') as f:
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

# takes a self object from a class and saves each component separately
# done this way because of memory considerations
def save_obj_components(data, dir, name = "data"):
    os.makedirs(os.path.join(dir, name))
    for key in data.__dict__.keys():
        try:
            with open(os.path.join(dir, name, key + '.pkl'), 'wb') as f:
                pickle.dump(data.__dict__[key], f, pickle.HIGHEST_PROTOCOL)
        except:
            print("Error saving: ", key)
            
def load_obj(sim_folder, name = "data"):
    with open(os.path.join(sim_folder, name + '.pkl'), 'rb') as f:
        return pickle.load(f)