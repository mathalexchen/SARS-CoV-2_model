import os
import shutil

def write_sim_folder(sim_dir, ext):
    sim_num = get_sim_folder(sim_dir)
    src_dir = "."
    dest_dir = os.path.join(sim_dir, sim_num)
    os.makedirs(dest_dir)
    copy_dir(src_dir, dest_dir, ext)
    return sim_num

# copy all files with the given extension from src_dir to dest_dir        
def copy_dir(src_dir, dest_dir, ext):
    _, _, filenames = next(os.walk(src_dir))
    for file in filenames:
        if file.split(".")[-1] in ext:
            shutil.copyfile(os.path.join(src_dir, file), os.path.join(dest_dir, file))

def get_sim_folder(dir):
    sim_num = 1
    while os.path.exists(os.path.join(dir, str(sim_num))):
        sim_num += 1
    return str(sim_num)

def get_subdirs(dir):
    return [x[0] for x in os.walk(dir)]
        
