import numpy as np
cimport numpy as np

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

DTYPEint = np.int
ctypedef np.int_t DTYPEint_t

cimport cython

@cython.boundscheck(False) 
@cython.wraparound(False)  

def get_cell_coord(y, z, cell_diam):
    cell_coord_y = int(round(y/cell_diam))
    cell_coord_z = int(round(z/cell_diam))
    return cell_coord_y, cell_coord_z

def simulate_infectivity(infected_dict, cell_coord_y, cell_coord_z, infectible_perc):
    if np.random.rand() < infectible_perc:
        infected_dict[cell_coord_y][cell_coord_z] = 1
    else:
        infected_dict[cell_coord_y][cell_coord_z] = 0

def add_infectible_entry(infected_dict, cell_coord_y, cell_coord_z, infectible_perc):
    if cell_coord_y not in infected_dict:
        infected_dict[cell_coord_y] = {}
    
    if cell_coord_z not in infected_dict[cell_coord_y]:
        simulate_infectivity(infected_dict, cell_coord_y, cell_coord_z, infectible_perc)

def infection_model(sim, infected_dict, DTYPEint_t num_virus):
    cdef DTYPE_t Dv = sim.Dv
    cdef DTYPE_t PCL_THICKNESS = sim.PCL_THICKNESS
    cdef DTYPE_t MUCUS_THICKNESS = sim.MUCUS_THICKNESS
    cdef DTYPE_t CELL_DIAM = sim.CELL_DIAM
    cdef DTYPE_t CIRCUM = sim.CIRCUM
    cdef DTYPE_t infectible_perc = sim.infectible_perc
    cdef DTYPE_t infection_prob = sim.infection_prob
    cdef DTYPEint_t T = sim.T
    cdef DTYPE_t v_adv = sim.v_adv
    cdef DTYPEint_t virus_at_interface = sim.virus_at_interface

    cdef int NUM_BLOCK = 100000 # number of random numbers in block to simulate
    cdef np.ndarray[DTYPE_t, ndim=2] randnums = np.random.randn(NUM_BLOCK, 3).astype(DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] randu = np.random.rand(NUM_BLOCK).astype(np.float32)
    
    cdef np.ndarray[DTYPEint_t, ndim=1] arrival_times = np.zeros(num_virus, dtype = DTYPEint)
    cdef np.ndarray[DTYPE_t, ndim=1] pos_x = np.zeros(num_virus, dtype = DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] pos_y = np.zeros(num_virus, dtype = DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] pos_z = np.zeros(num_virus, dtype = DTYPE)
    cdef int num = 0
    cdef int numu = 0
    cdef DTYPE_t x, y, z
    cdef int i, j
    cdef int cell_coord_y, cell_coord_z
        
    for i in range(num_virus):
        x = virus_at_interface
        y = 0
        z = 0
        for j in range(1, T + 1):
            if num == randnums.shape[0]: # generate new random numbers
                randnums = np.random.randn(randnums.shape[0], 3).astype(DTYPE)
                num = 0

            x += (2*Dv*1)**0.5*randnums[num,0]
            y += (2*Dv*1)**0.5*randnums[num,1]
            z += (2*Dv*1)**0.5*randnums[num,2]
            num += 1

            if x > MUCUS_THICKNESS + PCL_THICKNESS:
                x = 2*(MUCUS_THICKNESS + PCL_THICKNESS) - x
            
            # stopping condition is that cell is reached with infection prob and the cell is infectible
            if x < 0:
                cell_coord_y, cell_coord_z = get_cell_coord(y, z, CELL_DIAM)
                add_infectible_entry(infected_dict, cell_coord_y, cell_coord_z, infectible_perc)
                if infected_dict[cell_coord_y][cell_coord_z] > 0: 
                    numu += 1
                    if numu == len(randu): # generate new random numbers
                        randu = np.random.rand(NUM_BLOCK).astype(np.float32)
                        numu = 0
                        
                    if randu[numu] < infection_prob:
                        arrival_times[i] = j
                        break
            x = abs(x)
            
            if x > PCL_THICKNESS:
                y += v_adv # advection um/s

        pos_x[i] = x
        pos_y[i] = cell_coord_y
        pos_z[i] = cell_coord_z
    
    return pos_x, pos_y, pos_z, arrival_times, infected_dict