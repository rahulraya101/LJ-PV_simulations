import numpy as np
import espressomd
import math
required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

from espressomd import thermostat


#system units temperature in K, time in picoseconds ps, energy kT, length in nm.
#system parameters 150 particles, volume in nm**3, pressure in bar, density in nm**-3.

n_part = 150
volume  = 300.
box_l = volume**(1./3.)
system = espressomd.system(box_l=[box_l]*3)
density = (n_part / volume)
#intaraction parameters (non bounded intaractions, repulsive lennard jones)

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 2.5 * lj_sig
lj_cap = 20       #by artificially capping the forces, it is possible to simulate a system with overlaps.                  #by gradually raising the cap value, possible overlaps become unfavorable, and the system equilibrates to an overlap-free configuration.
time_step = 0.01

temperature = 250
gamma = 1.0
system.time_step = time_step
system.cell_system.skin = 0.4

#integration parameters

system.integrator.set_vv()
system.thermostat.set_langevin(kt=temperature, gamma=gamma, seed=42)


# warmup integration (with capped lj potential)
warm_steps = 100
warm_n_times = 30
# do the warmup until the particles have at least the distance min_dist
min_dist = 0.9*lj_sig

#integration

int_steps = 1000
int_n_times = 20

#interaction setup

system.non_bonded_inter[0,0].lennard_jones.set_params(
        epsilon=lj_eps, sigma=lj_sig,
        cutoff=lj_cut, shift="auto")
system.force_cap = lj_cap

print("lj-parameters:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# particle setup #either write a loop or use an (n_part x 3) array for positions. use np.random.random() to generate random numbers.
for i in range(n_part):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l)

system.analysis.dist_to(0)

print("simulate {} particles in a cubic simulation box of length {} at density {}."
      .format(n_part, box_l, density).strip())
print("Interactions:\n")
act_min_dist = system.analysis.min_dist()
print("Start with minimal distance {}".format(act_min_dist))

system.cell_system.max_num_cells = 2744

# Warmup Integration 

#open Observable file
obs_file = open("pylj_liquid.obs", "w")
obs_file.write("# Time\tE_tot\tE_kin\tE_pot\n")

print("""
Start warmup integration:
At maximum {} times {} steps
Stop if minimal distance is larger than {}
""".strip().format(warm_n_times, warm_steps, min_dist))


# Warmup Integration Loop
i = 0
while (i < warm_n_times and act_min_dist < min_dist):
    #print("run %d at time=%f " % (i, system.time))

    system.integrator.run(steps=warm_steps)
    # Warmup criterion
    act_min_dist = system.analysis.min_dist()
    i += 1   
# pressure = system.analysis.pressure()
#print(pressure)

# write parameter file

# polyBlockWrite "$name$ident.set" {box_l time_step skin} ""
set_file = open("pylj_liquid.set", "w")
set_file.write("box_l %s\ntime_step %s\nskin %s\n" %
               (box_l, system.time_step, system.cell_system.skin))

#Integration
for i in range(int_n_times):
   # print("run %d at time=%f " % (i, system.time))

    system.integrator.run(steps=int_steps)
    P = system.analysis.pressure()['total']
    obs_file.write("%f\n" % P)




#    pressure = (system.time, system.analysis.pressure()['total'])
#   obs_file.write("%s\n" % pressure)

#system.time = 0 #reset system timer
#pressure = np.zeros((int_n_times))


# write end configuration
xyz_file = open("pylj_liquid.xyz", "w")
xyz_file.write("{ time %f } \n { box_l %f }\n" % (system.time, box_l))
xyz_file.write("{ particles {id pos type} }")
for i in range(N_part):
    xyz_file.write("%s\n" % system.part[i].pos)
    # id & type not working yet

obs_file.close()
set_file.close()
xyz_file.close()

# terminate program
print("\nFinished.")
