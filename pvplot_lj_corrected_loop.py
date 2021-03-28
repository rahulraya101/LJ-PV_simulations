import numpy as np
import matplotlib.pyplot as plt
import espressomd
from espressomd.io.writer import vtf
import math
required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)
import sys
from espressomd import thermostat

#system is isothermal-isochoric ensemble (NVT)
#system units temperature in K, time in picoseconds ps, energy kT, length in nm.
#system parameters 150 particles, volume in nm**3, pressure in bar, density in nm**-3.

N_part = 150
volume  = float(sys.argv[1]) #300. # nm3
box_l = volume**(1./3.) #/sigma_nm
system = espressomd.System(box_l=[box_l]*3)
Density = (N_part / volume)
#intaraction parameters (non bounded intaractions, repulsive lennard jones)

temperature = 250  #float(sys.argv[2])  #in kelvin
lj_eps = 119.8/temperature
lj_sig = 0.314   #in pystar for argon (from Sasha) eps/k_B(in K) = 111.84 and sig(in nm) = 0.3623  # (from book)sig(in nm) = 0.314, eps/k = 119.8K 
lj_cut = 2.5 * lj_sig
lj_cap = 20       #By artificially capping the forces, it is possible to simulate a system with overlaps.
#By gradually raising the cap value, possible overlaps become unfavorable, and the system equilibrates to an overlap-free configuration.
time_step = 0.01

gamma = 1.0
system.time_step = time_step
system.cell_system.skin = 0.4

#integration parameters

system.integrator.set_vv()
system.thermostat.set_langevin(kT = temperature/300, gamma=gamma, seed=42)   #300K is room temperature. using this here for considering given temperature values which is 250K

# warmup integration (with capped LJ potential)
warm_steps = 1000
warm_n_times = 50
# do the warmup until the particles have at least the distance min_dist
min_dist = 0.9*lj_sig

#integration setup

int_steps = 10000
int_n_times = 100

system.non_bonded_inter[0,0].lennard_jones.set_params(
        epsilon=lj_eps, sigma=lj_sig,
        cutoff=lj_cut, shift="auto")
system.force_cap = lj_cap

print("LJ-parameters:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# Particle setup #Either write a loop or use an (N_PART x 3) array for positions. Use np.random.random() to generate random numbers.
for i in range(N_part):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l)


print("Simulate {} particles in a cubic simulation box of length {} at density {}."
      .format(N_part, box_l, Density).strip())
print("Interactions:\n")
act_min_dist = system.analysis.min_dist()
print("Start with minimal distance {}".format(act_min_dist))

#open Observable file
obs_file = open("pvplot_lj.obs", "w")
obs_file.write("# Time\tE_tot\tE_kin\tE_pot\tP_tot\n")
obs_file1 = open("pvplot_lj_avg.obs", "w")
obs_file1.write("# Time\tp_avg\tp_std\te_avg\te_std\n")
print("""
Start warmup integration:
At maximum {} times {} steps
Stop if minimal distance is larger than {}
""".strip().format(warm_n_times, warm_steps, min_dist))


# Warmup Integration Loop
i = 0
while (i < warm_n_times and act_min_dist < min_dist):
    system.integrator.run(steps=warm_steps)
    act_min_dist = system.analysis.min_dist()
    i += 1   
#open set file
set_file = open("pvplot_lj.set", "w")
set_file.write("box_l %s\ntime_step %s\nskin %s\n" %
               (box_l, system.time_step, system.cell_system.skin))

#Integration
p_samples = []
kT = 1.38e-23 * temperature   #now in jolues
to_bars = kT / (1e5 * 1e-27 ) #pas to bar 1e5

times, p_tots, e_tots, e_kins, e_pots, = [], [], [], [], [], 

for i in range(int_n_times):
    good = True
    system.integrator.run(steps=int_steps)
    time = system.time
    energy = system.analysis.energy()
    e_tot = energy['total']
    e_kin = energy['kinetic']
    e_pot = energy['non_bonded']
    p_tot = system.analysis.pressure()['total']
   # obs_file.write("%f %f %f %f %f\n" % (time, e_tot , e_kin , e_pot , p_tot * to_bars,))
    

    if e_tot > 300: good = False
    if p_tot > 200: good = False
    if good: obs_file.write("%f %f %f %f %f\n" % (time, e_tot , e_kin , e_pot , p_tot * to_bars,))
    if good:
        times.append(time)
        p_tots.append(p_tot * to_bars)
        e_tots.append(e_tot)
        e_kins.append(e_kin)
        e_pots.append(e_pot)

p_avg = np.mean(p_tots)
e_avg = np.mean(e_tots)
p_std = np.std(p_tots)
e_std = np.std(e_tots)
obs_file1.write("%f %f %f %f\n" % (p_avg , p_std , e_avg , e_std))

xyz_file = open("pvplot_lj.xyz", "w")
xyz_file.write("{ time %f } \n { box_l %f }\n" % (system.time, box_l))
xyz_file.write("{ particles {id pos type} }\n")
for i in range(N_part):
    current_position_array = system.part[i].pos
    good_coords = []
    for coord in current_position_array:
        while(coord < 0 or coord > box_l):
            if coord > box_l: coord -= box_l
            if coord < 0: coord += box_l
        good_coords.append(coord)

    xyz_file.write("%s\n" % good_coords)

obs_file.close()
obs_file1.close()
set_file.close()
xyz_file.close()

plt.figure(figsize=(10, 6))
plt.subplot(211)
plt.plot(times, p_tots)
plt.legend()
#plt.xlabel('t' , size=16)
plt.ylabel('p' , size=16)
plt.title("p_tot-t plot")
plt.savefig("P_t_{}.pdf".format(volume))

#plt.figure(figsize=(10, 6))
plt.subplot(212)
plt.plot(times, e_tots)
plt.legend()
plt.xlabel('t' , size=16)
plt.ylabel('e' , size=16)
plt.title("e_tot-t plot")
plt.savefig("e_t_{}.pdf".format(volume))



plt.show()


# terminate program
print("\nFinished.")
