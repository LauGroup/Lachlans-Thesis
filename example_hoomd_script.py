import sys
sys.path.append("/home/chem/lada0992/hoomd-blue/3.8.1/pylibs")
import hoomd

import random
import itertools
import math
import gsd.hoomd
import numpy as np
import encapsulinv2

import os

### Example script used for a monte carlo diffusion simulation in hoomd


### first step is constructing the class for the encapsulin cage. Details on the method for constructing the model can be found in the 'encapsulinv2.py' comments
encshell = encapsulinv2.Encapsulin()
encshell.read_faces()
encshell.make_points(20,2,0.009469,0.1,0.1)

### HOOMD uses 'snapshots' which describe the state of the simulation system. Snapshots are used for initialisation and for recording frames of trajectories
### First the snapshot is initialised, along with the number of particles (3840 diffusing particles + 1 encapsulin + 1 LBT), and the simulation box (edge length 200 nm)
### Particle types are designated.

snapshot=gsd.hoomd.Snapshot()
snapshot.particles.N = 3842
snapshot.configuration.box = [200,200,200,0,0,0]
snapshot.particles.types=['LBT','Na','Enc']

### Particle ids are assigned to the simulation particles. 0 is the LBT/binding sphere, 1 is the diffusing particle (labelled Na by mistake, radius corresponds to Tb), 2 is the encapsulin

snapshot.particles.typeid = np.zeros(snapshot.particles.N,dtype='i')

snapshot.particles.typeid[0] = 0
snapshot.particles.typeid[1:3841] = 1
snapshot.particles.typeid[3841] = 2

### Charges assigned. For the simulations included in thesis no charge interactions were considered. In HPMC the charges don't affect the simulation unless a custom potential using the charges is provided

snapshot.particles.charge = np.zeros(snapshot.particles.N,dtype='f')

snapshot.particles.charge[0] = 0
snapshot.particles.charge[1:3841] = 0
snapshot.particles.charge[3841] = 0


###Particle positions are designated. Each particle gets a random position, after which the positions for the LBT and encapsulin are reset on the origin

snapshot.particles.position = np.zeros(shape=(snapshot.particles.N,3),dtype='f')
with np.nditer(snapshot.particles.position,['multi_index']) as positer:
    for i in positer:
        snapshot.particles.position[positer.multi_index] = random.randint(-100,100)

snapshot.particles.position[0,:] = [0.0,0.0,0.0]

snapshot.particles.position[3841,:] = [0.0,0.0,0.0]

### particles get diameters assigned. In HPMC the diameters designated here don't get used in the simulation, they're only for displaying in viewing software. 
### A radius of 0 for encapuslins means the coordinates/faces from the encapsulin class can be written to an .obj and used directly in gsd (viewing software) without any headaches with scaling

snapshot.particles.diameter = np.zeros(shape=(snapshot.particles.N),dtype='f')

snapshot.particles.diameter[0] = 0
snapshot.particles.diameter[1:3841] = 0.24
snapshot.particles.diameter[3841] = 0

### Timestep set to 0 and the snapshot is written to a file

snapshot.configuration.step = 0

with gsd.hoomd.open('230610pt11_snapshot.gsd',mode='wb') as f:
    f.append(snapshot)

### The simulation is created by designated the simulation device (CPU in this case, v3 might support HPMC with JIT on the GPU now but this wasn't the case when I started)
### The snapshot is loaded

cpu = hoomd.device.CPU()
sim = hoomd.Simulation(device=cpu,seed=13)
sim.create_state_from_gsd(filename = '230610pt11_snapshot.gsd')

### The HPMC integrator is designated. The integrator defines the geometry of all the particles in the simulation. In this case Polyhedron is used so that the encapsulin points and faces calcuated earlier can be used
### Spheres are made by designating a single vertex at the particle origin, a single face of all vertices, and designating a sweep radius that defines the actual particle radius used in HPMC
### The encapsulin shape is made based on the previously determined vertices and faces

mc = hoomd.hpmc.integrate.Polyhedron(translation_move_probability=1)

mc.shape['Enc'] = dict(vertices=encshell.shell_points,faces=encshell.faces)
mc.shape['Na'] = dict(vertices=[(0,0,0)],faces=[[0,0,0]],sweep_radius=0.12)
mc.shape['LBT'] = dict(vertices=[(0,0,0)],faces=[[0,0,0]])

### The move probabilities are assigned. d corresponds to translational move step lengths, a to rotational moves. Only the diffusing particles are allowed to move in the simulation, and rotations are turned off
### due to the particles being spherical

mc.d['Na'] = 1.09
mc.a['Na'] = 0
mc.d['Enc'],mc.a['Enc'],mc.d['LBT'],mc.a['LBT'] = 0,0,0,0

### Interaction potential between particles are defined. No electrostatic interactions are included (holdover from when I was trying to get them to work)
### The potential is a pairwise interaction between LBT/interaction sphere (type id = 0) and the diffusing particle (type id  = 1). If the sum of the ids of two particles = 1 (ie one is diffusing and one is lbt),
### an energy of -40 is returned, which is scaled based on the approximate binding energy of the LBT-Tb interaction with Kd 57 nM. Essentially its to make binding irreversible.

### The potential is square-well for simplicity of counting particles, and because its assumed a binding interaction between a peptide and ion in solution (where electrostatic interactions are quite short range) would
### be similar enough.


electrostatic = """if (type_i == 0 and type_j == 1 or type_i == 1 and type_j == 0) {
                            return -(40.0);
                            }
                        else
                            return 0.0;
                """
patch = hoomd.hpmc.pair.user.CPPPotential(r_cut=10, code=electrostatic,param_array=[])
mc.pair_potential = patch

### the integrator for the simulation is set to the one created above, including particle interactions

sim.operations.integrator = mc

### simulation is initialised so energy can be logged

sim.run(0)

class Energy():

    def __init__(self,sim):
            self.sim = sim

    @property
    def get_energy(self):
        return (sim.operations.integrator.pair_potential.energy)

logger = hoomd.logging.Logger(['scalar'])
patch_e = Energy(sim)

logger[('Energy','patch_energy')] = (patch_e,'get_energy','scalar')
logger.add(sim,quantities=['timestep','walltime'])

### output file is designated

gsd_writer = hoomd.write.GSD(filename='230610pt11_traj.gsd',
                             trigger=hoomd.trigger.Periodic(100),
                             mode='ab')
sim.operations.writers.append(gsd_writer)
gsd_writer.log = logger

print(sim.operations.integrator.pair_potential.energy)


### Simulation is run for 100000 steps

sim.run(100000)
