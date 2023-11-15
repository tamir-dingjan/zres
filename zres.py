from collections import defaultdict
import pandas as pd
import numpy as np
import MDAnalysis as mda
import os
import logging

class Zres:
    def __init__(self, structfile):
        
        if not os.path.isfile(structfile):
            logging.warning("Couldn't access file: %s" % structfile)
            return
        
        self.structfile = structfile
        self.zres = defaultdict(list)
    
    def run(self):
        # load structure
        try:
            struct = mda.Universe(self.structfile)
        except:
            logging.error("Couldn't load structure file %s" % self.structfile)
            return

        # identify transmembrane region
        outer = struct.select_atoms("resname DUM and name O and around 5 protein")
        inner = struct.select_atoms("resname DUM and name N and around 5 protein")

        # Catch cases without both membranes defined in OPM
        if len(outer) == 0 or len(inner) == 0:
            logging.warning("Structure file missing membrane information: %s" % self.structfile)
            return

        outer_z = np.mean(struct.coord[outer.indices][:,2])
        inner_z = np.mean(struct.coord[inner.indices][:,2])
        zmax = np.max((outer_z, inner_z))
        zmin = np.min((outer_z, inner_z))

        # Add border region buffers of 10 Angstrom?
        buffer_size = 10
        #self.outer_z = self.outer_z + buffer_size
        #self.inner_z = self.inner_z - buffer_size

        # Store Z-positions of C-alphas within the transmembrane region plus buffer
        
        tmregion = struct.select_atoms(f"protein and name CA and prop z <= {zmax + buffer_size} and prop z >= {zmin - buffer_size}")

        # tmregion = calphas[np.where(np.logical_and((self.struct.xyz[0, calphas, 2] > np.min((self.outer_z, self.inner_z))),
        #                                            (self.struct.xyz[0, calphas, 2] < np.max((self.outer_z, self.inner_z)))))]
        
        # Shift Z positions relative to membrane center
        z_relative_to_center = struct.coord[tmregion.indices][:,2] - ((zmax + zmin)/2)

        for atom, z in zip(tmregion, z_relative_to_center):
            self.zres[atom.resname].append(z)

        # report as a pandas dataframe
        return pd.DataFrame(dict([(k,pd.Series(v)) for k,v in self.zres.items()])).melt(var_name="residue").dropna() 
