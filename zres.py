from collections import defaultdict
import pandas as pd
import requests
import mdtraj as md
import numpy as np

class Zres:
    def __init__(self, pdb_code):
        self.pdb_code = pdb_code.lower()
        self.url = "https://opm-assets.storage.googleapis.com/pdb/" + self.pdb_code + ".pdb"
        self.structfile = pdb_code + ".pdb"
        self.zres = defaultdict(list)

    def run(self):
        # download the structure from OPM, removing END lines
        r = requests.get(self.url)
        # Check if OPM has this file
        if r.status_code == 404:
            print("OPM entry does not exist: %s" % self.pdb_code)
            return
        with open(self.structfile, 'w') as outfile:
            for line in r.text.split("\n"):
                if not "END" in line:
                    outfile.write(line+"\n")
        
        # load into MDTraj
        try:
            self.struct = md.load(self.structfile)
        except:
            print("ERROR: Couldn't load file %s" % self.structfile)

        # identify transmembrane region
        outer = self.struct.topology.select("resname DUM and name O")
        inner = self.struct.topology.select("resname DUM and name N")

        # Check if the membranes are aligned on the Z axis
        # OPM structures should be already aligned
        for border in [outer, inner]:
            if not (np.min(self.struct.xyz[0, border, 2]) == np.max(self.struct.xyz[0, border, 2])):
                print("Membrane boundary not aligned to Z-axis")

        self.outer_z = np.mean(self.struct.xyz[0, outer, 2])
        self.inner_z = np.mean(self.struct.xyz[0, inner, 2])
        self.memb_z = np.max((self.outer_z, self.inner_z)) - np.min((self.outer_z, self.inner_z))

        # Store Z-positions of C-alphas within the transmembrane region
        calphas = self.struct.topology.select("name CA")
        tmregion = calphas[np.where(np.logical_and((self.struct.xyz[0, calphas, 2] > np.min((self.outer_z, self.inner_z))),
                                                   (self.struct.xyz[0, calphas, 2] < np.max((self.outer_z, self.inner_z)))))]
        # Remap Z positions to range 0-to-1
        norm_z = (self.struct.xyz[0, tmregion, 2] - np.min((self.outer_z, self.inner_z))) / self.memb_z

        for atom, z in zip(tmregion, norm_z):
            self.zres[self.struct.topology.atom(atom).residue.name].append(z)

        # report as a pandas dataframe
        self.report = pd.DataFrame.from_dict({x:[self.zres[x]] for x in self.zres.keys()})
        return self.report
