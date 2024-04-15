import logging
import zres
import os
import sys
import argparse
import MDAnalysis as mda
import pandas as pd


logging.basicConfig(filename="log-aligned_analysis", filemode="w", level=logging.DEBUG)

# Argument parsing
ap = argparse.ArgumentParser()

ap.add_argument(
    "-a", "--alignment", required=True, help="FASTA-formatted alignment file"
)
ap.add_argument(
    "-s",
    "--select",
    required=True,
    help="residue position in the alignment to begin analysis from",
)

args = ap.parse_args()
logging.debug(vars(args))

if not os.path.isfile(args.alignment):
    sys.exit("Couldn't find alignment file: ", args.alignment)

try:
    select_from_residue = int(args.select)
except:
    sys.exit("Couldn't interpret selection as a residue number: ", args.select)


# Load in alignment file and cropping position
residue_threshold_lookup = {}
with open(args.alignment, "r") as alignment:
    id = ""
    read_associated_sequence = False
    for line in alignment.readlines():
        if ">" in line:
            if "|" in line:
                read_associated_sequence = True
                id = line.split("|")[1]
            continue
        else:
            if not read_associated_sequence:
                continue
            # Get the residue number corresponding to the given slignment position
            # AlphaFold structures begin from resid 1
            res_count = 0
            for seq_i, seq in enumerate(line):
                if seq_i == select_from_residue:
                    residue_threshold_lookup[id] = res_count
                    break
                if seq != "-":
                    res_count += 1
            read_associated_sequence = False


# Collect structures and remove any not found on disk
structures = [x for x in residue_threshold_lookup.keys()]

logging.debug(f"Structures collected from the alignment: {structures}")

for i in structures:
    if not os.path.isfile(i + "out.pdb"):
        structures.remove(i)

logging.debug(f"Found in total {len(structures)} aligned structures.")

valid_results = []

membrane_positions = []

# For each structure
for i in structures:

    # Load the structurefile
    structfile = i + "out.pdb"

    if not os.path.isfile(structfile):
        logging.debug("Structure file not found: ", structfile)
        continue

    membrane_boundaries = zres.Zres(structfile).get_membrane_boundaries()
    if not membrane_boundaries is None:
        membrane_positions.append(membrane_boundaries)

    struct = mda.Universe(structfile)

    # Get this structure's selection from which to begin the analysis

    residue_threshold = residue_threshold_lookup[i]

    # Set the desired portion of the in-memory structure to analyse using the "segid" attribute using the alignment
    selected_segment_label = "B"
    selected_segment = struct.add_Segment(segid=selected_segment_label)
    struct.residues[residue_threshold:].segments = selected_segment

    # Analyse residue Z-positions
    analysis = zres.Zres(structfile).run_on_modified_structure(
        struct, segid=selected_segment_label, buffer=0
    )
    if not analysis is None:
        analysis.index = [structfile] * len(analysis)
        valid_results.append(analysis)

if len(valid_results) != 0:
    results = pd.concat(valid_results)
    membrane = pd.concat(membrane_positions)

    # Save out summary
    results.to_csv("results.csv")
    logging.info("Results saved to 'results.csv'.")

    membrane.to_csv("membrane.csv")
    logging.info("Membrane positions saved to 'membrane.csv'.")
else:
    logging.info("No valid results from structure analysis.")
