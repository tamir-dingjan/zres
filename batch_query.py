import glob
import subprocess
import sys
import os
import requests
import logging
import zres
import argparse

import pandas as pd
import numpy as np

logging.basicConfig(
    filename="log",
    filemode='w',
    level=logging.DEBUG
)

# Construct argument parser
ap = argparse.ArgumentParser()

ap.add_argument("-q", "--query", required=True, help="tab-seperated query file containing UniProt IDs")
ap.add_argument("-n", "--nterm", required=True, help="N-terminal topology. Must be either 'in' or 'out'")
ap.add_argument("--crop", required=True, action=argparse.BooleanOptionalAction)

args = ap.parse_args()
logging.debug(vars(args))

# Read the Uniprot results 
if os.path.isfile(args.query):
    query = pd.read_csv(args.query, sep="\t")
else:
    sys.exit("Couldn't find query file: ", args.query)

if args.nterm in ["in", "out"]:
    nterm_topol = args.nterm
else:
    sys.exit("\nN-terminal topology must be either 'in' or 'out'. \nUsage: $> python batch_query query.csv (in|out)\n")


def is_desired_tm_portion(x, tm_start, tm_end):
    if not "ATOM" in x:
        return True
    else:
        try:
            resnum = int(x[22:26])
        except:
            logging.warning("Couldn't read residue number. Bailing on this file.")
            return False
        
        if (resnum > tm_start) and (resnum < tm_end):
            return True


def fetch_structures():
    # Fetch the structure files
    local_structures = []
    tm_region_buffer = 10

    for entry_id, tm in zip(query["Entry"].to_list(), query["Transmembrane"].to_list()):
        if (tm == np.nan) or (pd.isna(tm)):
            logging.warning("No TM location information for entry: %s" % entry_id)
            continue
        
        try:
            locs = tm.split()[1].split(";")[0]
            tm_start = int(locs.split(".")[0]) - tm_region_buffer
            tm_end = int(locs.split(".")[-1]) + tm_region_buffer
        except:
            logging.warning("Couldn't process TM location for entry: %s" % entry_id)
            logging.warning("Entry: %s, Transmembrane field: %s" % (entry_id, tm))
            continue
        
        structfile = entry_id + ".pdb"
        if os.path.isfile(structfile):
            logging.info("Found file: %s" % structfile)
            local_structures.append(structfile)
            continue
        else:
            url = "https://alphafold.ebi.ac.uk/files/AF-" + entry_id + "-F1-model_v4.pdb"
            r = requests.get(url)
            if r.status_code == 404:
                logging.warning("Couldn't find entry at provider's URL: %s" %  url)
                continue    
            try:
                with open(structfile, 'w') as outfile:
                    for line in r.text.split("\n"):
                        if args.crop:
                            if (not "END" in line) and (is_desired_tm_portion(line, tm_start, tm_end)):
                                outfile.write(line+"\n")
                        else:
                            if (not "END" in line):
                                outfile.write(line+"\n")

                logging.info("Downloaded file: %s" % structfile)
                local_structures.append(structfile)
            except:
                logging.critical("Couldn't write file to disk: %s" % structfile)
    return local_structures


def run_ppm(local_structures):
    # Build the PPM input
    logging.info("Writing PPM input file.")
    with open("ppm.inp", 'w') as outfile:
        outfile.write("1\n")
        for structure in local_structures:
            outfile.write(" 0 ERm " + nterm_topol.ljust(3) + " " + structure + "\n")

    # Run PPM
    logging.info("Now running PPM...")
    # To make the STDIN and STDOUT piping simpler, use a shell wrapper
    with open("run_ppm.sh", 'w') as outfile:
        outfile.write("./immers<ppm.inp>ppm.out")
    subprocess.run(
        ["bash","run_ppm.sh"],
        cwd=os.getcwd()
    )
    logging.info("PPM complete.")


def analyse_structures():
    # Use the Zres module to analyse the membrane-embedded structures 
    result_structures = glob.glob(os.getcwd()+"/*out.pdb")
    logging.info("Found %i PPM result structures." % len(result_structures))

    valid_results = []
    valid_structures = []

    # Collect zres results
    for i in result_structures:
        analysis = zres.Zres(i).run()
        if not analysis is None:
            analysis.index = [i]*len(analysis)
            valid_results.append(analysis)
            valid_structures.append(i)

    if len(valid_results) != 0:
        results = pd.concat(valid_results)
        
        # Save out summary
        results.to_csv('results.csv')
        logging.info("Results saved to 'results.csv'.")
    else:
        logging.info("No valid results from structure analysis.")

if __name__ == "__main__":
    local_structures = fetch_structures()
    run_ppm(local_structures)
    analyse_structures()
