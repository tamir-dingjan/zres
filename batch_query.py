import glob
import subprocess
import sys
import os
import requests
import logging
import zres

import pandas as pd
import numpy as np

logging.basicConfig(
    filename="log",
    filemode='w',
    level=logging.DEBUG
)

# Read the Uniprot results 
if len(sys.argv) == 3:
    query = pd.read_csv(sys.argv[1], sep="\t")
    nterm_topol = sys.argv[2]
    if not (nterm_topol in ["in", "out"]):
        sys.exit("\nN-terminal topology must be either 'in' or 'out'. \nUsage: $> python batch_query query.csv (in|out)\n")
else:
    logging.debug(sys.argv)
    sys.exit("\nUsage: $> python batch_query.py query.csv (in|out)\n")


def fetch_structures():
    # Fetch the structure files
    local_structures = []
    for entry_id in query["Entry"].to_list():
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
                        if not "END" in line:
                            outfile.write(line+"\n")
                logging.info("Downloaded file: %s" % structfile)
                local_structures.append(structfile)
            except:
                logging.critical("Couldn't write file to disk: %s" % structfile)
    return local_structures


# TODO: Crop the structure files to the TM portion identified in the query results


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
            valid_results.append(analysis)
            valid_structures.append(i)

    if len(valid_results) != 0:
        results = pd.concat(valid_results)
        results.index = valid_structures

        # Save out summary
        results.to_csv('results.csv')
        logging.info("Results saved to 'results.csv'.")
    else:
        logging.info("No valid results from structure analysis.")

if __name__ == "__main__":
    local_structures = fetch_structures()
    run_ppm(local_structures)
    analyse_structures()
