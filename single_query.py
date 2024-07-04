import sys
import os
import logging
import zres
import argparse

logging.basicConfig(
    filename="log",
    filemode='w',
    level=logging.DEBUG
)

# Construct argument parser
ap = argparse.ArgumentParser()

ap.add_argument("-f", "--file", required=True, help="PDB file to analyse. Must contain PPM-formatted bilayer dummy atoms.")

args = ap.parse_args()
logging.debug(vars(args))

def analyse_structure(file):
    
    if not os.path.isfile(file):
        sys.exit("Couldn't find PDB file: ", file)

    outfile = os.path.join(os.path.dirname(file), os.path.splitext(os.path.basename(file))[0])
        
    # Use the Zres module to analyse the membrane-embedded structures 
    logging.info("Input file specified: %s " % file)
    
    analysis = zres.Zres(file).run()
    
    if analysis is None:
        logging.info("No valid results from structure analysis.")
        sys.exit()
    else:
        analysis.to_csv(outfile+".csv")
        zres.Zres(file).get_membrane_boundaries().to_csv(outfile+"_membrane.csv")
        logging.info("Results saved to {%s}" % outfile+".csv")

if __name__ == "__main__":
    analyse_structure(args.file)
