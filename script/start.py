import os
import configparser

from scripts_package.genotype_formatting import genotype_binarization, genotype_from_binarised_txt
from scripts_package.phenotype_formatting import phenotype_formating, convert_grin_to_csv

def run_processes(config_path):
    config = configparser.ConfigParser()
    config.read(config_path)

    for section in config.sections():
        directory_path = config[section]['path']  # Ensure path is correctly formatted without extra quotes
        try:
            condition = config[section]['condition']
        except:
            condition = None
        
        os.chdir(directory_path)  # Change working directory
        try:
            config[section]['phenotype_categories']
            genotype_from_binarised_txt.main()
        except:
            # Running the genotype binarization script
            genotype_binarization.main(condition)
            
        # Convert GRIN to CSV
        convert_grin_to_csv.main(config[section])

        # # Phenotype Formatting
        phenotype_formating.main(config[section])

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run processes based on configuration.")
    parser.add_argument("config", help="Path to the configuration INI file.")
    args = parser.parse_args()


    run_processes(args.config)
