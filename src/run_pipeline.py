import argparse
import logging
import os

from src.processing_functions.logging import setup_logger
from processing_functions.generic_process import create_dirs
from processing_functions.yaml_config import process_yaml_file

# TODO: MIRAR LAS COLUMNAS OPCIONALES Y OBlIGATORIAS DE CADA ELEMENTO


def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description='Python script to run combined NGS pipelines.')
    
    # Add argument for log level (required)
    parser.add_argument('--project', help='Name of the project. It has to be the same as the \
                        folder name located in the projects folder.')
    parser.add_argument('--verbose', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='WARNING',
                        help='Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)')
    parser.add_argument('--yaml', default='config.yaml',
                        help='Configuration yaml file.')
    parser.add_argument('--samplesheet', default='samplesheet.csv',
                        help='File with information of fastq files and other parameters for use.')
    


    # Parse command-line arguments
    args = parser.parse_args()

    log_file = f"results/{args.project}/{args.project}.log"


    # Set up logger based on command-line arguments
    setup_logger(getattr(logging, args.verbose), log_file)



    # Create work and results dir for the project
    create_dirs(args=args)
    
    

    # Parse config.yaml
    # TODO: COMPROBAR OUTPUTS
    list_samplesheets, list_dbs_to_download, list_pipeline_commands = process_yaml_file(yaml_file=f"projects/{args.project}/{args.yaml}")


    create_samplesheets(list_samplesheets)


    download_DBs(list_dbs_to_download)


    run_pipeline(list_pipeline_commands)





if __name__ == '__main__':
    main()