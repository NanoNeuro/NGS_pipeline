import logging

def setup_logger(log_level, log_file=None):
    # Set up logging configuration
    logging.basicConfig(level=log_level,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename=log_file)
