import pandas as pd
import numpy as np

from .yaml_config import NF_CONFIG_PIPELINES


DICT_ACCEPTED_COLUMNS = {'rnaseq': ['sample', 'fastq_1', 'fastq_2', 'strandedness'],
                         'scrnaseq': ['sample', 'fastq_1', 'fastq_2', 'protocol', 'expected_cells', 'seq_center'],
                         'smrnaseq': ['sample', 'fastq_1', 'fastq_2'], 
                         'circrna': ['sample', 'fastq_1', 'fastq_2'],
                         'circdna': ['sample', 'fastq_1', 'fastq_2'],
                         'taxprofiler': ['sample', 'fastq_1', 'fastq_2', 'strandedness']}




# TODO CHECKEAR QUE EL INPUT Y EL OUTPUT DE LOS FASTQS INTERMEDIARIOS ES EL CORRECTO
def create_samplesheets(master_samplesheet_path, project, yaml_dict, list_samplesheets):
    """
    This function is used to create the samplesheets required fo each process of any project.
    The way this function works is by loading the master samplesheet and, for any process, create a subset of the
    samplesheet.
        - If the process name is not located in the master samplesheet, it throws an error.
        - If the columns required for the process type are not present, it throws an error.

    Then it creates the process-specific samplesheets in the work directory. 
    
    In the case of taxprofiler, where the location of the 1st and 2nd mapping is performed and dependent on the user, 
    this function will create the spreadsheets accordingly. This is because the sample names in the original fastqs
    may not be the same as in the processed ones. 
    """



    master_samplesheet = pd.read_csv(master_samplesheet_path, sep=',', dtype=str)

    for process_name, pipeline, path_to_process_csv in list_samplesheets:
        samplesheet_process = master_samplesheet.loc[[process_name in i.split(';') for i in master_samplesheet['process'].values]]
        samplesheet_process = samplesheet_process.loc[:, [i for i in DICT_ACCEPTED_COLUMNS[pipeline] if i in samplesheet_process.columns]]
    
        if pipeline in NF_CONFIG_PIPELINES:
            pass

        elif pipeline == 'taxprofiler':
            list_samples = samplesheet_process['sample'].drop_duplicates().values

            path_out_1 = f"results/<PROJECT>/{pipeline}/1st_map/star_salmon/unmapped/<SAMPLE>.unmapped_<STRAND>.fastq.gz"
            samplesheet_process_map_1 = create_samplesheet_map(list_samples, samplesheet_process, path_out_1, project)

            path_out_2 = f"results/<PROJECT>/{pipeline}/2nd_map/<SAMPLE>.unmapped.fastq.<STRAND>.gz"
            samplesheet_process_map_2 = create_samplesheet_map(list_samples, samplesheet_process, path_out_2, project)
            

            samplesheet_process_map_1.to_csv(path_to_process_csv.replace('.csv', '_map_1.csv'), index=None)
            samplesheet_process_map_2.to_csv(path_to_process_csv.replace('.csv', '_map_2.csv'), index=None)


        samplesheet_process.to_csv(path_to_process_csv, index=None)







def create_samplesheet_map(list_samples, samplesheet_process, path_out, project):
    samplesheet_process_map = pd.DataFrame({
                    'sample':list_samples, 
                    'fastq_1': [path_out.replace('<PROJECT>', project).replace('<SAMPLE>', sample).replace('<STRAND>', '1') for sample in list_samples],
                    'fastq_2': [path_out.replace('<PROJECT>', project).replace('<SAMPLE>', sample).replace('<STRAND>', '2') for sample in list_samples]
                    })

    # If there are samples that have nan fastqs, we need to set these values as false.
    if 'fastq_2' not in samplesheet_process.columns:
        samplesheet_process_map['fastq_2'] = np.NaN
    else:
        samplesheet_process_drop = samplesheet_process.drop_duplicates(subset=['sample'])
        isnan = pd.isnull(samplesheet_process_drop['fastq_2']).values
        samplesheet_process_map.loc[isnan, 'fastq_2'] = np.NaN

    return samplesheet_process_map