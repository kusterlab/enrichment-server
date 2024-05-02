from pathlib import Path
import json
import pickle

import pandas as pd
from kstar import helpers, calculate, mapping, config


def run_kstar(filepath: Path) -> Path:
    output_dir = filepath.parent
    input_json = json.load(open(filepath))
    input_df = pd.DataFrame.from_dict(input_json)
    data_columns = [col for col in input_df if
                    col not in ['Sequence', 'Uniprot_Accession'] and not col.startswith('KSTAR')]
    # KSTAR automatically recognizes data columns if you prepend them with 'data:'
    input_df = input_df.rename({col: f'data:{col}' for col in data_columns}, axis=1)
    # Map the data to KinPred
    Path.mkdir(output_dir / 'MAPPED_DATA', parents=True, exist_ok=True)
    mapping_log = helpers.get_logger('mapping_log', output_dir / 'MAPPED_DATA' / 'mapping_log.log')
    mapDict = {'peptide': 'Sequence', 'accession_id': 'Uniprot_Accession'}
    exp_mapper = mapping.ExperimentMapper(experiment=input_df,
                                          columns=mapDict,
                                          logger=mapping_log)
    # Now perform the activity scoring
    Path.mkdir(output_dir / 'RESULTS', parents=True, exist_ok=True)
    activity_log = helpers.get_logger('activity_log', output_dir / 'RESULTS' / 'activity_log.log')
    # Load the pickles containing the 50 pruned networks for s/t kinases
    networks = {}
    networks['ST'] = pickle.load(open(config.NETWORK_ST_PICKLE, "rb"))
    networks['Y'] = pickle.load(open(config.NETWORK_Y_PICKLE, "rb"))
    kinact_dict = calculate.enrichment_analysis(exp_mapper.experiment, activity_log, networks,
                                                phospho_types=['ST', 'Y'],
                                                # We already filtered for regulations, so 0 is an acceptable threshold
                                                agg='mean', threshold=0,
                                                # We expect kinase inhibition, so check for values smaller than the threshold
                                                greater=False, PROCESSES=4)
    # Post Process and convert into JSON
    st_result_dict = kinact_dict['ST'].activities.rename({
        # Trim away the 'data:'
        col: col[5:] for col in kinact_dict['ST'].activities.columns
    }, axis=1).reset_index(names='Kinase').to_dict(orient='records')

    y_result_dict = kinact_dict['Y'].activities.rename({
        # Trim away the 'data:'
        col: col[5:] for col in kinact_dict['Y'].activities.columns
    }, axis=1).reset_index(names='Kinase').to_dict(orient='records')

    result = {'ST': st_result_dict, 'Y': y_result_dict}

    output_json = output_dir / f'kstar_result.json'
    with open(output_json, 'w') as outfile:
        json.dump(result, outfile)

    return output_json
