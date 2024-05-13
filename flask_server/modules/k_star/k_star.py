from pathlib import Path
import json
import pickle

import numpy as np
import pandas as pd
import psite_annotation as pa
from kstar import helpers, calculate, mapping, config


def run_kstar(filepath: Path) -> Path:
    output_dir = filepath.parent
    input_json = json.load(open(filepath))
    input_df = pd.DataFrame.from_dict(input_json)
    data_columns = [col for col in input_df if
                    col not in ['Modified sequence', 'Proteins']]

    # We need to convert the sequences into +/-7 flanking format with modified residues in lowercase
    input_df = pa.addPeptideAndPsitePositions(input_df, '../db/Phosphosite_seq.fasta', pspInput=True,
                                              context_left=7, context_right=7, retain_other_mods=True)

    input_df['Uniprot_Accession'] = input_df['Matched proteins'].apply(lambda prot: prot.split(';')[0])

    input_df['Sequence'] = input_df['Site sequence context'].apply(lambda prot: prot.split(';')[0])

    input_df = input_df[['Uniprot_Accession', 'Sequence'] + data_columns]

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

    # Test if there is enough evidence to perform ST and/or Y enrichment, only then perform it
    result = dict(ST=[], Y=[])
    for phospho_type in ['ST', 'Y']:
        for direction in ['up', 'down']:
            kinact = calculate.KinaseActivity(exp_mapper.experiment,
                                              activity_log,
                                              phospho_type=phospho_type)
            threshold_test = kinact.test_threshold(agg='mean',
                                                   threshold=0,
                                                   greater=(direction == 'up'),
                                                   return_evidence_sizes=True)

            if threshold_test.min() > 0:
                kinact_dict = calculate.enrichment_analysis(exp_mapper.experiment, activity_log, networks,
                                                            phospho_types=[phospho_type],
                                                            # We already filtered for regulations, so 0 is an acceptable threshold
                                                            agg='mean', threshold=0,
                                                            # We expect kinase inhibition, so check for values smaller than the threshold
                                                            greater=(direction == 'up'), PROCESSES=4)

                result_df = np.log10(kinact_dict[phospho_type].activities) * (-1 if direction == 'up' else 1)

                # Post Process and convert into JSON
                result_list = result_df.rename({
                    # Trim away the 'data:'
                    col: col[5:] for col in kinact_dict[phospho_type].activities.columns
                }, axis=1).reset_index(names='Kinase').to_dict(orient='records')

                result[phospho_type] += result_list

    output_json = output_dir / f'kstar_result.json'
    with open(output_json, 'w') as outfile:
        json.dump(result, outfile)

    return output_json
