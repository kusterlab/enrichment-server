import subprocess
import json
from pathlib import Path
import requests
import time
import numpy as np
import pandas as pd
import py4cytoscape as p4c

PHONEMES_PKN = Path('../db/phonemesPKN.csv')
PHONEMES_KSN = Path('../db/phonemesKSN.csv')
PHONEMES_PKN_KSN = Path('../db/phonemes_PKN_KSN.csv')
CYTOSCAPE_PATH = '../Cytoscape_v3.10.1/Cytoscape'


def preprocess_phonemes(filepath: Path) -> Path:
    output_dir = filepath.parent
    input_json = json.load(open(filepath))
    Path.mkdir(output_dir, parents=True, exist_ok=True)
    output_prefix = output_dir / f'{filepath.stem}'

    sites_df = pd.DataFrame.from_dict(input_json['sites'])
    idcolumn = 'Site' if 'Site' in sites_df else 'id'
    experiment_columns = [col for col in sites_df.columns if col != idcolumn]

    phonemes_pkn_ksn = pd.read_csv(PHONEMES_PKN_KSN)
    all_pkn_ksn_nodes = set(np.concatenate(
        [phonemes_pkn_ksn['source'].values, phonemes_pkn_ksn['target'].values]))

    sites_df = sites_df[sites_df['Site'].apply(lambda site: site in all_pkn_ksn_nodes)]

    # Preprocess so that each experiment gets an own 'sites' and 'targets' file.
    # An additional output file contains the list of experiments
    for experiment in experiment_columns:
        target_series = pd.concat([
            pd.Series(data=1, index=input_json['targets'][experiment]['up']),
            pd.Series(data=-1, index=input_json['targets'][experiment]['down'])
        ])

        sites_df[[idcolumn, experiment]].dropna().to_csv(f'{output_prefix}_{experiment}_sites.csv', index=False)
        target_series.to_csv(f'{output_prefix}_{experiment}_targets.csv', index=True, header=False)

    # Create a file containing the experiment names
    with open(str(output_prefix) + '_experiments.csv', 'w') as experiment_outfile:
        experiment_outfile.write(','.join(experiment_columns))

    return output_prefix


def run_phonemes(file_prefix: Path) -> Path:
    output_dir = file_prefix.parent

    experiment_file = output_dir / 'input_experiments.csv'
    with open(experiment_file) as infile:
        experiments = infile.read().split(',')

    for experiment in experiments:
        output_path = output_dir / f'{experiment}_phonemes_out.sif'
        subprocess.run(["Rscript",
                        "run_phonemes.R",
                        f'{file_prefix}_{experiment}_sites.csv',
                        f'{file_prefix}_{experiment}_targets.csv',
                        output_path
                        ])

    return output_dir


def run_cytoscape(phonemes_outputfolder: Path) -> Path:
    cytoscape = subprocess.Popen(CYTOSCAPE_PATH)
    # Wait until cytoscape is ready
    not_found = True
    while not_found:
        try:
            p4c.cytoscape_ping()
            not_found = False
        except (
                requests.exceptions.RequestException,
                requests.exceptions.ConnectionError,
                requests.exceptions.HTTPError):
            time.sleep(.5)

    experiment_file = phonemes_outputfolder / 'input_experiments.csv'
    with open(experiment_file) as infile:
        experiments = infile.read().split(',')

    for experiment in experiments:
        filepath = phonemes_outputfolder / f'{experiment}_phonemes_out.sif'
        phonemes_df = pd.read_csv(filepath)
        network = phonemes_df.rename(
            columns={"Node1": "source", "Node2": "target", "Weight": "weight", "Sign": "sign"})

        nodes_values = network[["source", "target"]].values.ravel()
        nodes_data = pd.unique(nodes_values)
        nodes = pd.DataFrame(data=nodes_data, columns=["id"])
        p4c.create_network_from_data_frames(nodes=nodes, edges=network, collection="phonemes2cytoscape",
                                            # TODO: Put experiment name
                                            title=experiment)
        p4c.layout_network()
        output_path = phonemes_outputfolder / f'{experiment}_cytoscape_out.cx'
        p4c.export_network(filename=str(output_path), type='CX', overwrite_file=True)

    cytoscape.kill()
    return phonemes_outputfolder


def create_pathway_skeleton(cytoscape_outputfolder: Path) -> Path:
    experiment_file = cytoscape_outputfolder / 'input_experiments.csv'
    with open(experiment_file) as infile:
        experiments = infile.read().split(',')

    pathway_skeleton_list = []
    for experiment in experiments:
        filepath = cytoscape_outputfolder / f'{experiment}_cytoscape_out.cx'
        with open(filepath) as infile:
            cx_json = json.load(infile)

        # For some reason this is in a singleton-list format, so we extract all of those nested keys into one single json object
        cx_unnested = {key: val for singleton in cx_json for key, val in singleton.items()}
        nodes_dict = {entry['@id']: entry['n'] for entry in cx_unnested['nodes']}
        skeleton_json = {'pathway': {'name': experiment},
                         'nodes': [{
                             'id': node['node'],
                             'geneNames': [nodes_dict[node['node']]],
                             'type': 'gene_protein',
                             'x': node['x'],
                             # Trial and error showed we need to downscale the y-axis a bit
                             'y': node['y'] * 0.75} for node in cx_unnested['cartesianLayout']],
                         'links': [{
                             'id': f"relation-{relation['@id']}",
                             'sourceId': relation['s'],
                             'targetId': relation['t'],
                             'types': ['other']
                         }
                             for relation in cx_unnested['edges']]}
        #Make sure all 'y's are positive by getting the minimum and subtracting it from all
        miny = min([node['y'] for node in skeleton_json['nodes']])
        for node in skeleton_json['nodes']:
            node['y'] -= miny - 100 #The additional 100 makes sure the minimum y is 100, so should be well within the visible range.

        pathway_skeleton_list.append(skeleton_json)

    output_path = cytoscape_outputfolder / f'json_skeletons.json'
    with open(output_path, 'w') as outfile:
        json.dump(pathway_skeleton_list, outfile)

    return output_path
