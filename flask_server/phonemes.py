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

    sites_df = pd.DataFrame.from_dict(input_json['sites'])
    idcolumn = 'Site' if 'Site' in sites_df else 'id'
    # TODO: Iterate over experiments, perform PHONEMeS for each separately - create an input file with two experiments
    experiment_columns = [col for col in sites_df.columns if col != idcolumn]

    phonemes_pkn_ksn = pd.read_csv(PHONEMES_PKN_KSN)
    all_pkn_ksn_nodes = set(np.concatenate(
        [phonemes_pkn_ksn['source'].values, phonemes_pkn_ksn['target'].values]))

    sites_df = sites_df[sites_df['Site'].apply(lambda site: site in all_pkn_ksn_nodes)]

    target_series = pd.concat([
        pd.Series(data=1, index=input_json['targets']['up']),
        pd.Series(data=-1, index=input_json['targets']['down'])
    ])

    # TODO: Here you should iterate over experiments and create a file for each one
    # Don't forget to drop NAs if a site appears in one experiment but not in another
    Path.mkdir(output_dir, parents=True, exist_ok=True)
    output_prefix = output_dir / f'{filepath.stem}'
    sites_df.to_csv(str(output_prefix) + '_sites.csv', index=False)
    target_series.to_csv(str(output_prefix) + '_targets.csv', index=True, header=False)

    return output_prefix


def run_phonemes(file_prefix: Path) -> Path:
    # TODO: Iterate experiments again
    output_dir = file_prefix.parent
    output_path = output_dir / f'phonemes_out.sif'
    subprocess.run(["Rscript",
                    "run_phonemes.R",
                    str(file_prefix) + '_sites.csv',
                    str(file_prefix) + '_targets.csv',
                    output_path
                    ])
    return output_path


def run_cytoscape(filepath: Path) -> Path:
    phonemes_df = pd.read_csv(filepath)

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

    network = phonemes_df.rename(
        columns={"Node1": "source", "Node2": "target", "Weight": "weight", "Sign": "sign"})

    nodes_values = network[["source", "target"]].values.ravel()
    nodes_data = pd.unique(nodes_values)
    nodes = pd.DataFrame(data=nodes_data, columns=["id"])
    p4c.create_network_from_data_frames(nodes=nodes, edges=network, collection="phonemes2cytoscape",
                                        # TODO: Put experiment name
                                        title='Experiment01')
    p4c.layout_network()
    output_dir = filepath.parent
    output_path = output_dir / f'cytoscape_out.cx'
    p4c.export_network(filename=str(output_path), type='CX', overwrite_file=True)
    return output_path


def create_pathway_skeleton(filepath: Path) -> Path:
    # TODO: One final time: Iterate experiments
    with open(filepath) as infile:
        cx_json = json.load(infile)

    # For some reason this is in a singleton-list format, so we extract all of those nested keys into one single json object
    cx_unnested = {key: val for singleton in cx_json for key, val in singleton.items()}
    nodes_dict = {entry['@id']: entry['n'] for entry in cx_unnested['nodes']}
    skeleton_json = {'pathway': {'name': 'Experiment01'},
                     'nodes': [{
                         'id': node['node'],
                         'geneNames': [nodes_dict[node['node']]],
                         'type': 'gene_protein',
                         'x': node['x'],
                         # Trial and error showed we need to downscale the y-axis a bit
                         'y': node['y']*0.75} for node in cx_unnested['cartesianLayout']],
                     'links': [{
                         'id': f"relation-{relation['@id']}",
                         'sourceId': relation['s'],
                         'targetId': relation['t'],
                         'types': ['other']
                     }
                         for relation in cx_unnested['edges']]}
    output_dir = filepath.parent
    output_path = output_dir / f'json_skeleton.json'
    with open(output_path, 'w') as outfile:
        json.dump(skeleton_json, outfile)

    return output_path
