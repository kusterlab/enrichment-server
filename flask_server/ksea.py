from pathlib import Path
import json
import pandas as pd
import numpy as np
import kinact


def perform_ksea(filepath: Path, session_id: str, dataset_name: str) -> Path:
    input_json = json.load(open(filepath))
    input_df = pd.DataFrame.from_dict(input_json)
    input_df.set_index('Site', inplace=True)

    adjacency_matrix = pd.read_csv('../db/psp_kinase_substrate_adjacency_matrix.csv', skiprows=1).set_index('p_site')

    ksea_results = []
    for experiment in input_df:
        scores, p_values = kinact.ksea.ksea_mean(
            data_fc=input_df[experiment],
            interactions=adjacency_matrix,
            mP=input_df.values.mean(),
            delta=input_df.values.std())
        res = pd.DataFrame({f'Score ({experiment})': scores, f'-log(p) ({experiment})': -np.log10(p_values)})
        ksea_results.append(res)
    ksea_results_df = pd.concat(ksea_results, axis=1)
    ksea_results_df.index.name = 'Gene'
    output_dir = Path('..') / session_id / dataset_name
    Path.mkdir(output_dir, parents=True, exist_ok=True)
    output_json = output_dir / f'{filepath.stem}_ksea.json'
    ksea_results_df.reset_index().to_json(path_or_buf=output_json,
                                          orient='records',
                                          indent=1  # For DEBUG
                                          )
    return output_json
