from pathlib import Path
import json
import subprocess
import pandas as pd
import kinact


def preprocess_ksea(filepath: Path) -> Path:
    output_dir = filepath.parent
    input_json = json.load(open(filepath))
    input_df = pd.DataFrame.from_dict(input_json)

    # Todo: Merge with preprocess_ssgsea function
    idcolumn = 'Site' if 'Site' in input_df else 'id'

    def abs_max_signed(group):
        idx = group.abs().idxmax()
        return group.loc[idx] if pd.notna(idx) else float('nan')

    experiment_columns = [col for col in input_df.columns if col != idcolumn]
    input_df_grouped = input_df.groupby(idcolumn)
    unique_values = [input_df_grouped[exp].apply(abs_max_signed) for exp in experiment_columns]
    input_df = pd.DataFrame(unique_values).T.reset_index()

    Path.mkdir(output_dir, parents=True, exist_ok=True)

    output_csv = output_dir / f'{filepath.stem}.csv'
    input_df.to_csv(output_csv, index=False)
    return output_csv


def run_rokai(filepath: Path) -> Path:
    output_path = filepath.parent / f'rokai_result.csv'
    subprocess.run(["Rscript",
                    "modules/ksea/run_rokai.R",
                    str(filepath),
                    str(output_path)])
    return output_path


def perform_ksea(filepath: Path) -> Path:
    input_df = pd.read_csv(filepath)
    input_df.set_index('Site', inplace=True)

    adjacency_matrix = pd.read_csv('../db/psp_kinase_substrate_adjacency_matrix.csv', skiprows=1).set_index('p_site')

    ksea_results = []
    for experiment in input_df:
        try:
            scores, p_values = kinact.ksea.ksea_mean(
                data_fc=input_df[experiment],
                interactions=adjacency_matrix,
                mP=input_df[experiment].mean(),
                delta=input_df[experiment].std())
            res = pd.DataFrame({f'Score ({experiment})': scores, f'adj p-val ({experiment})': p_values})
            ksea_results.append(res)
        except ZeroDivisionError:
            continue
    output_json = filepath.parent / f'ksea_result.json'

    if len(ksea_results) == 0:
        with open(output_json, 'w') as o:
            o.write(str(ksea_results))
    else:
        ksea_results_df = pd.concat(ksea_results, axis=1)
        ksea_results_df.index.name = 'Gene'

        ksea_results_df.reset_index().to_json(path_or_buf=output_json,
                                              orient='records',
                                              # indent=1  # For DEBUG
                                              )
    return output_json
