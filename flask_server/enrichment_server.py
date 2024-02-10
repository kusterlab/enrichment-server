# How to send a request: curl -X POST -F file=@<input_file> -F session_id=ABCDEF12345  -F dataset_name=FooBar http://127.0.0.1:5000/<route>
import subprocess
from pathlib import Path
import json

import werkzeug.wrappers
from werkzeug.utils import secure_filename
from flask import Flask, request, send_file, jsonify
import flask.wrappers
import pandas as pd
from cmapPy.pandasGEXpress import parse_gct

app = Flask(__name__)


@app.route('/', methods=['GET'])
def get_status() -> flask.wrappers.Response:
    return jsonify(status=200)


@app.route('/ssgsea/<string:ssgsea_type>', methods=['POST'])
def handle_ssgsea_request(ssgsea_type) -> werkzeug.wrappers.Response | str:
    valid_ssgsea_types = ['ssc', 'gc', 'gcr']

    if ssgsea_type not in valid_ssgsea_types:
        return f"Invalid 'ssgsea_type'. Allowed values are {', '.join(valid_ssgsea_types)}"

    # Check if the POST request has the file part
    if 'file' not in request.files:
        return 'Error: No file part in the request.\n'

    file = request.files['file']
    form = request.form
    required_parameters = ['session_id', 'dataset_name']
    for param in required_parameters:
        if param not in form:
            return f'Error: parameter {param} not specified.\n'

    # If the user does not select a file, the browser might submit an empty part without a filename
    if file.filename == '':
        return 'Error: No input file given.\n'

    filepath = Path(secure_filename(file.filename))
    file.save(filepath)

    # Preprocess the json input into a gct file
    ssgsea_input = preprocess_ssgsea(filepath, form['session_id'], form['dataset_name'], ssgsea_type == 'gc')

    ssgsea_combined_output = run_ssgsea(ssgsea_input, ssgsea_type)

    # TODO: Delete everything when done to not accumulate files?
    return send_file(postprocess_ssgsea(ssgsea_combined_output), as_attachment=True)


def preprocess_ssgsea(filepath: Path, session_id: str, dataset_name: str, type_is_gc) -> Path:
    output_dir = Path('..') / session_id / dataset_name
    Path.mkdir(output_dir, parents=True, exist_ok=True)
    input_json = json.load(open(filepath))
    input_df = pd.DataFrame.from_dict(input_json)

    # If it's a non-redundant gene-centric ssGSEA, we need to eliminate duplicates
    if type_is_gc:
        def abs_max_signed(group):
            idx = group.abs().idxmax()
            return group.loc[idx] if pd.notna(idx) else float('nan')

        experiment_columns = [col for col in input_df.columns if col != 'id']
        input_df_grouped = input_df.groupby('id')
        unique_values = [input_df_grouped[exp].apply(abs_max_signed) for exp in experiment_columns]
        input_df = pd.DataFrame(unique_values).T.reset_index()

    # There is a method cmapPy.pandasGEXpress.write_gct,
    # but I could not get it to run
    # (it tries to use DataFrame.at[] with ranges and I couldn't find a pandas version where this works)
    # So I wrote a minimal version myself
    input_gct_file = output_dir / 'ssgsea_input.gct'
    with open(input_gct_file, 'w') as f:
        # Write Version and dims
        f.write("#1.3\n")
        f.write(f"{input_df.shape[0]}\t{input_df.shape[1] - 1}\t0\t0\n")
        # Write data df
        input_df.to_csv(f, sep='\t', index=False)
    return input_gct_file


def run_ssgsea(filepath: Path, ssgsea_type) -> Path:
    output_dir = filepath.parent
    Path.mkdir(output_dir, parents=True, exist_ok=True)
    output_prefix = output_dir / f'ssgsea_{ssgsea_type}_out'

    match ssgsea_type:
        case 'ssc':
            database = "../ssGSEA2.0/db/ptmsigdb/ptm.sig.db.all.flanking.human.v2.0.0.gmt"
        case 'gc' | 'gcr':
            database = "../db/c2.cp.kegg+wp.v2023.2.Hs.symbols.gmt"
        case _:
            # This is unreachable, the parent function would have errored out already
            database = None

    subprocess.run(["Rscript",
                    "../ssGSEA2.0/ssgsea-cli.R",
                    "-i", str(Path('..') / 'flask_server' / filepath),
                    "-o", str(output_prefix),
                    "-d", database,
                    "-w", "0.75",
                    "-e", "FALSE",
                    ])
    return Path(str(output_prefix) + '-combined.gct')


def postprocess_ssgsea(output_gct) -> Path:
    gct_parsed = parse_gct.parse(output_gct)
    experiment_names = gct_parsed.data_df.columns
    gct_df_joined = gct_parsed.row_metadata_df[[f'Signature.set.overlap.percent.{exp}' for exp in experiment_names]
                                               + [f'fdr.pvalue.{exp}' for exp in experiment_names]
                                               ].join(gct_parsed.data_df).reset_index()

    gct_df_joined.columns = (['Signature ID']
                             + [f'Overlap [%] ({exp})' for exp in experiment_names]
                             + [f'Adj. p-value ({exp})' for exp in experiment_names]
                             + [f'Score ({exp})' for exp in experiment_names])

    output_json = output_gct.parent / f'{output_gct.stem}.json'
    gct_df_joined.to_json(path_or_buf=output_json, orient='records',
                          indent=1  # For DEBUG
                          )
    return output_json


if __name__ == '__main__':
    app.run(debug=True)
