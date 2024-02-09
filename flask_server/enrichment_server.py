# How to send a request: curl -X POST -F file=@minimal_input.gct -F session_id=ABCDEF12345  -F dataset_name=FooBar http://127.0.0.1:5000/
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


@app.route('/enrichment/ptm-sea', methods=['POST'])
def handle_ptmsea_request() -> werkzeug.wrappers.Response | str:
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

    #Preprocess the json input into a gct file
    ptmsea_input = preprocess_ptmsea(filepath, form['session_id'], form['dataset_name'])

    ptmsea_combined_output = run_ptmsea(ptmsea_input)

    # TODO: Delete everything when done to not accumulate files?
    return send_file(postprocess_ptmsea(ptmsea_combined_output), as_attachment=True)


def preprocess_ptmsea(filepath: Path, session_id: str, dataset_name: str) -> Path:
    output_dir = Path('..') / session_id / dataset_name
    Path.mkdir(output_dir, parents=True, exist_ok=True)
    input_json = json.load(open(filepath))
    input_df = pd.DataFrame.from_dict(input_json).rename({'flanking': 'id'}, axis=1)
    # There is a method cmapPy.pandasGEXpress.write_gct,
    # but I could not get it to run
    # (it tries to use DataFrame.at[] with ranges and I couldn't find a pandas version where this works)
    # So I wrote a minimal version myself
    input_gct_file = output_dir / 'ptmsea_input.gct'
    with open(input_gct_file, 'w') as f:
        # Write Version and dims
        f.write("#1.3\n")
        f.write(f"{input_df.shape[0]}\t{input_df.shape[1] - 1}\t0\t0\n")
        # Write data df
        input_df.to_csv(f, sep='\t', index=False)
    return input_gct_file


def run_ptmsea(filepath: Path) -> Path:
    output_dir = filepath.parent
    Path.mkdir(output_dir, parents=True, exist_ok=True)
    output_prefix = output_dir / 'ptmsea_out'
    subprocess.run(["Rscript",
                    "../ssGSEA2.0/ssgsea-cli.R",
                    "-i", str(Path('..') / 'flask_server' / filepath),
                    "-o", str(output_prefix),
                    "-d", "../ssGSEA2.0/db/ptmsigdb/ptm.sig.db.all.flanking.human.v2.0.0.gmt",
                    "-w", "0.75",
                    "-e", "FALSE",
                    ])
    return Path(str(output_prefix) + '-combined.gct')


def postprocess_ptmsea(output_gct) -> Path:
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
                          # indent=1  # For DEBUG
                          )
    return output_json


if __name__ == '__main__':
    app.run(debug=True)
