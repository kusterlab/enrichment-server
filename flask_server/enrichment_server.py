# How to send a request: curl -X POST -F file=@minimal_input.gct -F session_id=ABCDEF12345  -F dataset_name=FooBar http://127.0.0.1:5000/
import flask.wrappers
import werkzeug.wrappers
from werkzeug.utils import secure_filename
import subprocess
from pathlib import Path
from flask import Flask, request, send_file, jsonify
import pandas as pd
import numpy as np
from cmapPy.pandasGEXpress import parse_gct, concat

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

    ptmsea_combined_output = run_ptmsea(filepath, form['session_id'], form['dataset_name'])

    # TODO: Delete everything when done to not accumulate files?
    return send_file(postprocess_ptmsea(ptmsea_combined_output), as_attachment=True)


def run_ptmsea(filepath: Path, session_id: str, dataset_name: str) -> Path:
    output_dir = Path('..') / session_id / dataset_name
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
