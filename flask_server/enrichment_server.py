# How to send a request: curl -X POST -F file=@minimal_input.gct -F session_id=ABCDEF12345  -F dataset_name=FooBar http://127.0.0.1:5000/
from werkzeug.utils import secure_filename
import subprocess
from pathlib import Path
from flask import Flask, request, send_file, jsonify

app = Flask(__name__)


@app.route('/', methods=['GET'])
def get_status():
    return jsonify(status=200)


@app.route('/enrichment/ptm-sea', methods=['POST'])
def handle_ptmsea_request():
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

    # TODO: Post-process result and return as JSON
    # TODO: Delete everything when done to not accumulate files?
    return send_file(ptmsea_combined_output, as_attachment=True)


def run_ptmsea(filepath: Path, session_id: str, dataset_name: str):
    output_dir = Path('..') / session_id / dataset_name
    Path.mkdir(output_dir, parents=True, exist_ok=True)
    output_prefix = output_dir / 'ptmsea_out'
    R = subprocess.run(["Rscript",
                        "../ssGSEA2.0/ssgsea-cli.R",
                        "-i", str(Path('..') / 'flask_server' / filepath),
                        "-o", str(output_prefix),
                        "-d", "../ssGSEA2.0/db/ptmsigdb/ptm.sig.db.all.flanking.human.v2.0.0.gmt",
                        "-w", "0.75",
                        "-e", "FALSE",
                        ])
    return str(output_prefix) + '-combined.gct'


if __name__ == '__main__':
    app.run(debug=True)
