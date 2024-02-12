# How to send a request:
# curl -X POST -F file=@<input_file> -F session_id=ABCDEF12345  -F dataset_name=FooBar http://127.0.0.1:5000/<route>
from pathlib import Path

import werkzeug.wrappers
from werkzeug.utils import secure_filename
from flask import Flask, request, send_file, jsonify
import flask.wrappers
import ssgsea
import ksea

app = Flask(__name__)


@app.route('/', methods=['GET'])
def get_status() -> flask.wrappers.Response:
    return jsonify(status=200)


@app.route('/ssgsea/<string:ssgsea_type>', methods=['POST'])
def handle_ssgsea_request(ssgsea_type) -> werkzeug.wrappers.Response | str:
    valid_ssgsea_types = ['ssc', 'gc', 'gcr']

    if ssgsea_type not in valid_ssgsea_types:
        return f"Invalid 'ssgsea_type'. Allowed values are {', '.join(valid_ssgsea_types)}"

    post_request_processed = process_post_request(request)

    if len(post_request_processed) == 1:
        return post_request_processed

    filepath, session_id, dataset_name = post_request_processed

    # Preprocess the json input into a gct file
    ssgsea_input = ssgsea.preprocess_ssgsea(filepath, session_id, dataset_name, ssgsea_type == 'gc')

    ssgsea_combined_output = ssgsea.run_ssgsea(ssgsea_input, ssgsea_type)

    # TODO: Delete everything when done to not accumulate files?
    return send_file(ssgsea.postprocess_ssgsea(ssgsea_combined_output), as_attachment=True)


@app.route('/ksea', methods=['POST'])
@app.route('/ksea/<string:ksea_type>', methods=['POST'])
def handle_ksea_request(ksea_type=None) -> werkzeug.wrappers.Response | str:
    post_request_processed = process_post_request(request)

    if len(post_request_processed) == 1:
        return post_request_processed

    filepath, session_id, dataset_name = post_request_processed

    preprocessed_filepath = ksea.preprocess_ksea(filepath, session_id, dataset_name)
    if ksea_type == 'rokai':
        preprocessed_filepath = ksea.run_rokai(preprocessed_filepath)

    ksea_result = ksea.perform_ksea(preprocessed_filepath)
    return send_file(ksea_result, as_attachment=True)


def process_post_request(post_request: werkzeug.Request) -> tuple[Path, str, str] | str:
    # Check if the POST request has the file part
    if 'file' not in post_request.files:
        return 'Error: No file part in the request.\n'

    file = post_request.files['file']
    form = post_request.form
    required_parameters = ['session_id', 'dataset_name']
    for param in required_parameters:
        if param not in form:
            return f'Error: parameter {param} not specified.\n'

    # If the user does not select a file, the browser might submit an empty part without a filename
    if file.filename == '':
        return 'Error: No input file given.\n'

    filepath = Path(secure_filename(file.filename))
    file.save(filepath)
    return filepath, form['session_id'], form['dataset_name']


if __name__ == '__main__':
    app.run(debug=True)
