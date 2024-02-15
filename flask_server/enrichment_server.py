# How to send a request:
# curl -X POST -F file=@<input_file> -F session_id=ABCDEF12345  -F dataset_name=FooBar http://127.0.0.1:1234/<route>
from pathlib import Path
from urllib.parse import urlparse
import werkzeug.wrappers
from werkzeug.utils import secure_filename
from flask import Flask, request, send_file, jsonify, make_response
import flask.wrappers
import ssgsea
import ksea

app = Flask(__name__)


@app.route('/', methods=['GET'])
def get_status() -> flask.wrappers.Response:
    return send_response(jsonify(status=200))


# TODO: In the second route, the ssgsea_type actually can only be ssc. Can I enforce this?
@app.route('/ssgsea/<string:ssgsea_type>', methods=['POST'])
@app.route('/ssgsea/<string:ssgsea_type>/<string:ssc_input_type>', methods=['POST'])
def handle_ssgsea_request(ssgsea_type, ssc_input_type='flanking') -> werkzeug.wrappers.Response | str:
    valid_ssgsea_types = ['ssc', 'gc', 'gcr']

    if ssgsea_type not in valid_ssgsea_types:
        return f"Invalid 'ssgsea_type'. Allowed values are {', '.join(valid_ssgsea_types)}"

    valid_ssc_input_types = ['flanking', 'uniprot']
    if ssc_input_type not in valid_ssc_input_types:
        return f"Invalid 'ssc_input_type'. Allowed values are {', '.join(valid_ssc_input_types)}"

    post_request_processed = process_post_request(request)

    if type(post_request_processed) is str:
        return post_request_processed

    filepath = post_request_processed

    # Preprocess the json input into a gct file
    ssgsea_input = ssgsea.preprocess_ssgsea(filepath, ssgsea_type == 'gc')

    ssgsea_combined_output = ssgsea.run_ssgsea(ssgsea_input, ssgsea_type, ssc_input_type)

    # TODO: Delete everything when done to not accumulate files? Or: Keep the final output file, and use it as a cache?
    return send_response(send_file(ssgsea.postprocess_ssgsea(ssgsea_combined_output), as_attachment=True))


@app.route('/ksea', methods=['POST'])
@app.route('/ksea/<string:ksea_type>', methods=['POST'])
def handle_ksea_request(ksea_type=None) -> werkzeug.wrappers.Response | str:
    post_request_processed = process_post_request(request)

    if type(post_request_processed) is str:
        return post_request_processed

    filepath = post_request_processed

    preprocessed_filepath = ksea.preprocess_ksea(filepath)
    if ksea_type == 'rokai':
        preprocessed_filepath = ksea.run_rokai(preprocessed_filepath)

    ksea_result = ksea.perform_ksea(preprocessed_filepath)
    # TODO: Delete everything when done to not accumulate files? Or: Keep the final output file, and use it as a cache?
    return send_response(send_file(ksea_result, as_attachment=True))


def process_post_request(post_request: werkzeug.Request) -> Path | str:
    request_url = urlparse(request.base_url)

    form = post_request.form
    required_parameters = ['session_id', 'dataset_name']
    for param in required_parameters:
        if param not in form:
            return f'Error: parameter {param} not specified.\n'

    output_dir = Path('..') / secure_filename(form['session_id']) / secure_filename(
        form['dataset_name'] + request_url.path.replace('/', '_'))
    Path.mkdir(output_dir, parents=True, exist_ok=True)
    input_filepath = output_dir / 'input.json'
    # Check if the POST request has the file part, and else if it has the data part
    if 'file' in post_request.files and post_request.files['file'].filename != '':
        file = post_request.files['file']
        file.save(input_filepath)
    elif 'data' in post_request.form:
        # TODO: This variant is not tested yet
        with open(input_filepath, 'w') as o:
            o.write(post_request.form['data'])
    else:
        return "Error: You must either provide the input data " \
               + "as a JSON string (-F data=<JSON_String>) or as a file (-F file=@<Filepath>).\n"

    return input_filepath


def send_response(result) -> flask.Response:
    response = make_response(result)
    # TODO: I added this for cross-origin resource sharing, but is it unsafe?
    # Maybe using flask-cors (https://flask-cors.readthedocs.io/en/latest/)
    response.headers.add('Access-Control-Allow-Origin', '*')
    return response


if __name__ == '__main__':
    # TODO: What is the best way to switch debug off for deployment?
    app.run(debug=True, host='0.0.0.0', port=4321)
