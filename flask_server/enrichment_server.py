# How to send a request:
# curl -X POST -F file=@<input_file> -F session_id=ABCDEF12345  -F dataset_name=FooBar http://127.0.0.1:1234/<route>
from pathlib import Path
import shutil
from urllib.parse import urlparse
import werkzeug.wrappers
from werkzeug.utils import secure_filename
from flask import Flask, request, send_file, jsonify, make_response
import flask.wrappers
import ssgsea
import ksea
import phonemes
import motif_enrichment

app = Flask(__name__)

VERSION = '0.0.3'


@app.route('/', methods=['GET'])
def get_status() -> flask.wrappers.Response:
    return send_response(jsonify(status=200, version=VERSION))


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
    ssgsea_input = ssgsea.preprocess_ssgsea(filepath, ssgsea_type != 'gcr')

    ssgsea_combined_output = ssgsea.run_ssgsea(ssgsea_input, ssgsea_type, ssc_input_type)

    return send_response(send_file(ssgsea.postprocess_ssgsea(ssgsea_combined_output), as_attachment=False),
                         filepath.parent)


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
    return send_response(send_file(ksea_result, as_attachment=False), filepath.parent)


@app.route('/phonemes', methods=['POST'])
def handle_phonemes_request() -> werkzeug.wrappers.Response | str:
    post_request_processed = process_post_request(request)

    if type(post_request_processed) is str:
        return post_request_processed

    filepath = post_request_processed

    preprocessed_filepath = phonemes.preprocess_phonemes(filepath)

    phonemes_result = phonemes.run_phonemes(preprocessed_filepath)
    cytoscape_result = phonemes.run_cytoscape(phonemes_result)
    pathway_skeletons_json = phonemes.create_pathway_skeleton(cytoscape_result)

    return send_response(send_file(pathway_skeletons_json, as_attachment=False), filepath.parent)


@app.route('/motif_enrichment', methods=['POST'])
def handle_motif_enrichment_request() -> werkzeug.wrappers.Response | str:
    post_request_processed = process_post_request(request)

    if type(post_request_processed) is str:
        return post_request_processed

    filepath = post_request_processed
    motif_enrichment_result = motif_enrichment.run_motif_enrichment(filepath)

    return send_response(send_file(motif_enrichment_result, as_attachment=False), filepath.parent)


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


def send_response(result: werkzeug.wrappers.Response, output_folder=None) -> flask.Response:
    response = make_response(result)
    # TODO: I added this for cross-origin resource sharing, but is it unsafe?
    # Maybe using flask-cors (https://flask-cors.readthedocs.io/en/latest/)
    response.headers.add('Access-Control-Allow-Origin', '*')
    if output_folder:
        shutil.rmtree(output_folder)
    return response


if __name__ == '__main__':
    # # DEBUG: Limit memory usage
    # import resource
    # import sys
    #
    #
    # def memory_limit_half():
    #     """Limit max memory usage to half."""
    #     soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    #     # Convert KiB to bytes, and divide in two to half
    #     resource.setrlimit(resource.RLIMIT_AS, (int(get_memory() * 1024 / 2), hard))
    #
    #
    # def get_memory():
    #     with open('/proc/meminfo', 'r') as mem:
    #         free_memory = 0
    #         for i in mem:
    #             sline = i.split()
    #             if str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
    #                 free_memory += int(sline[1])
    #     return free_memory  # KiB
    #
    #
    # # GUBED
    # memory_limit_half()

    # TODO: What is the best way to switch debug off for deployment?
    app.run(debug=True, host='0.0.0.0', port=4321)
