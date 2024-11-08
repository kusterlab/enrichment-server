# How to send a request:
# curl -X POST -F file=@<input_file> -F session_id=ABCDEF12345  -F dataset_name=FooBar http://127.0.0.1:1234/<route>
import os
from pathlib import Path
import shutil
import json
from urllib.parse import urlparse
import werkzeug.wrappers
from werkzeug.utils import secure_filename
from flask import Flask, request, send_file, jsonify, make_response
import flask.wrappers
import logging
import sys

from modules.ssgsea import ssgsea
from modules.ksea import ksea
from modules.phonemes import phonemes
from modules.motif_enrichment import motif_enrichment
from modules.kea3 import kea3
from modules.k_star import k_star

VERSION = '0.1.3'


def setup_logger():
    global LOGGER
    LOGGER = logging.getLogger('enrichment_server')
    LOGGER.setLevel(logging.INFO)

    file_handler = logging.FileHandler('enrichment_server_logfile.log')
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    file_handler.setFormatter(formatter)
    LOGGER.addHandler(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    LOGGER.addHandler(console_handler)

    # Redirect print statements to the logger
    class LoggerWriter:
        def __init__(self, level):
            self.level = level

        def write(self, message):
            if message.strip():  # Avoid logging empty messages
                self.level(message)

        def flush(self):
            pass  # No-op for compatibility with sys.stdout/err

    # Replace stdout and stderr with logger
    sys.stdout = LoggerWriter(LOGGER.info)
    sys.stderr = LoggerWriter(LOGGER.error)


def create_app():
    setup_logger()
    app = Flask(__name__)
    print('App created.')
    return app


app = create_app()


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

    post_request_processed = process_post_request(request, f'ssGSEA ({ssgsea_type.upper()})')

    if type(post_request_processed) is str:
        return post_request_processed

    filepath = post_request_processed

    # Preprocess the json input into a gct file
    ssgsea_input = ssgsea.preprocess_ssgsea(filepath, ssgsea_type != 'gcr')

    ssgsea_combined_output = ssgsea.run_ssgsea(ssgsea_input, ssgsea_type, ssc_input_type)

    return send_response(postprocess_request_response(ssgsea.postprocess_ssgsea(ssgsea_combined_output),
                                                      f'ssGSEA ({ssgsea_type.upper()})'),
                         filepath.parent)


@app.route('/ksea', methods=['POST'])
@app.route('/ksea/<string:ksea_type>', methods=['POST'])
def handle_ksea_request(ksea_type=None) -> werkzeug.wrappers.Response | str:
    post_request_processed = process_post_request(request, 'KSEA' if not ksea_type else 'RoKAI+KSEA')

    if type(post_request_processed) is str:
        return post_request_processed

    filepath = post_request_processed

    preprocessed_filepath = ksea.preprocess_ksea(filepath)
    if ksea_type == 'rokai':
        preprocessed_filepath = ksea.run_rokai(preprocessed_filepath)

    ksea_result = ksea.perform_ksea(preprocessed_filepath)
    return send_response(postprocess_request_response(
        ksea_result, 'KSEA' if not ksea_type else 'RoKAI+KSEA'),
        filepath.parent)


@app.route('/phonemes', methods=['POST'])
def handle_phonemes_request() -> werkzeug.wrappers.Response | str:
    post_request_processed = process_post_request(request, 'PHONEMeS')

    if type(post_request_processed) is str:
        return post_request_processed

    filepath = post_request_processed

    preprocessed_filepath = phonemes.preprocess_phonemes(filepath)

    phonemes_result = phonemes.run_phonemes(preprocessed_filepath)
    cytoscape_result = phonemes.run_cytoscape(phonemes_result)
    pathway_skeletons_json = phonemes.create_pathway_skeleton(cytoscape_result)

    return send_response(postprocess_request_response(pathway_skeletons_json, 'PHONEMeS'), filepath.parent)


@app.route('/motif_enrichment', methods=['POST'])
def handle_motif_enrichment_request() -> werkzeug.wrappers.Response | str:
    post_request_processed = process_post_request(request, 'Motif Enrichment')

    if type(post_request_processed) is str:
        return post_request_processed

    filepath = post_request_processed
    motif_enrichment_result = motif_enrichment.run_motif_enrichment(filepath)

    return send_response(postprocess_request_response(motif_enrichment_result, 'Motif Enrichment'), filepath.parent)


@app.route('/kea3', methods=['POST'])
def handle_kea3_request() -> werkzeug.wrappers.Response | str:
    post_request_processed = process_post_request(request, 'KEA3')

    if type(post_request_processed) is str:
        return post_request_processed

    filepath = post_request_processed
    kea3_result = kea3.run_kea3_api(filepath)

    return send_response(postprocess_request_response(kea3_result, 'KEA3'), filepath.parent)


@app.route('/kstar', methods=['POST'])
def handle_kstar_request() -> werkzeug.wrappers.Response | str:
    post_request_processed = process_post_request(request, 'KSTAR')

    if type(post_request_processed) is str:
        return post_request_processed

    filepath = post_request_processed
    kstar_result = k_star.run_kstar(filepath)

    return send_response(postprocess_request_response(kstar_result, 'KSTAR'), filepath.parent)


def process_post_request(post_request: werkzeug.Request, method: str) -> Path | str:
    print(f"{method} request received.")
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


def postprocess_request_response(result_path: Path, method: str) -> werkzeug.wrappers.Response:
    result_raw = json.load(open(result_path))
    result_with_log = {'Log': {'Version': VERSION}, 'Result': result_raw}
    with open(result_path, 'w') as outfile:
        json.dump(result_with_log, outfile)
    print(f"{method} analysis completed successfully.")
    return send_file(result_path, as_attachment=False)


def send_response(result: werkzeug.wrappers.Response, output_folder=None) -> flask.Response:
    response = make_response(result)
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

    app.run(debug=os.getenv("PRODUCTION", '0') != '1', host='0.0.0.0', port=int(os.getenv("PORT", '4321')))
