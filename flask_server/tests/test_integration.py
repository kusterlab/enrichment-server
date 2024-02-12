from pathlib import Path
import json
import pytest
from enrichment_server import app as application


@pytest.fixture()
def app():
    yield application


@pytest.fixture()
def client(app):
    return app.test_client()


class TestClass:
    session_id = 'TESTSESSION'
    dataset_name = None
    input_json = None
    actual_result = None
    expected_result = None

    def evaluate_ssgsea(self):
        assert len(self.actual_result) == len(self.expected_result) and all(
            res['Signature ID'] == exp['Signature ID'] and
            res['Overlap [%] (Experiment01)'] == exp['Overlap [%] (Experiment01)'] and
            res['Overlap [%] (Experiment02)'] == exp['Overlap [%] (Experiment02)']
            for res, exp in zip(self.actual_result, self.expected_result))

    def test_get_status(self, client):
        response = client.get('/')
        # You can do either of the following
        assert response.status == '200 OK', response.status
        assert response.json == {'status': 200}, response.json

    def test_ssgsea_ssc(self, client):
        self.input_json = Path('../fixtures/ptm-sea/input/input.json')
        self.dataset_name = 'ptmsea_test'

        response = client.post('/ssgsea/ssc', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)
        expected_result_file = Path('../fixtures/ptm-sea/expected_output/output_indented.json')
        self.expected_result = json.load(open(expected_result_file))
        self.evaluate_ssgsea()

    def test_ssgsea_gc(self, client):
        self.input_json = Path('../fixtures/ssgsea/input/input.json')
        self.dataset_name = 'ssgsea_gc_test'

        response = client.post('/ssgsea/gc', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)
        expected_result_file = Path('../fixtures/ssgsea/expected_output/output_gc.json')
        self.expected_result = json.load(open(expected_result_file))
        self.evaluate_ssgsea()

    def test_ssgsea_gcr(self, client):
        self.input_json = Path('../fixtures/ssgsea/input/input.json')
        self.dataset_name = 'ssgsea_gcr_test'

        response = client.post('/ssgsea/gcr', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)
        expected_result_file = Path('../fixtures/ssgsea/expected_output/output_gcr.json')
        self.expected_result = json.load(open(expected_result_file))
        self.evaluate_ssgsea()

    def test_ksea(self, client):
        self.input_json = Path('../fixtures/ksea/input/input.json')
        self.dataset_name = 'ksea_test'

        response = client.post('/ksea', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)
        expected_result_file = Path('../fixtures/ksea/expected_output/output_ksea.json')
        self.expected_result = json.load(open(expected_result_file))
        assert self.actual_result == self.expected_result

    def test_ksea_rokai(self, client):
        self.input_json = Path('../fixtures/ksea/input/input.json')
        self.dataset_name = 'ksea_rokai_test'

        response = client.post('/ksea/rokai', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)
        expected_result_file = Path('../fixtures/ksea/expected_output/output_ksea_rokai.json')
        self.expected_result = json.load(open(expected_result_file))
        assert self.actual_result == self.expected_result
