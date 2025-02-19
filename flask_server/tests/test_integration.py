from pathlib import Path
import json
import pytest
from enrichment_server import app as application, VERSION


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
            res['Overlap (Experiment01)'] == exp['Overlap (Experiment01)'] and
            res['Overlap (Experiment02)'] == exp['Overlap (Experiment02)'] and
            res['Percent Overlap (Experiment01)'] == exp['Percent Overlap (Experiment01)'] and
            res['Percent Overlap (Experiment02)'] == exp['Percent Overlap (Experiment02)']

            for res, exp in zip(self.actual_result, self.expected_result))

    def evaluate_ksea(self):
        assert len(self.actual_result) == len(self.expected_result) and all(
            res['Gene'] == exp['Gene'] and
            all(round(res[f'Score (Experiment_{i})'], 5) == round(exp[f'Score (Experiment_{i})'], 5) for i in
                [1, 2, 3]) and
            all(round(res[f'adj p-val (Experiment_{i})'], 5) == round(exp[f'adj p-val (Experiment_{i})'], 5) for i in
                [1, 2, 3])
            for res, exp in zip(self.actual_result, self.expected_result))

    def evaluate_phonemes(self):
        assert len(self.actual_result) == len(self.expected_result) and all(
            res['pathway'] == exp['pathway'] and
            set([node['geneNames'][0] for node in res['nodes']]) == set(
                [node['geneNames'][0] for node in exp['nodes']]) and
            set([node['uniprotAccs'][0] for node in res['nodes']]) == set(
                [node['uniprotAccs'][0] for node in exp['nodes']]) and
            len(res['links']) == len(exp['links'])
            for res, exp in zip(self.actual_result, self.expected_result)
        )

    def evaluate_motif_enrichment(self):
        assert len(self.actual_result) == len(self.expected_result) and all(
            res == exp for res, exp in zip(self.actual_result, self.expected_result)
        )

    def evaluate_kea3(self):
        for key in self.expected_result.keys():
            for ranktype in 'MeanRank', 'TopRank':
                assert len(self.actual_result[key][ranktype]) == len(self.expected_result[key][ranktype])
                for rank_actual, rank_expected in zip(self.actual_result[key][ranktype], self.expected_result[key][ranktype]):
                    assert rank_actual['TF'] == rank_expected['TF']
                    assert rank_actual['Score'] == rank_expected['Score']

    def evaluate_kstar(self):
        for phospho_type in ['ST', 'Y']:
            assert len(self.actual_result[phospho_type]) == len(self.expected_result[phospho_type])
            for actual_elem, expected_elem in zip(self.actual_result[phospho_type], self.expected_result[phospho_type]):
                assert actual_elem['Kinase'] == expected_elem['Kinase']
                assert round(actual_elem['Experiment01'], 5) == round(expected_elem['Experiment01'], 5)
                assert round(actual_elem['Experiment02'], 5) == round(expected_elem['Experiment02'], 5)

    def test_get_status(self, client):
        response = client.get('/')
        # You can do either of the following
        assert response.status == '200 OK', response.status
        assert response.json == {'status': 200, 'version': VERSION}, response.json

    def test_ssgsea_ssc_flanking(self, client):
        self.input_json = Path('../fixtures/ptm-sea/input/input_flanking.json')
        self.dataset_name = 'ptmsea_test'

        response = client.post('/ssgsea/ssc/flanking', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ptm-sea/expected_output/output_flanking.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ssgsea()

    def test_ssgsea_ssc_uniprot(self, client):
        self.input_json = Path('../fixtures/ptm-sea/input/input_uniprot.json')
        self.dataset_name = 'ptmsea_test'

        response = client.post('/ssgsea/ssc/uniprot', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ptm-sea/expected_output/output_uniprot.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ssgsea()

    def test_ssgsea_gc(self, client):
        self.input_json = Path('../fixtures/ssgsea/input/input.json')
        self.dataset_name = 'ssgsea_gc_test'

        response = client.post('/ssgsea/gc', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ssgsea/expected_output/output_gc.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ssgsea()

    def test_ssgsea_gcr(self, client):
        self.input_json = Path('../fixtures/ssgsea/input/input.json')
        self.dataset_name = 'ssgsea_gcr_test'

        response = client.post('/ssgsea/gcr', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ssgsea/expected_output/output_gcr.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ssgsea()

    def test_ksea(self, client):
        self.input_json = Path('../fixtures/ksea/input/input.json')
        self.dataset_name = 'ksea_test'

        response = client.post('/ksea', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ksea/expected_output/output_ksea.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ksea()

    def test_ksea_rokai(self, client):
        self.input_json = Path('../fixtures/ksea/input/input.json')
        self.dataset_name = 'ksea_rokai_test'

        response = client.post('/ksea/rokai', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ksea/expected_output/output_ksea_rokai.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ksea()

    def test_motif_enrichment(self, client):
        self.input_json = Path('../fixtures/motif_enrichment/input/input.json')
        self.dataset_name = 'motif_enrichment_test'

        response = client.post('/motif_enrichment', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/motif_enrichment/expected_output/output.json')
        self.expected_result = json.load(open(expected_result_file))
        self.evaluate_motif_enrichment()

    def test_kea3(self, client):
        self.input_json = Path('../fixtures/kea3/input/input.json')
        self.dataset_name = 'kea3_test'

        response = client.post('/kea3', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/kea3/expected_output/output.json')
        self.expected_result = json.load(open(expected_result_file))
        self.evaluate_kea3()

    def test_kstar(self, client):
        self.input_json = Path('../fixtures/kstar/input/input.json')
        self.dataset_name = 'kstar_test'

        response = client.post('/kstar', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/kstar/expected_output/output.json')
        self.expected_result = json.load(open(expected_result_file))
        self.evaluate_kstar()

#Run PHONEMeS last because it takes the longest
    def test_phonemes(self, client):
        self.input_json = Path('../fixtures/phonemes/input/input.json')
        self.dataset_name = 'phonemes_test'

        response = client.post('/phonemes', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/phonemes/expected_output/json_skeletons.json')
        self.expected_result = json.load(open(expected_result_file))
        self.evaluate_phonemes()