from pathlib import Path
import json
import pytest
import pandas as pd
from enrichment_server import app as application, VERSION
import io


@pytest.fixture()
def app():
    yield application


@pytest.fixture()
def client(app):
    return app.test_client()


class TestClass:
    session_id = 'TESTSESSION'
    dataset_name = None
    input_file = None
    parameters_toml = None
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

    def evaluate_ssgsea_csv(self):
        pd.testing.assert_frame_equal(self.actual_result[
                                          ['Signature ID', 'Overlap (Experiment01)', 'Overlap (Experiment02)',
                                           'Percent Overlap (Experiment01)', 'Percent Overlap (Experiment02)']],
                                      self.expected_result[
                                          ['Signature ID', 'Overlap (Experiment01)', 'Overlap (Experiment02)',
                                           'Percent Overlap (Experiment01)', 'Percent Overlap (Experiment02)']])

    def evaluate_go_enrichment(self):
        assert len(self.actual_result) == len(self.expected_result) and all(
            res['GO_ID'] == exp['GO_ID'] and
            res['Intersection_Size (Experiment01)'] == exp['Intersection_Size (Experiment01)'] and
            ((res['Intersection (Experiment02)'] is None and exp['Intersection (Experiment02)'] is None) or
             set(res['Intersection (Experiment02)'].split(',')) ==
             set(exp['Intersection (Experiment02)'].split(','))) and
            res['neg_log10_adjusted_p_value (Experiment03)'] == exp['neg_log10_adjusted_p_value (Experiment03)']

            for res, exp in zip(self.actual_result, self.expected_result))

    def evaluate_ksea(self):
        assert len(self.actual_result) == len(self.expected_result) and all(
            res['Gene'] == exp['Gene'] and
            all(round(res[f'Score (Experiment_{i})'], 5) == round(exp[f'Score (Experiment_{i})'], 5) for i in
                [1, 2, 3]) and
            all(round(res[f'adj p-val (Experiment_{i})'], 5) == round(exp[f'adj p-val (Experiment_{i})'], 5) for i in
                [1, 2, 3])
            for res, exp in zip(self.actual_result, self.expected_result))

    def evaluate_ksea_csv(self):
        pd.testing.assert_series_equal(self.actual_result['Gene'], self.expected_result['Gene'])
        pd.testing.assert_series_equal(self.actual_result['Score (Experiment_1)'],
                                       self.expected_result['Score (Experiment_1)'], check_exact=False,
                                       atol=1e-4)
        pd.testing.assert_series_equal(self.actual_result['adj p-val (Experiment_2)'],
                                       self.expected_result['adj p-val (Experiment_2)'], check_exact=False,
                                       atol=1e-4)
        pd.testing.assert_series_equal(
            self.actual_result['Overlap (Experiment_3)'].apply(lambda s: json.loads(s.replace("'", '"'))).apply(
                set),
            self.expected_result['Overlap (Experiment_3)'].apply(
                lambda s: json.loads(s.replace("'", '"'))).apply(set))

    def evaluate_rokai(self):
        assert len(self.actual_result) == len(self.expected_result) and all(
            res['Gene'] == exp['Gene'] and
            all(round(res[f'Activity (Quantification{suffix})'], 5) == round(exp[f'Activity (Quantification{suffix})'],
                                                                             5) for suffix in
                ['', '.2']) and
            all(round(res[f'FDR (Quantification{suffix})'], 5) == round(
                exp[f'FDR (Quantification{suffix})'], 5) for suffix in
                ['', '.2'])
            for res, exp in zip(self.actual_result, self.expected_result))

    def evaluate_rokai_csv(self):
        pd.testing.assert_series_equal(self.actual_result['Gene'], self.expected_result['Gene'])
        pd.testing.assert_series_equal(self.actual_result['Activity (Experiment_1)'],
                                       self.expected_result['Activity (Experiment_1)'], check_exact=False,
                                       atol=1e-4)
        pd.testing.assert_series_equal(self.actual_result['FDR (Experiment_3)'],
                                       self.expected_result['FDR (Experiment_3)'], check_exact=False,
                                       atol=1e-4)
        pd.testing.assert_series_equal(self.actual_result['ZScore (Experiment_3)'],
                                       self.expected_result['ZScore (Experiment_3)'], check_exact=False,
                                       atol=1e-4)

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

    def evaluate_motif_enrichment_csv(self):
        pd.testing.assert_frame_equal(self.actual_result, self.expected_result)

    def evaluate_kea3(self):
        for key in self.expected_result.keys():
            for ranktype in 'MeanRank', 'TopRank':
                assert len(self.actual_result[key][ranktype]) == len(self.expected_result[key][ranktype])
                for rank_actual, rank_expected in zip(self.actual_result[key][ranktype],
                                                      self.expected_result[key][ranktype]):
                    assert rank_actual['TF'] == rank_expected['TF']
                    assert rank_actual['Score'] == rank_expected['Score']

    def evaluate_kstar_json(self):
        for phospho_type in ['ST', 'Y']:
            assert len(self.actual_result[phospho_type]) == len(self.expected_result[phospho_type])
            for actual_elem, expected_elem in zip(self.actual_result[phospho_type], self.expected_result[phospho_type]):
                assert actual_elem['Kinase'] == expected_elem['Kinase']
                if 'up (Experiment01)' in actual_elem:
                    assert round(actual_elem['up (Experiment01)'], 5) == round(expected_elem['up (Experiment01)'], 5)
                else:
                    assert round(actual_elem['down (Experiment02)'], 5) == round(expected_elem['down (Experiment02)'],
                                                                                 5)

    def evaluate_kstar_csv(self):
        pd.testing.assert_frame_equal(self.actual_result, self.expected_result,
                                      check_exact=False, atol=1e-4)

    def test_get_status(self, client):
        response = client.get('/')
        # You can do either of the following
        assert response.status == '200 OK', response.status
        assert response.json == {'status': 200, 'version': VERSION}, response.json

    def test_ssgsea_ssc_flanking(self, client):
        self.input_file = Path('../fixtures/ptm-sea/input/input_flanking_hsa.json')
        self.dataset_name = 'ptmsea_test'

        response = client.post('/ssgsea/ssc/flanking', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ptm-sea/expected_output/output_flanking_hsa.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ssgsea()

    def test_ssgsea_ssc_flanking_w_parameters(self, client):
        self.input_file = Path('../fixtures/ptm-sea/input/input_flanking_hsa.json')
        self.parameters_toml = Path('../fixtures/ptm-sea/input/parameters.toml')
        self.dataset_name = 'ptmsea_test_w_parameters'

        response = client.post('/ssgsea/ssc/flanking', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb'),
            "parameters": self.parameters_toml.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ptm-sea/expected_output/output_flanking_w_parameters.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ssgsea()

    def test_ssgsea_ssc_uniprot(self, client):
        self.input_file = Path('../fixtures/ptm-sea/input/input_uniprot_hsa.json')
        self.dataset_name = 'ptmsea_test'

        response = client.post('/ssgsea/ssc/uniprot', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ptm-sea/expected_output/output_uniprot_hsa.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ssgsea()

    def test_ssgsea_ssc_flanking_csv(self, client):
        self.input_file = Path('../fixtures/ptm-sea/input/input_flanking_hsa.txt')
        self.dataset_name = 'ptmsea_test_csv'

        response = client.post('/ssgsea/ssc/flanking', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })
        self.actual_result = pd.read_csv(io.BytesIO(response.data), sep='\t')
        expected_result_file = Path('../fixtures/ptm-sea/expected_output/ptmsea_output_hsa.txt')
        self.expected_result = pd.read_csv(expected_result_file, sep='\t')
        self.evaluate_ssgsea_csv()

    def test_ssgsea_gc(self, client):
        self.input_file = Path('../fixtures/ssgsea/input/input_hsa.json')
        self.dataset_name = 'ssgsea_gc_test'

        response = client.post('/ssgsea/gc', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ssgsea/expected_output/output_gc_hsa.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ssgsea()

    def test_ssgsea_gcr(self, client):
        self.input_file = Path('../fixtures/ssgsea/input/input_hsa.json')
        self.dataset_name = 'ssgsea_gcr_test'

        response = client.post('/ssgsea/gcr', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ssgsea/expected_output/output_gcr_hsa.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ssgsea()

    def test_ssgsea_mouse_gc(self, client):
        self.input_file = Path('../fixtures/ssgsea/input/input_mmu.txt')
        self.dataset_name = 'ssgsea_mouse_gc_test'

        response = client.post('/ssgsea/gc', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb'),
            "organism": 'mmu'
        })

        self.actual_result = pd.read_csv(io.BytesIO(response.data), sep='\t')
        expected_result_file = Path('../fixtures/ssgsea/expected_output/output_gc_mmu.txt')
        self.expected_result = pd.read_csv(expected_result_file, sep='\t')
        self.evaluate_ssgsea_csv()

    def test_ssgsea_mouse_gcr(self, client):
        self.input_file = Path('../fixtures/ssgsea/input/input_mmu.txt')
        self.dataset_name = 'ssgsea_mouse_gcr_test'

        response = client.post('/ssgsea/gcr', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb'),
            "organism": 'mmu'
        })

        self.actual_result = pd.read_csv(io.BytesIO(response.data), sep='\t')
        expected_result_file = Path('../fixtures/ssgsea/expected_output/output_gcr_mmu.txt')
        self.expected_result = pd.read_csv(expected_result_file, sep='\t')
        self.evaluate_ssgsea_csv()

    def test_ssgsea_mouse_ssc_flanking(self, client):
        self.input_file = Path('../fixtures/ptm-sea/input/input_flanking_mmu.txt')
        self.dataset_name = 'ptmsea_mouse_test'

        response = client.post('/ssgsea/ssc/flanking', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb'),
            "organism": 'mmu'
        })

        self.actual_result = pd.read_csv(io.BytesIO(response.data), sep='\t')
        expected_result_file = Path('../fixtures/ptm-sea/expected_output/ptmsea_output_mmu.txt')
        self.expected_result = pd.read_csv(expected_result_file, sep='\t')
        self.evaluate_ssgsea_csv()

    def test_go_enrichment_human(self, client):
        self.input_file = Path('../fixtures/go_enrichment/input/input_hsa.json')
        self.dataset_name = 'go_enrichment_test'

        response = client.post('/go_enrichment', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/go_enrichment/expected_output/go_enrichment_result_hsa.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_go_enrichment()

    def test_go_enrichment_w_custom_background(self, client):
        self.input_file = Path('../fixtures/go_enrichment/input/input_w_background.json')
        self.dataset_name = 'go_enrichment_w_custom_background_test'

        response = client.post('/go_enrichment', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path(
            '../fixtures/go_enrichment/expected_output/go_enrichment_w_custom_background_result.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_go_enrichment()

    def test_go_enrichment_bad_organism(self, client):
        self.input_file = Path('../fixtures/go_enrichment/input/input_hsa.json')
        self.dataset_name = 'go_enrichment_bad_organism_test'

        response = client.post('/go_enrichment', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb'),
            'organism': 'lol'
        })

        assert response.status_code == 400

    def test_go_enrichment_mouse(self, client):
        self.input_file = Path('../fixtures/go_enrichment/input/input_mmu.json')
        self.dataset_name = 'go_enrichment_mouse_test'

        response = client.post('/go_enrichment', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb'),
            "organism": "mmu"
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/go_enrichment/expected_output/go_enrichment_result_mmu.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_go_enrichment()

    def test_ksea_json(self, client):
        self.input_file = Path('../fixtures/ksea/input/input_hsa.json')
        self.dataset_name = 'ksea_test_json'

        response = client.post('/ksea', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ksea/expected_output/output_ksea.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ksea()

    def test_ksea_csv(self, client):
        self.input_file = Path('../fixtures/ksea/input/input_hsa.csv')
        self.dataset_name = 'ksea_test_csv'

        response = client.post('/ksea', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        self.actual_result = pd.read_csv(io.BytesIO(response.data), sep='\t')
        expected_result_file = Path('../fixtures/ksea/expected_output/output_ksea_hsa.txt')
        self.expected_result = pd.read_csv(expected_result_file, sep='\t')
        self.evaluate_ksea_csv()

    def test_ksea_w_parameters(self, client):
        self.input_file = Path('../fixtures/ksea/input/input_hsa.json')
        self.parameters_toml = Path('../fixtures/ksea/input/parameters.toml')
        self.dataset_name = 'ksea_test_w_parameters'

        response = client.post('/ksea', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb'),
            "parameters": self.parameters_toml.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ksea/expected_output/output_ksea_w_parameters.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ksea()

    def test_ksea_mouse(self, client):
        self.input_file = Path('../fixtures/ksea/input/input_mmu.csv')
        self.dataset_name = 'ksea_test_mouse'

        response = client.post('/ksea', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb'),
            "organism": 'mmu'
        })

        self.actual_result = pd.read_csv(io.BytesIO(response.data), sep='\t')
        expected_result_file = Path('../fixtures/ksea/expected_output/output_ksea_mmu.txt')
        self.expected_result = pd.read_csv(expected_result_file, sep='\t')
        self.evaluate_ksea_csv()

    def test_ksea_rokai(self, client):
        self.input_file = Path('../fixtures/ksea/input/input_hsa.json')
        self.dataset_name = 'ksea_rokai_test'

        response = client.post('/ksea/rokai', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb'),
            'organism': 'hsa'
        })

    def test_ksea_rokai_w_parameters(self, client):
        self.input_file = Path('../fixtures/ksea/input/input_hsa.json')
        self.parameters_toml = Path('../fixtures/ksea/input/parameters_ksea+rokai.toml')
        self.dataset_name = 'ksea_rokai_test_w_parameters'

        response = client.post('/ksea/rokai', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb'),
            'organism': 'hsa',
            "parameters": self.parameters_toml.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/ksea/expected_output/output_ksea_rokai_w_parameters.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_ksea()

    def test_rokai_hsa_json(self, client):
        self.input_json = Path('../fixtures/rokai/input/input_hsa.json')
        self.dataset_name = 'rokai_test'

        response = client.post('/rokai', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/rokai/expected_output/rokai_result_hsa.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_rokai()

    def test_rokai_mmu_csv(self, client):
        self.input_json = Path('../fixtures/rokai/input/input_mmu.csv')
        self.dataset_name = 'rokai_test_mmu'

        response = client.post('/rokai', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb'),
            'organism': 'mmu'
        })

        self.actual_result = pd.read_csv(io.BytesIO(response.data), sep='\t')
        expected_result_file = Path('../fixtures/rokai/expected_output/rokai_result_mmu.txt')
        self.expected_result = pd.read_csv(expected_result_file, sep='\t')
        self.evaluate_rokai_csv()

    def test_rokai_w_parameters(self, client):
        self.input_json = Path('../fixtures/rokai/input/input_hsa.json')
        self.dataset_name = 'rokai_test_w_parameters'
        self.parameters_toml = Path('../fixtures/rokai/input/parameters.toml')

        response = client.post('/rokai', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_json.open('rb'),
            "parameters": self.parameters_toml.open('rb')

        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/rokai/expected_output/rokai_result_w_parameters.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_rokai()

    def test_motif_enrichment(self, client):
        self.input_file = Path('../fixtures/motif_enrichment/input/input.json')
        self.dataset_name = 'motif_enrichment_test'

        response = client.post('/motif_enrichment', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/motif_enrichment/expected_output/output.json')
        self.expected_result = json.load(open(expected_result_file))
        self.evaluate_motif_enrichment()

    def test_motif_enrichment_w_parameters(self, client):
        self.input_file = Path('../fixtures/motif_enrichment/input/input.json')
        self.parameters_toml = Path('../fixtures/motif_enrichment/input/parameters.toml')
        self.dataset_name = 'motif_enrichment_test_w_parameters'

        response = client.post('/motif_enrichment', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb'),
            "parameters": self.parameters_toml.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/motif_enrichment/expected_output/output_w_parameters.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_motif_enrichment()

    def test_motif_enrichment_csv(self, client):
        self.input_file = Path('../fixtures/motif_enrichment/input/input.txt')
        self.dataset_name = 'motif_enrichment_test_csv'

        response = client.post('/motif_enrichment', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })
        self.actual_result = pd.read_csv(io.BytesIO(response.data), sep='\t')
        expected_result_file = Path('../fixtures/motif_enrichment/expected_output/output.txt')
        self.expected_result = pd.read_csv(expected_result_file, sep='\t')
        self.evaluate_motif_enrichment_csv()

    def test_kea3(self, client):
        self.input_file = Path('../fixtures/kea3/input/input.json')
        self.dataset_name = 'kea3_test'

        response = client.post('/kea3', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/kea3/expected_output/output.json')
        self.expected_result = json.load(open(expected_result_file))
        self.evaluate_kea3()

    def test_kstar_json(self, client):
        self.input_file = Path('../fixtures/kstar/input/input.json')
        self.dataset_name = 'kstar_test_json'

        response = client.post('/kstar', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/kstar/expected_output/output.json')
        self.expected_result = json.load(open(expected_result_file))['Result']
        self.evaluate_kstar_json()

    def test_kstar_csv(self, client):
        self.input_file = Path('../fixtures/kstar/input/input.csv')
        self.dataset_name = 'kstar_test_csv'

        response = client.post('/kstar', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        self.actual_result = pd.read_csv(io.BytesIO(response.data), sep='\t')
        expected_result_file = Path('../fixtures/kstar/expected_output/output.txt')
        self.expected_result = pd.read_csv(expected_result_file, sep='\t')
        self.evaluate_kstar_csv()

    def test_kstar_w_parameters(self, client):
        self.input_file = Path('../fixtures/kstar/input/input.csv')
        self.parameters_toml = Path('../fixtures/kstar/input/parameters.toml')

        self.dataset_name = 'kstar_test_w_parameters'

        response = client.post('/kstar', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb'),
            "parameters": self.parameters_toml.open('rb')
        })

        self.actual_result = pd.read_csv(io.BytesIO(response.data), sep='\t')
        expected_result_file = Path('../fixtures/kstar/expected_output/output_w_parameters.txt')
        self.expected_result = pd.read_csv(expected_result_file, sep='\t')
        self.evaluate_kstar_csv()

    # Run PHONEMeS last because it takes the longest
    def test_phonemes(self, client):
        self.input_file = Path('../fixtures/phonemes/input/input.json')
        self.dataset_name = 'phonemes_test'

        response = client.post('/phonemes', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        self.actual_result = json.loads(response.data)['Result']
        expected_result_file = Path('../fixtures/phonemes/expected_output/json_skeletons.json')
        self.expected_result = json.load(open(expected_result_file))
        self.evaluate_phonemes()

    def test_nofile(self, client):
        self.dataset_name = 'no_file'

        response = client.post('/phonemes', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name
        })
        assert response.status_code == 400

    def test_forbidden_csv_inputs(self, client):
        # Use KSEA Input, but it could be any file ending in .csv
        self.input_file = Path('../fixtures/ksea/input/input_hsa.csv')
        self.dataset_name = 'forbidden_csv_test'

        # Call three endpoints that have no csv
        phonemes_response = client.post('/phonemes', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        assert phonemes_response.status_code == 400

        # Call three endpoints that have no csv
        kea3_response = client.post('/kea3', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        assert kea3_response.status_code == 400

        # Call three endpoints that have no csv
        go_enrichment_response = client.post('/go_enrichment', data={
            "session_id": self.session_id,
            "dataset_name": self.dataset_name,
            "file": self.input_file.open('rb')
        })

        assert go_enrichment_response.status_code == 400
