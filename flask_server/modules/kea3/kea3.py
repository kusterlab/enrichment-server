import json
import requests
from pathlib import Path

KEA3_URL = 'https://amp.pharm.mssm.edu/kea3/api/enrich/'


def run_kea3_api(filepath: Path) -> Path:
    input_json = json.load(open(filepath))
    result = dict()
    for experiment in input_json.keys():
        payload = {'gene_set': input_json[experiment], 'query_name': experiment}
        response = requests.post(KEA3_URL, data=json.dumps(payload))
        if not response.ok:
            raise Exception(f'Error in Enrichment of Experiment {experiment}')
        response_json = json.loads(response.text)
        result[experiment] = {
            'MeanRank': response_json['Integrated--meanRank'],
            'TopRank': response_json['Integrated--topRank']
        }

    output_json = filepath.parent / f'kea3_result.json'
    with open(output_json, 'w') as outfile:
        json.dump(result, outfile)

    return output_json
