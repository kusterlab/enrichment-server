import json
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
from typing import TypedDict, List, Dict
from requests.structures import CaseInsensitiveDict as CIDict


ORGANISM_TO_ANNOTATION_FILE = \
    {'hsa': '../db/go/gprofiler_hsapiens.GO.name.gmt',
     'mmu': '../db/go/gprofiler_mmusculus.GO.name.gmt'}

SUPPORTED_ORGANISMS: list[str] = list(ORGANISM_TO_ANNOTATION_FILE.keys())
GOAnnotation = TypedDict('GOAnnotation', {'ID': str, 'Name': str, 'Genes': list[str]})


def load_go_annotations(organism: SUPPORTED_ORGANISMS) -> list[GOAnnotation]:
    go_terms_file = Path(ORGANISM_TO_ANNOTATION_FILE.get(organism))
    go_terms = []
    with open(go_terms_file, 'r') as infile:
        line = infile.readline().strip()
        while line:
            linesplit = line.strip().split('\t')
            go_terms.append({'ID': linesplit[0], 'Name': linesplit[1], 'Genes': linesplit[2:]})
            line = infile.readline()

    return go_terms


def run_gene_ontology_fisher_test(query: set, annotation: GOAnnotation, background: set, experiment_name: str) -> dict:
    annotation_genes = set(annotation['Genes'])
    annotation_query_overlap = background.intersection(query.intersection(annotation_genes))
    # Check if we can speed it up with this:
    if not annotation_query_overlap:
        return None
    a = len(annotation_query_overlap)
    b = len(background.intersection(annotation_genes - query))
    c = len((background - annotation_genes).intersection(query))
    d = len((background - annotation_genes) - query)
    odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]])
    return {f'GO_ID': annotation['ID'],
            f'Name': annotation['Name'],
            f'Term_Size': len(background.intersection(annotation_genes)),
            f'p_value ({experiment_name})': p_value,
            f'Intersection_Size ({experiment_name})': len(annotation_query_overlap),
            f'Intersection ({experiment_name})': ",".join(annotation_query_overlap)}


def run_go_enrichment(filepath: Path, organism: SUPPORTED_ORGANISMS) -> Path:
    go_terms = load_go_annotations(organism)
    input_json = CIDict(json.load(open(filepath)))
    if 'background' in input_json:
        background = set(input_json.pop('background'))
    else:
        # If no Custom Background was provided, use all genes that appear in GO
        background = set(gene for go_term in go_terms for gene in go_term['Genes'])

    results = []
    for experiment_name in input_json.keys():
        query = set(input_json[experiment_name])
        all_go_enrichment_results = [run_gene_ontology_fisher_test(query, annotation, background, experiment_name) for
                                     annotation in go_terms]
        go_enrichment_results_df = pd.DataFrame([res for res in all_go_enrichment_results if res is not None])
        go_enrichment_results_df[f'adjusted_p_value ({experiment_name})'] = stats.false_discovery_control(
            ps=go_enrichment_results_df[f'p_value ({experiment_name})'],
            method='bh')
        go_enrichment_results_df[f'neg_log10_adjusted_p_value ({experiment_name})'] = np.round(
            -np.log10(go_enrichment_results_df[f'adjusted_p_value ({experiment_name})']),
            4)
        go_enrichment_results_df = go_enrichment_results_df.set_index(['GO_ID', 'Name', 'Term_Size'])
        results.append(
            go_enrichment_results_df.drop(
                [f'p_value ({experiment_name})', f'adjusted_p_value ({experiment_name})'],
                axis=1))
    result_df = pd.concat(results, axis=1)

    output_json = filepath.parent / f'go_enrichment_result.json'
    result_df.reset_index().to_json(path_or_buf=output_json,
                                    orient='records')

    return output_json
