# TODO: Give FB Credit for the logic

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import psite_annotation as pa
from scipy import interpolate
from scipy.stats import fisher_exact
import statsmodels.api as sm

PHOSPHOSITE_FASTA = '../db/Phosphosite_seq.fasta'
ODDS_PATH = '../db/kinase_library/Motif_Odds_Ratios.txt'
QUANTILE_MATRIX_PATH = '../db/kinase_library/Kinase_Score_Quantile_Matrix.txt'


def run_motif_enrichment(filepath: Path) -> Path:
    ## Load the ODD ratios
    ODDS = pd.read_csv(ODDS_PATH, sep='\t', index_col=['Kinase', 'Position', 'AA'])

    ODDS = ODDS['Odds Ratio'].to_dict()

    ## Load the qunatiles
    QUANTILE_MATRIX = pd.read_csv(QUANTILE_MATRIX_PATH, sep='\t', index_col='Score').T
    QUANTILES = {}
    for kinase, q in QUANTILE_MATRIX.iterrows():
        QUANTILES[kinase] = build_ecdf(scores=q.index, quantiles=q.values)

    input_json = json.load(open(filepath))
    input_df = pd.DataFrame.from_dict(input_json)
    experiment_columns = [col for col in input_df.columns if col not in ['Modified sequence', 'Proteins']]

    # TODO: Maybe add more lowercase letters in context for multiphospho-peptides
    input_df = pa.addPeptideAndPsitePositions(input_df, PHOSPHOSITE_FASTA, pspInput=True, context_left=5,
                                              context_right=5)

    # Explode for multiple phosphos becomming individual rows
    input_df['Site sequence context'] = input_df['Site sequence context'].str.split(';')
    input_df['Site weight'] = 1 / input_df['Site sequence context'].apply(len)
    input_df = input_df.explode('Site sequence context').reset_index(drop=True)
    input_df['Site weight'] = input_df['Site weight'] / input_df.groupby('Site sequence context')[
        'Site weight'].transform('size')

    # Filter missed matches
    mask = (input_df['Site sequence context'] == '')
    input_df = input_df[~mask]

    enrichment_dfs = []
    for experiment in experiment_columns:
        input_experiment = input_df[['Site sequence context', 'Site weight', experiment]].rename(
            {experiment: 'Regulation'}, axis=1).dropna().copy()

        enrichment_df_experiment = motif_enrichment_analysis(input_experiment, P=ODDS, Q=QUANTILES, site_weights=True)
        enrichment_df_experiment = correct_for_multipletesting(enrichment_df_experiment)
        enrichment_df_experiment.columns = [f'{col} ({experiment})' for col in enrichment_df_experiment.columns]
        enrichment_dfs.append(enrichment_df_experiment)

    output_json = filepath.parent / f'motif_enrichment_result.json'
    pd.concat(enrichment_dfs, axis=1).reset_index(names='Kinase').to_json(path_or_buf=output_json,
                                                                          orient='records',
                                                                          # indent=1
                                                                          )
    return output_json


def build_ecdf(scores, quantiles):
    """
    Build a empirical cumulative distribution function (ecdf) based on precomputed scores and quantiles

    Parameter
    ---------
    scores : array like
        input scores
    quantiles : array like
        input qunatiles

    Return
    ------
    ecdf function that maps x->cumprop
    """
    # Return the ecdf as function
    f = interpolate.interp1d(scores, quantiles, kind='linear')
    return f


def motif_enrichment_analysis(df, P, Q, top_n=15, threshold=-np.inf, threshold_type='percentile',
                              sort_type='percentile', site_weights=False):
    """
    Motif enrichment according to the Johnson paper DOI: 10.1038/s41586-022-05575-3 as default values.
    More options by selecting different sortings, thresholds, site_weights for counting, etc..
    This is not identical to Phosphosite plus !!!

    Input
    -----
    df: DataFrame with 'Regulation' & 'SITE_+/-5_AA' column
    P : ODDS Matrix
    Q : QUANTILE Matrix

    Returns
    -------
    enrichment DataFrame with results
    """
    # Filter for correct regulation
    regulation_types = {'down', 'up', 'not'}
    df = df[df['Regulation'].isin(regulation_types)]

    ## Annotate the Sites with the best kinases
    cols = ['Top Motif Kinases', 'Top Motif Scores', 'Top Motif Percentiles', 'Top Motif Totals']

    motif = df['Site sequence context'].apply(find_upstream_kinase, Q=Q, P=P, top_n=top_n, threshold=threshold,
                                              threshold_type=threshold_type, sort_type=sort_type)
    motif = pd.DataFrame(motif.tolist(), index=motif.index, columns=cols)
    motif['Regulation'] = df['Regulation'].copy()

    # initialize the 'Site weight' column and if customized overwrite from input df
    if (not site_weights) or ('Site weight' not in df):
        motif['Site weight'] = 1
    else:
        motif['Site weight'] = df['Site weight'].copy()

    ## Count total down, not, up sites
    total_regulations = motif.groupby('Regulation')['Site weight'].sum().to_dict()
    total_regulations = pd.Series({reg: total_regulations.get(reg, 0) for reg in regulation_types})

    ## For each site, count if a Kinase was present with the correct regulation category
    counter = defaultdict(lambda: {'down': 0, 'up': 0, 'not': 0})
    for i, row in motif.iterrows():
        regulation = row['Regulation']
        site_weight = row['Site weight']
        for kinase in row['Top Motif Kinases'].split(';'):
            counter[kinase][regulation] += site_weight
    kinase_counts = pd.DataFrame(counter).T

    ## Make the fisher exact test for each kinase
    enrichment = {}
    for kinase, counts in kinase_counts.iterrows():
        enrichment[kinase] = kinase_motif_enrichment_test(counts, total_regulations)
    enrichment = pd.DataFrame(enrichment, index=['-Log10 p_value', 'Log2 Enrichment']).T
    enrichment = pd.merge(kinase_counts, enrichment, left_index=True, right_index=True)
    return enrichment


def find_upstream_kinase(seq, Q, P, top_n=15, threshold=-np.inf, threshold_type='percentile', sort_type='percentile'):
    """
    Score all kinases against input sequence based on Q-Matrix and P-Matrix.
    Percentile is the standard metic according to Johnson et al. and has the best perfromace in my hands as well.

    Input
    -----
    seq : str
        seqeunce to score
    Q : dict
        Quantile matrix Q
    P : dict
        Probabilty matrix P
    top_n : int
        considers the top n ranked kinases only, default 15.
    threshold: float
        filters for total value bigger than threshold, default -np.inf.
    threshold_type: float
        specifies which metric should be used for filtering (0=score, 1=percentile, 2=score*percentile), default is percentile [1]
    sort_type: float
        specifies which metric should be used for ranking (0=score, 1=percentile, 2=score*percentile), default is percentile [1]

    Returns
    -------
    (kinases, scores, percentiles, totals)
    each is a string with semicolon sorted values
    """
    # Map the different parameter options
    str_to_int_map = {'score': 0, 'percentile': 1, 'total': 2, }
    if threshold_type not in str_to_int_map:
        raise ValueError('threshold_type')
    if sort_type not in str_to_int_map:
        raise ValueError('sort_type')
    threshold_type = str_to_int_map[threshold_type]
    sort_type = str_to_int_map[sort_type]

    # Fast return for pY
    if seq[len(seq) // 2] == 'y':
        return ('', '', '', '')
    # score all kinases
    result = {}
    for kinase in Q.keys():
        s = score(seq, kinase, P, motif_size=5, sig_digits=3)
        if not np.isfinite(s):
            continue
        q = round(float(Q[kinase](s)), 3)
        t = round(s * q, 3)
        result[kinase] = (s, q, t)
    # Sort by sort_type, then filter threshold_type > threshold, and then take the topN of the list
    out = sorted(result.items(), key=lambda item: item[1][sort_type], reverse=True)
    out = [(k, tpl) for k, tpl in out if tpl[threshold_type] > threshold]
    out = out[:top_n]
    # Transofrom to ; seperated lists
    kinases = ';'.join(k for k, tpl in out)  # kinase
    scores = ';'.join(str(tpl[0]) for k, tpl in out)  # [0] score
    percentiles = ';'.join(str(tpl[1]) for k, tpl in out)  # [1] percentile
    totals = ';'.join(str(tpl[2]) for k, tpl in out)  # [2] total = s*q
    return (kinases, scores, percentiles, totals)


def score(seq, kinase, P, motif_size=5, sig_digits=4):
    """
    Score the motife based on the AA positional ODDS matrix given a sequence and a kinase
    """
    assert len(seq) == (2 * motif_size + 1)
    score = []
    for i, aa in enumerate(seq):
        pos = i - motif_size
        score.append(P.get((kinase, pos, aa), np.nan))
    return round(np.log2(np.nanprod(score)), sig_digits)


def kinase_motif_enrichment_test(counts, total_regulations):
    """
    Motif enrichment according to the Johnson paper DOI: 10.1038/s41586-022-05575-3
    This is not identical to Phosphosite plus !!!

    counts : pd.Series with counts for down up not for one kinase {down:N, up:N, not:N}
    total_regulations : pd.Series with counts for down up not for all sites {down:Total_N, up:Total_N, not:Total_N}
    """
    # Build contingency table
    table_down = construct_contingency_tabel('down', counts, total_regulations)
    table_up = construct_contingency_tabel('up', counts, total_regulations)

    # Fisher exact test
    result_down = fisher_exact(table_down, alternative='greater')
    result_up = fisher_exact(table_up, alternative='greater')

    # up direction wins
    if result_down.statistic < result_up.statistic:
        return -np.log10(result_up.pvalue), calculate_log_odds_ratio(table_up)
    # down direction wins
    if result_down.statistic > result_up.statistic:
        return -np.log10(result_down.pvalue), -1 * calculate_log_odds_ratio(table_down)
    # No one wins (normally both p-value == 1.0)
    return -np.log10((result_down.pvalue + result_up.pvalue) / 2), 0.0


def construct_contingency_tabel(cat, counts, totals):
    """
    Makes a 2x2 contingency tabel

                     Regulation:
                   True  |  False
      Kinase: ||=========|=========||
        True  ||    a    |    b    ||
        ------||---------|---------||
        False ||    c    |    d    ||
              ||=========|=========||
    """
    # Calculate the elements
    mask = (counts.index == cat)
    a = counts[mask].sum()
    b = counts[~mask].sum()
    mask = (totals.index == cat)
    c = totals[mask].sum() - a
    d = totals[~mask].sum() - b
    # construct the table
    table = np.array([[a, b], [c, d]])
    table = haldane_correction(table)
    return table


def haldane_correction(a):
    # Add 0 every where if a 0 is present in the contingency tabale a
    if (a == 0).any():
        return a + 1
    return a


def calculate_log_odds_ratio(table):
    # Calculate the log2 odds ratio based on 2x2 contingency tabel
    odds_ratio = (table[(0, 0)] * table[(1, 1)]) / (table[(0, 1)] * table[(1, 0)])
    return np.log2(odds_ratio)


def correct_for_multipletesting(df):
    p_vals = 10 ** (-1 * df['-Log10 p_value'])
    reject, p_vals_corrected, _, _ = sm.stats.multipletests(p_vals, alpha=0.05, method='fdr_bh')
    df['-Log10 p_value adjusted'] = - np.log10(p_vals_corrected)
    return df
