"""This script retrieves the latest enzyme-substrate database from omnipath
and converts it into an adjacency matrix that can directly be used by kinact.
It is not run automatically because it takes quite long and database is relatively constant.
The date of the last execution is stored above the header of the data frame
"""

from datetime import date
from pypath import omnipath
import pandas as pd
import numpy as np

"""This is taken from the pypath-omnipath manual 
(https://github.com/saezlab/pypath/blob/master/docs/source/notebooks/manual.ipynb) """
es = omnipath.db.get_db('enz_sub')
es.make_df()
""" This is adapted from the get_kinase_targets() function in kinact """
es.df['p_site'] = es.df['substrate'].astype(str) + '_' + es.df['residue_type'].astype(str) + es.df[
    'residue_offset'].astype(str)
# Restrict the data on phosphorylation and dephosphorylation modification
ptms_omnipath_phospho = es.df.where(es.df['modification'].str.contains('phosphorylation')).dropna()
# Add value column
ptms_omnipath_phospho['value'] = 1
# all_sources = set(source for sources_str in ptms_omnipath_phospho['sources'] for source in sources_str.split(';'))
sources = ['PhosphoSite_MIMP', 'PhosphoSite_ProtMapper']
ptms_omnipath_phospho = ptms_omnipath_phospho[[len(set(s.split(';')).intersection(sources)) >= 1
                                               for s in ptms_omnipath_phospho['sources']]]
# Change value to negative value for dephosphorylation events
# This doesn't do anything right now, because the selected databases don't contain any dephosphorylation events.
ptms_omnipath_phospho.loc[ptms_omnipath_phospho['residue_type'].str.startswith('de'), 'value'] = -1
# Create pivot-table from pandas function
adjacency_matrix = pd.pivot_table(ptms_omnipath_phospho, values='value', index='p_site', columns='enzyme')
id_mapping = pd.read_csv('../id_conversion.txt')
fr = 'uniprot'
to = 'gene_name'
mapped_column_names = [id_mapping.loc[np.where(id_mapping[fr] == s)[0], to].values.tolist()[0]
                       if len(id_mapping.loc[np.where(id_mapping[fr] == s)[0], to].values.tolist()) == 1 else np.nan
                       for s in adjacency_matrix.columns]
adjacency_matrix.columns = mapped_column_names
with open('../psp_kinase_substrate_adjacency_matrix.csv', 'w') as f:
    f.write(f'### Database Version: {date.today()}\n')
    adjacency_matrix.to_csv(f)
