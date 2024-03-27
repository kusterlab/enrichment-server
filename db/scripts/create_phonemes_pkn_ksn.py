import pandas as pd
phonemes_pkn = pd.read_csv('../phonemesPKN.csv')
phonemes_ksn = pd.read_csv('../phonemesKSN.csv')
phonemes_pkn_ksn = pd.concat([phonemes_pkn, phonemes_ksn]).drop_duplicates()
phonemes_pkn_ksn.to_csv('../phonemes_PKN_KSN.csv', index=False)