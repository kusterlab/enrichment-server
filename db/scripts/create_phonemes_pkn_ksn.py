phonemes_pkn = pd.read_csv('../phonemesPKN.csv', index_col=0)
phonemes_ksn = pd.read_csv('../phonemesKSN.csv', index_col=0)
phonemes_pkn_ksn = pd.concat([phonemes_pkn, phonemes_ksn]).drop_duplicates()
phonemes_pkn_ksn.to_csv('../db/phonemes_PKN_KSN.csv')