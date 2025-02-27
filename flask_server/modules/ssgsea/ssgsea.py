import subprocess
from pathlib import Path
import json
import pandas as pd
from cmapPy.pandasGEXpress import parse_gct

R_SPECIAL_CHARACTERS_MAPPING = str.maketrans({elem: '.' for elem in [
    " ",  # Space
    "+",  # Plus
    "-",  # Minus
    "*",  # Multiplication
    "/",  # Division
    "^",  # Exponentiation
    "#",  # Pound sign
    "%",  # Percentage
    "&",  # Ampersand
    "!",  # Exclamation mark
    "?",  # Question mark
    ",",  # Comma
    ".",  # Period (when at the start of a name)
    ";",  # Semicolon
    ":",  # Colon
    "@",  # At sign
    "~",  # Tilde
    "\\",  # Backslash
    '"',  # Double quotation mark
    "'",  # Single quotation mark
    "(",  # Left parenthesis
    ")",  # Right parenthesis
    "<",  # Left angle bracket
    ">",  # Right angle bracket
    "{",  # Left curly bracket
    "}",  # Right curly bracket
    "[",  # Left square bracket
    "]",  # Right square bracket
    "$",  # Dollar sign
    "=",  # Equal sign
    "|",  # Pipe
    "`"  # Backtick
]})


def preprocess_ssgsea(filepath: Path, type_isnot_gcr) -> Path:
    output_dir = filepath.parent
    input_json = json.load(open(filepath))
    input_df = pd.DataFrame.from_dict(input_json)

    idcolumn = 'Site' if 'Site' in input_df else 'id'

    # If it's a non-redundant gene-centric ssGSEA, we need to eliminate duplicates
    if type_isnot_gcr:
        def abs_max_signed(group):
            idx = group.abs().idxmax()
            return group.loc[idx] if pd.notna(idx) else float('nan')

        experiment_columns = [col for col in input_df.columns if col != idcolumn]
        input_df_grouped = input_df.groupby(idcolumn)
        unique_values = [input_df_grouped[exp].apply(abs_max_signed) for exp in experiment_columns]
        input_df = pd.DataFrame(unique_values).T.reset_index()

    # There is a method cmapPy.pandasGEXpress.write_gct,
    # but I could not get it to run
    # (it tries to use DataFrame.at[] with ranges and I couldn't find a pandas version where this works)
    # So I wrote a minimal version myself
    input_gct_file = output_dir / 'ssgsea_input.gct'
    with open(input_gct_file, 'w') as f:
        # Write Version and dims
        f.write("#1.3\n")
        f.write(f"{input_df.shape[0]}\t{input_df.shape[1] - 1}\t0\t0\n")
        # Write data df
        input_df.to_csv(f, sep='\t', index=False)
    return input_gct_file


def run_ssgsea(filepath: Path, ssgsea_type, ssc_input_type) -> Path:
    output_dir = filepath.parent
    output_prefix = output_dir / f'ssgsea_{ssgsea_type}_out'

    match ssgsea_type:
        case 'ssc':
            if ssc_input_type == 'flanking':
                database = "../ssGSEA2.0/db/ptmsigdb/ptm.sig.db.all.flanking.human.v2.0.0.gmt"
            elif ssc_input_type == 'uniprot':
                database = "../ssGSEA2.0/db/ptmsigdb/ptm.sig.db.all.uniprot.human.v2.0.0.gmt"
        case 'gc' | 'gcr':
            database = "../db/c2.cp.kegg+wp.v2023.2.Hs.symbols.gmt"

    subprocess_output = subprocess.run(["Rscript",
                             "../ssGSEA2.0/ssgsea-cli.R",
                             "-i", str(Path('..') / 'flask_server' / filepath),
                             "-o", str(output_prefix),
                             "-d", database,
                             "-w", "0.75",
                             "-e", "FALSE",
                             ],
                            capture_output=True, text=True)
    print(subprocess_output.stdout)
    print(subprocess_output.stderr)
    return Path(str(output_prefix) + '-combined.gct')


def postprocess_ssgsea(output_gct: Path) -> Path:
    output_json = output_gct.parent / f'{output_gct.stem}_result.json'
    if not output_gct.exists():
        with open(output_json, 'w') as o:
            o.write('[]')
    else:
        gct_parsed = parse_gct.parse(output_gct)
        experiment_names = gct_parsed.data_df.columns
        # Sanitize experiment names - R turns several characters into dots
        experiment_names = [name.translate(R_SPECIAL_CHARACTERS_MAPPING) for name in experiment_names]
        gct_df_joined = gct_parsed.row_metadata_df[[f'Signature.set.overlap.percent.{exp}' for exp in experiment_names]
                                                   + [f'fdr.pvalue.{exp}' for exp in experiment_names]
                                                   + [f'Signature.set.overlap.{exp}' for exp in experiment_names]
                                                   ].join(gct_parsed.data_df).reset_index()

        gct_df_joined.columns = (['Signature ID']
                                 + [f'Percent Overlap ({exp})' for exp in experiment_names]
                                 + [f'adj p-val ({exp})' for exp in experiment_names]
                                 + [f'Overlap ({exp})' for exp in experiment_names]
                                 + [f'Score ({exp})' for exp in experiment_names])

        gct_df_joined.to_json(path_or_buf=output_json, orient='records',
                              # indent=1  # For DEBUG
                              )
    return output_json
