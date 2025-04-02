# Enrichment Server

Developed and Maintained by Julian Müller (julian2.mueller@tum.de).

## Usage

The Enrichment Server is currently running here: https://enrichment.kusterlab.org/main_enrichment-server/  
The currently implemented services are described below.  
All services work with _Homo sapiens_ data. Most services additionally work with _Mus musculus_ data 
(this depends on the availability of an annotation database for the respective algorithm).
The organism can be specified using the parameter `organism` (use `'hsa'` for _H. sapiens_ and `'mmu'` for _M. musculus_; 
if not specified, `'hsa'` is assumed).   
Examples are given using `curl` (https://curl.se/), but you can use any software that can send a POST request,
e.g. `httr` if you're using R or `requests` for Python.  

### Input Format
All endpoints accept the input in JSON format and also return JSON data.  
Most endpoints additionally accept CSV input and output.

### Parameter File
For some endpoints, you can also supply additional parameters to customize the execution of the respective algorithm.
Please specify them in a file `parameters.toml` in the format `<key> = <value>` (each parameter on a separate line).  
If you do not supply a parameters file, default values will be used (all tunable parameters are described below;
for details, please always refer to the respective publication).


## Description of Endpoints
<details>  
<summary> <b>PTM Signature Enrichment Analysis</b>
</summary>

<i>Description</i>

PTM-Centric Enrichment Analysis using the PTM Signature Database (PTMSigDB).
Basically a GSEA that is Single-Site-Centric (ssc).

<i>Endpoint</i>

`/ssgsea/ssc`  (JSON + CSV)

<i>Reference</i>

Code: https://github.com/broadinstitute/ssGSEA2.0  
Publication: https://www.mcponline.org/article/S1535-9476(20)31860-0/fulltext

<i>Input</i>

1. `.../ssc/flanking`: A list of PTM sites surrounded by their +-7 flanking sequence, and their expression in each
   experiment.
   E.g.:

JSON:
```
 [...,
 {
  "id":"ALLQLDGTPRVCRAA-p",
  "Experiment01": 15.7046003342,
  "Experiment02": 12.9784002304
 },
 ...]
```

CSV:

```
id	Experiment01	Experiment02
ALLQLDGTPRVCRAA-p 15.7046003342  12.9784002304
...
```

2. `.../ssc/uniprot`: Alternatively, encode the sites as a list of Uniprot identifiers and site positions:
   E.g.:

```
 [...,
 {
  "id":"Q96MK2;T832-p",
  "Experiment01":15.7046003342,
  "Experiment02":12.9784002304
 },
 ...]
```

<i>Parameters</i>  
- `norm`: Type of Sample Normalization. Possible values are `'rank', 'log', 'log.rank', 'none'` (default: `'rank'`)
- `weight`: Weight of quantitative values (default: `0.75`)
- `correl`: Correlation Type. Possible values are `'rank', 'z.score', 'symm.rank'` (default: `'z.score'`)
- `test`: Test statistic. Possible values are `'area.under.RES', 'Kolmogorov-Smirnov'` (default: `'area.under.RES'`)
- `score`: Score type. Possible values are `'ES', 'NES'` (default: `'NES'`)
- `perm`: Number of permutations for p-value estimation (default: `1000`)
- `minoverlap`: Minimal overlap between signature and data set (default: `10`)

<i>Supported Organisms</i>  
`hsa`, `mmu`

<i>Example Command</i>

`curl -X POST -F file=@fixtures/ptm-sea/input/input_flanking_mmu.json
-F parameters=@fixtures/ptm-sea/input/parameters.toml
-F organism=mmu
https://enrichment.kusterlab.org/main_enrichment-server/ssgsea/ssc/flanking
-o output_ptmsea_flanking_mmu.json`

`curl -X POST -F file=@fixtures/ptm-sea/input/input_uniprot.json
-F organism=hsa
https://enrichment.kusterlab.org/main_enrichment-server/ssgsea/ssc/uniprot
-o output_ptmsea_uniprot.json`
</details>  



<details>  
<summary> <b>Gene-Centric Pathway Enrichment Analysis</b>
</summary>

<i>Description</i>  
Basically a GSEA against a database of pathway signatures.
We use the same algorithm as for PTM-SEA (ssGSEA),
but with the MSigDB database instead of PTMSigDB (https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp).
This means when using this endpoint on a PTM datasets, the site-specific information cannot be used (data has to be
collapsed to gene level).  
We use the KEGG and Wikipathways signatures only
(running against the entire MSigDB would take a long time and is strongly discouraged by the creators).

<i>Endpoint</i>

`/ssgsea/gc` (JSON + CSV)

<i>Reference</i>

Code: https://github.com/broadinstitute/ssGSEA2.0  
Publication: https://www.mcponline.org/article/S1535-9476(20)31860-0/fulltext

<i>Input</i>

A list of gene symbols, and their expression in each experiment.

E.g.:  
JSON:

```
 [...,
 {
  "id":"PSEN1",
  "Experiment01":10.0033998489,
  "Experiment02":14.6499004364
 },
 ...]
```

CSV:

```
id	Experiment01	Experiment02
PSEN1 10.0033998489  14.6499004364
...
```

<i>Parameters</i>  
Same as for PTM-SEA.   

<i>Supported Organisms</i>  
`hsa`, `mmu`

<i>Example Command</i>

`curl -X POST -F file=@fixtures/ssgsea/input/input.json
https://enrichment.kusterlab.org/main_enrichment-server/ssgsea/gc
-o output_gc.json`
</details>


<details>  
<summary> <b>Gene-Centric Redundant Pathway Enrichment Analysis</b>
</summary>

<i>Description</i>  
The only difference to gene-centric enrichment is that genes are repeatedly counted for each regulated site in the data.
It was shown in Krug et al. 2019 that while not performing as good as PTM-level enrichment,
this works better than only counting each gene with regulated sites once, regardless of the number of regulated sites.
Since gene-centric signatures are more comprehensive than site-centric signatures (e.g., they cover all human
WikiPathways and KEGG pathways),
it poses a good compromise between the two approaches.

<i>Endpoint</i>

`/ssgsea/gcr`  (JSON + CSV)

<i>Reference</i>

Code: https://github.com/broadinstitute/ssGSEA2.0  
Publication: https://www.mcponline.org/article/S1535-9476(20)31860-0/fulltext

<i>Input</i>

Identical to Non-Redundant Gene-Centric PEA.  
E.g.:

JSON:

```
 [...,
 {
  "id":"PSEN1",
  "Experiment01":10.0033998489,
  "Experiment02":14.6499004364
 },
 ...]
```

CSV:

```
id	Experiment01	Experiment02
PSEN1 10.0033998489  14.6499004364
...
```

<i>Parameters</i>  
Same as for PTM-SEA.   

<i>Supported Organisms</i>  
`hsa`, `mmu`


<i>Example Command</i>

`curl -X POST -F file=@fixtures/ssgsea/input/input.json
https://enrichment.kusterlab.org/main_enrichment-server/ssgsea/gcr
-o output_gcr.json`
</details>

<details>  
<summary> <b>GO Enrichment</b>
</summary>

<i>Description</i>  
This endpoints performs Gene Ontology (GO) Enrichment Analysis using Fisher's exact test.
It is a more straightforward method than ssGSEA, and requires only a set of differentially regulated genes as input.
Each gene set in the GO database is then tested for overrepresentation.  
In order to use this approach with PTM datasets, you need to map the peptide/site-level information to gene-level
information.
This approach also cannot take into account the direction, fold change, or significance of each regulation;
therefore we advise you to use PTM-SEA or ssGSEA instead, if you have more than just a list of regulated genes.  
As database, we use the GO annotations provided by g:Profiler.

<i>Endpoint</i>

`/go_enrichment` (JSON)

<i>Reference</i>

Code:  Custom, for running the Fisher test we use
SciPy: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html  
Publication:  https://www.nature.com/articles/ng0500_25, https://academic.oup.com/genetics/article/224/1/iyad031/7068118

<i>Input</i>

A list of differentially regulated genes for each experiment in your dataset.
You can supply a custom background (e.g. all genes measured in your model system) using the "Background" key (
case-insensitive).
E.g.:

```
{
  "Experiment01": [
    "FOXM1",
    "SMAD9"
  ],
  "Experiment02": [
    "ZNF264",
    "TMPO",
    "ISL2"
  ],
  "Background": [
    "SLC25A25",
    "TMEM217",
    ...
  ]
}
```

The default background are all genes annotated in the GO database.

<i>Supported Organisms</i>  
`hsa`, `mmu`


<i>Example Command</i>

`curl -X POST -F file=@fixtures/go_enrichment/input/input.json
https://enrichment.kusterlab.org/main_enrichment-server/go_enrichment
-o output_go.json`

</details>


<details>  
<summary> <b>KSEA</b>
</summary>

<i>Description</i>  
KSEA uses phosphoproteomics data (usually fold changes) and prior knowledge on kinase-substrate relationships to infer
kinase activities.
There are multiple implementations for KSEA, we use the one from the `kinact` package,
which compares the mean fold change among the set of substrates of a kinase to an expected value.
The implementation is based on a publication by Casado et al. (see below).
The prior knowledge we use are the most recent kinase-substrate relationships from PhosphoSitePlus, retrieved using
Omnipath on 2024-02-11. If you're interested, you can find the code to update the database
in `db/scripts/update_ksea_es_db.py`.

<i>Endpoint</i>

`/ksea`, `/kinact`   (JSON + CSV)

<i>Reference</i>

Code:  https://github.com/saezlab/kinact  
Publication:  https://www.science.org/doi/10.1126/scisignal.2003573

<i>Input</i>

A list of phosphosites, encoded in the format `<Uniprot_Acc>_<Res><Position>`, and their expression in each experiment.
E.g.:

JSON:

```
 [...,
 {
  "Site":"O75822_S11",
  "Experiment_1":0.0,
  "Experiment_2":-0.002266224,
  "Experiment_3":0.0
 },
 ...]
```

CSV:

```
Site,Experiment_1,Experiment_2,Experiment_3
O75822_S11,0.0,-0.002266224,0.0
...
```

<i>Parameters</i>  
- `minimum_set_size`: Minimum overlap between the sites in the dataset and in the substrate set of a kinase (default: `5`)
- `median`: Use median instead of mean substrate Fold Change. Possible values are `true, false` (default: `false`)


<i>Supported Organisms</i>  
`hsa`, `mmu`

<i>Example Command</i>

`curl -X POST -F file=@fixtures/ksea/input/input.json
https://enrichment.kusterlab.org/main_enrichment-server/ksea
-o output_ksea.json`

</details>


<details>  
<summary> <b>RoKAI</b>
</summary>

<i>Description</i>  
This endpoint uses `RoKAI` to refine the phosphorylation profiles and subsequently infer kinase activities.
`RoKAI` has been shown to produce more robust results when combined with any kinase activity inference method (see the
publication by Yılmaz et al. below).
We use all 5 components of RoKAI's functional/structural neighbourhood network as information source (see Fig. 3 in the
publication).

<i>Endpoint</i>

`/rokai`  (JSON + CSV)

<i>Reference</i>

Code: https://github.com/serhan-yilmaz/RokaiApp  
Publication: https://www.nature.com/articles/s41467-021-21211-6

<i>Input</i>

Identical to KSEA.  
E.g.:

JSON:

```
 [...,
 {
  "Site":"O75822_S11",
  "Experiment_1":0.0,
  "Experiment_2":-0.002266224,
  "Experiment_3":0.0
 },
 ...]
```

CSV:

```
Site,Experiment_1,Experiment_2,Experiment_3
O75822_S11,0.0,-0.002266224,0.0
...
```

<i>Parameters</i>  
- `datanorm`: Normalization Strategy. Possible values are `'Normalized', 'Centered', 'Raw'` (default: `'Normalized'`)
- `signor`: Whether or not to use Signor kinase-substrate annotation in addition to PhosphoSitePlus. Possible values are `true, false` (default: `false`)
- `ppi`: Whether or not to include the protein-protein interactions from STRING in the RoKAI network. Possible values are `true, false` (default: `true`)
- `sd`: Whether or not to include the structural distance information from PTMcode in the RoKAI network. Possible values are `true, false` (default: `true`)
- `coev`: Whether or not to include the coevolution information from PTMcode in the RoKAI network. Possible values are `true, false` (default: `true`)

<i>Supported Organisms</i>  
`hsa`, `mmu`


<i>Example Command</i>

`curl -X POST -F file=@fixtures/rokai/input/input.json
https://enrichment.kusterlab.org/main_enrichment-server/rokai
-F organism=mmu
-o output_rokai.json`

</details>



<details>  
<summary> <b>KSEA with RoKAI</b>
</summary>

<i>Description</i>  
This endpoint uses only the first part of `RoKAI` to refine the phosphorylation profiles.
These are then further processed by `kinact` to perform KSEA.

<i>Endpoint</i>

`/ksea/rokai`  (JSON + CSV)

<i>Input</i>

Identical to KSEA.  
E.g.:

JSON:

```
 [...,
 {
  "Site":"O75822_S11",
  "Experiment_1":0.0,
  "Experiment_2":-0.002266224,
  "Experiment_3":0.0
 },
 ...]
```

CSV:

```
Site,Experiment_1,Experiment_2,Experiment_3
O75822_S11,0.0,-0.002266224,0.0
...
```

<i>Parameters</i>  
All parameters from both KSEA and RoKAI.


<i>Supported Organisms</i>  
`hsa`, `mmu`


<i>Example Command</i>

`curl -X POST -F file=@fixtures/ksea/input/input.json
https://enrichment.kusterlab.org/main_enrichment-server/ksea/rokai
-o output_ksea_rokai.json`

</details>



<details>
<summary> <b>PHONEMeS</b>
</summary>

<i>Description</i>

`PHONEMeS` uses a prior knowledge network of PPIs and Kinase-Substrate Relationships to reconstruct
a signaling network from a phosphoproteomics dataset and a set of perturbation targets.
The current version is a wrapper around the causal reasoning tool `CARNIVAL`.
Essentially it works by trimming away parts of the prior knowledge network until the resulting subnetwork
optimally explains the observed data.    
This endpoint first runs PHONEMeS on the input data and uses Cytoscape to set 2-D coordinates for the protein nodes.  
The _yFiles_ plugin (https://www.yworks.com/products/yfiles-layout-algorithms-for-cytoscape) is utilized to arrange the
graph in a hierarchic layout. The result is converted into JSON format and sent back to the User.
Note that the phosphosite nodes are trimmed away from the PHONEMeS result, only protein
nodes are returned.

<i>Endpoint</i>

`/phonemes` (JSON)

<i>Reference</i>

Code: https://github.com/saezlab/PHONEMeS  
Publication: https://pubs.acs.org/doi/full/10.1021/acs.jproteome.0c00958

<i>Input</i>

A list of targets, split by experiment and regulation direction, as well as a list of sites,
encoded in the format `<Uniprot_Acc>_<Res><Position>`, together with the expression of each site in each experiment.

E.g.:

```
{
  "targets": {
    "Experiment01": {
      "up": [
        "RICTOR"
      ],
      "down": [
        "EGFR",
        "MAPKAPK2"
      ]
    },
    "Experiment02": {
      "up": [
        "AHNAK",
        "MTOR"
      ],
      "down": [
        "AKT1S1"
      ]
    }
  },
    "sites":  [...,
       {
        "Site":"O75822_S11",
        "Experiment_1":0.0,
        "Experiment_2":-0.002266224,
        "Experiment_3":0.0
       },
 ...]
 }
```

<i>Supported Organisms</i>  
`hsa`


<i>Example Command</i>

`curl -X POST -F file=@fixtures/phonemes/input/input.json
https://enrichment.kusterlab.org/main_enrichment-server/phonemes
-o output_phonemes.json`

</details>


<details>  
<summary> <b>Motif Enrichment</b>
</summary>

<i>Description</i>

Performs a Kinase Motif Enrichment by making use of the Kinase Library (Johnson et al., Nature 2023).  
Position-specific scoring matrices are used to score the motif of each kinase against a phosphoproteomics dataset.  
The endpoint returns the enrichment values for every scored kinase motif.

<i>Endpoint</i>

`/motif_enrichment`  (JSON + CSV)

<i>Reference</i>

Code: https://kinase-library.phosphosite.org  
Publication: https://www.nature.com/articles/s41586-022-05575-3

<i>Input</i>

A list of modified sequences, the Uniprot accession number(s) of the proteins they reside on,
and for each experiment whether the peptide was up- or down-regulated.
E.g.:

JSON:
```
 [...,
  {
    "Modified sequence": "RDS(ph)ASYR",
    "Proteins": "A0A1X7SBZ2;A0A5H1ZRQ2;Q92841;Q92841-1;Q92841-2;Q92841-3",
    "Experiment01": "down",
    "Experiment02": "up"
  },
 ...]
```

CSV:
```
Modified sequence,Proteins,Experiment01,Experiment02
RDS(ph)ASYR,A0A1X7SBZ2;A0A5H1ZRQ2;Q92841;Q92841-1;Q92841-2;Q92841-3,down,up
...
```

<i>Parameters</i>   
- `top_n`: How many top kinases to consider (default: `15`)
- `threshold`: Filters for values at least this big (default: `-inf` (specify a numeric value))
- `threshold_type`: Metric to use for filtering. Possible values are `'score', 'percentile', 'total'` (default: `'percentile'`)
- `sort_type`: Metric to use for ranking. Possible values are `'score', 'percentile', 'total'` (default: `'percentile'`)

<i>Supported Organisms</i>  
`hsa`


<i>Example Command</i>

`curl -X POST -F file=@fixtures/motif_enrichment/input/input.json
https://enrichment.kusterlab.org/main_enrichment-server/motif_enrichment
-o output_motif_enrichment.json`

</details>


<details>  
<summary> <b>KEA3</b>
</summary>

<i>Description</i>

Performs Kinase Enrichment Analysis 3 (KEA3) enrichment.
KEA3 infers upstream kinases whose putative substrates are overrepresented
in a user-inputted list of proteins or differentially phosphorylated proteins.  
The endpoint calls the API of KEA3 and returns the `MeanRank` and `TopRank` tables of the query result.

<i>Endpoint</i>

`/kea3`  (JSON)

<i>Reference</i>

Code: https://maayanlab.cloud/kea3/templates/api.jsp  
Publication: https://academic.oup.com/nar/article/49/W1/W304/6279841

<i>Input</i>

A list of gene names representing proteins that are differentially regulated in each experiment.
E.g.:

```
{
  "Experiment01": [
    "FOXM1",
    "SMAD9"
  ],
    "Experiment02": [
    "ZNF264",
    "TMPO",
    "ISL2"
  ]
}
```

<i>Supported Organisms</i>  
`hsa`


<i>Example Command</i>

`curl -X POST -F file=@fixtures/kea3/input/input.json
https://enrichment.kusterlab.org/main_enrichment-server/kea3
-o output_kea3.json`

</details>

<details>  
<summary> <b>KSTAR</b>
</summary>

<i>Description</i>

Performs Kinase Activity Prediction using the KSTAR algorithm.  
Since KSTAR can only test for activity changes in one direction at a time, we only score down-regulations.  
As a threshold for retaining phosphorylation sites, we use a fixed value of 0, i.e., we retain all negative values.
Thus, the user needs to make sure to filter out non-significant regulations before using the endpoint.  
For reasons of performance, this endpoint only performs the hypergeometric tests for calculating enrichment scores
and p-values. The subsequent random analysis and Mann-Whitney-U test steps are omitted since they require significantly
more processing power and time.  

<i>Endpoint</i>

`/kstar` (JSON + CSV)

<i>Reference</i>

Code: https://github.com/NaegleLab/KSTAR    
Publication: https://www.nature.com/articles/s41467-022-32017-5  

<i>Input</i>

A list of modified sequences, the Uniprot accession number(s) of the proteins they reside on,
and for each experiment the expression value of the peptide.
E.g.:

JSON:
```
 [...,
 {
  "Modified sequence":"RS(ph)VGSDE",
  "Proteins":"C9JBX5;E9PAL7;P43307;P43307-2",
  "Experiment01":-1.2895137775,
  "Experiment02":-2.2462854621
 },
 ...]
```

CSV:
```
Modified sequence	Proteins	Experiment01	Experiment02
RS(ph)VGSDE C9JBX5;E9PAL7;P43307;P43307-2 -1.2895137775  -2.2462854621
...
```

<i>Parameters</i>   
- `agg`: How to aggregate sites that appear multiple times. Possible values are `'mean', 'max', 'min', 'count'` (default: `'mean'`)
- `threshold`: Cutoff for keeping a site as evidence  (default: `0`)

<i>Supported Organisms</i>  
`hsa`


<i>Example Command</i>

`curl -X POST -F file=@fixtures/kstar/input/input.json
https://enrichment.kusterlab.org/main_enrichment-server/kstar
-o output_kstar.json`

</details>

## Hosting
If you would like to host an instance of the Enrichment Server yourself, there are two preliminary steps: 

- Download **CPLEX** _(required for the PHONEMeS endpoint)_: Obtain a license for IBM CPLEX (free for academic use) 
and download it from here: https://www.ibm.com/products/ilog-cplex-optimization-studio.
Then, create a folder `CPLEX` at the top level of this repository and put the `cplex` executable file inside it.
- Download pre-pruned **NetworKIN** networks _(required for the KSTAR endpoint)_: Download the networks from here: 
https://figshare.com/articles/dataset/NETWORKS/14944305?file=28768155 and unzip it.
Then, create a folder `kstar` inside the `db` subfolder of this repository and copy the `NetworKIN` folder there.

If you don't want to make use of the PHONEMeS or KSTAR endpoint(s), you can also skip these steps.

Now you can just build and run the docker container:  
`docker build -t enrichment_server .`  
`docker run --network host enrichment_server`

## Reference
If the Enrichment Server is useful for your research, please cite the following publication:  

 **PTMNavigator: Interactive Visualization of Differentially Regulated Post-Translational Modifications in Cellular Signaling Pathways**  
Julian Müller, Florian P. Bayer, Mathias Wilhelm, Maximilian G. Schuh, Bernhard Kuster, Matthew The  
_Nature Communications_ 16:510 (2025); doi: https://doi.org/10.1038/s41467-024-55533-y