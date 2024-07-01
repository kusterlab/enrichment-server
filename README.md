# Enrichment Server

Developed and Maintained by Julian Müller (julian2.mueller@tum.de).

## Usage

(If you want to host an instance of the Enrichment Server yourself, please read the section **Hosting** below.)  
The Enrichment Server is currently running at this address: https://enrichment.kusterlab.org/main_enrichment-server/.
The currently implemented services are described below. You can use each one of them by sending a POST request
and attaching your input data in JSON format, as well as a session ID and a dataset name
(those are needed for PTMNavigator, you can use whatever - maybe I will implement defaults for that at some point).  
<b>Pro Tip:</b> If you are preparing your input data as a `pandas` data frame, an easy way to convert it into the
required input format
is using
`df.to_json(orient='records')`.

<details>  
<summary> <b>PTM Signature Enrichment Analysis</b>
</summary>

<i>Description</i>

PTM-Centric Enrichment Analysis using the PTM Signature Database (PTMSigDB).
Basically a GSEA that is Single-Site-Centric (ssc).

<i>Endpoint</i>

`/ssgsea/ssc`

<i>Reference</i>

Code: https://github.com/broadinstitute/ssGSEA2.0  
Publication: https://www.mcponline.org/article/S1535-9476(20)31860-0/fulltext

<i>Input</i>

1. `.../ssc/flanking`: A list of PTM sites surrounded by their +-7 flanking sequence, and their expression in each
   experiment.
   E.g.:

```
 [...,
 {
  "id":"ALLQLDGTPRVCRAA-p",
  "Experiment01": 15.7046003342,
  "Experiment02": 12.9784002304
 },
 ...]
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

<i>Example Command</i>

`curl -X POST -F file=@fixtures/ptm-sea/input/input_flanking.json
-F session_id=ABCDEF12345
-F dataset_name=ptm-sea http://10.152.171.101:4321/ssgsea/ssc/flanking
-o output_ptmsea_flanking.json`

`curl -X POST -F file=@fixtures/ptm-sea/input/input_uniprot.json
-F session_id=ABCDEF12345
-F dataset_name=ptm-sea http://10.152.171.101:4321/ssgsea/ssc/uniprot
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

`/ssgsea/gc`

<i>Reference</i>

Code: https://github.com/broadinstitute/ssGSEA2.0  
Publication: https://www.mcponline.org/article/S1535-9476(20)31860-0/fulltext

<i>Input</i>

A list of gene symbols, and their expression in each experiment.

E.g.:

```
 [...,
 {
  "id":"PSEN1",
  "Experiment01":10.0033998489,
  "Experiment02":14.6499004364
 },
 ...]
```

<i>Example Command</i>

`curl -X POST -F file=@fixtures/ssgsea/input/input.json
-F session_id=ABCDEF12345
-F dataset_name=genecentric http://10.152.171.101:4321/ssgsea/gc
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

`/ssgsea/gcr`

<i>Reference</i>

Code: https://github.com/broadinstitute/ssGSEA2.0  
Publication: https://www.mcponline.org/article/S1535-9476(20)31860-0/fulltext

<i>Input</i>

Identical to Non-Redundant Gene-Centric PEA.  
E.g.:

```
 [...,
 {
  "id":"PSEN1",
  "Experiment01":10.0033998489,
  "Experiment02":14.6499004364
 },
 ...]
```

<i>Example Command</i>

`curl -X POST -F file=@fixtures/ssgsea/input/input.json
-F session_id=ABCDEF12345
-F dataset_name=genecentricredundant http://10.152.171.101:4321/ssgsea/gcr
-o output_gcr.json`
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

`/ksea`

<i>Reference</i>

Code:  https://github.com/saezlab/kinact  
Publication:  https://www.science.org/doi/10.1126/scisignal.2003573

<i>Input</i>

E.g.:
A list of phosphosites, encoded in the format `<Uniprot_Acc>_<Res><Position>`, and their expression in each experiment.

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

<i>Example Command</i>

`curl -X POST -F file=@fixtures/ksea/input/input.json
-F session_id=ABCDEF12345
-F dataset_name=ksea http://10.152.171.101:4321/ksea
-o output_ksea.json`

</details>


<details>  
<summary> <b>KSEA with RoKAI</b>
</summary>

<i>Description</i>  
This endpoint uses `RoKAI` to refine the phosphorylation profiles before using `kinact` to perform KSEA.
`RoKAI` has been shown to produce more robust results when combined with any kinase activity inference method (see the
publication by Yılmaz et al. below).
We use all 5 components of RoKAI's functional/structural neighbourhood network as information source (see Fig. 3 in the
publication).

<i>Endpoint</i>

`/ksea/rokai`

<i>Reference</i>

Code: https://github.com/serhan-yilmaz/RokaiApp  
Publication: https://www.nature.com/articles/s41467-021-21211-6

<i>Input</i>

Identical to KSEA.  
E.g.:

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

<i>Example Command</i>

`curl -X POST -F file=@fixtures/ksea/input/input.json
-F session_id=ABCDEF12345
-F dataset_name=ksea_rokai http://10.152.171.101:4321/ksea/rokai
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

`/phonemes`

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

<i>Example Command</i>

`curl -X POST -F file=@fixtures/phonemes/input/input.json
-F session_id=ABCDEF12345
-F dataset_name=phonemes http://10.152.171.101:4321/phonemes
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

`/motif_enrichment`

<i>Reference</i>

Code: https://kinase-library.phosphosite.org  
Publication: https://www.nature.com/articles/s41586-022-05575-3

<i>Input</i>

A list of modified sequences, the Uniprot accession number(s) of the proteins they reside on,
and for each experiment whether the peptide was up- or down-regulated.
E.g.:

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

<i>Example Command</i>

`curl -X POST -F file=@fixtures/motif_enrichment/input/input.json
-F session_id=ABCDEF12345
-F dataset_name=motif_enrichment http://10.152.171.101:4321/motif_enrichment
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

`/kea3`

<i>Reference</i>

Code: https://maayanlab.cloud/kea3/templates/api.jsp  
Publication: https://academic.oup.com/nar/article/49/W1/W304/6279841

<i>Input</i>

A list of proteins for each experiment.
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
```

<i>Example Command</i>

`curl -X POST -F file=@fixtures/kea3/input/input.json
-F session_id=ABCDEF12345
-F dataset_name=kea3 http://10.152.171.101:4321/kea3
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

`/kstar`

<i>Reference</i>

Code: https://github.com/NaegleLab/KSTAR    
Publication: https://www.nature.com/articles/s41467-022-32017-5  

<i>Input</i>

A list of modified sequences, the Uniprot accession number(s) of the proteins they reside on,
and for each experiment the expression value of the peptide.
E.g.:

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

<i>Example Command</i>

`curl -X POST -F file=@fixtures/kstar/input/input.json
-F session_id=ABCDEF12345
-F dataset_name=kstar http://10.152.171.101:4321/kstar
-o output_kstar.json`

</details>

## Hosting
To run the Enrichment Server yourself, you need to obtain a few files that are not in this repository (either because they are too large for Github or require a license):

- Get a licensed copy of IBM CPLEX, which is required to run PHONEMeS (free for academic use): https://www.ibm.com/products/ilog-cplex-optimization-studio. 
Create a folder `CPLEX` in the same folder as this README file and place the executable file `cplex` in there. 
- To run KSTAR, obtain the NetworKIN pre-pruned network here: https://figshare.com/articles/dataset/NETWORKS/14944305. Create a folder `db/kstar` and place the `NetworKIN` folder there.

If you don't need the PHONEMeS and/or KSTAR endpoints, you can also skip these steps.  
Now, you can simply build and run the docker container. 