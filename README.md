# Enrichment Server

Developed and Maintained by Julian Müller (julian2.mueller@tum.de).

## Usage

The Enrichment Server is currently running internally on Atlas at Port 4321 (http://10.152.171.101:4321).
The currently implemented services are described below. You can use each one of them by sending a POST request
and attaching your input data in JSON format, as well as a session ID and a dataset name
(those are needed for PTMNavigator, you can use whatever - maybe I will implement defaults for that at some point).  
<b>Pro Tip:</b> If you are preparing your input data as a `pandas` data frame, an easy way to convert it into the required input format
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

1. `.../ssc/flanking`: A list of PTM sites surrounded by their +-7 flanking sequence, and their expression in each experiment.
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
This means when using this endpoint on a PTM datasets, the site-specific information cannot be used (data has to be collapsed to gene level).  
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
Since gene-centric signatures are more comprehensive than site-centric signatures (e.g., they cover all human WikiPathways and KEGG pathways), 
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
KSEA uses phosphoproteomics data (usually fold changes) and prior knowledge on kinase-substrate relationships to infer kinase activities.
There are multiple implementations for KSEA, we use the one from the `kinact` package, 
which compares the mean fold change among the set of substrates of a kinase to an expected value. 
The implementation is based on a publication by Casado et al. (see below). 
The prior knowledge we use are the most recent kinase-substrate relationships from PhosphoSitePlus, retrieved using Omnipath on 2024-02-11. If you're interested, you can find the code to update the database in `db/scripts/update_ksea_es_db.py`.    

<i>Endpoint</i>

`/ksea`

<i>Reference</i>

Code:  https://github.com/saezlab/kinact  
Publication:  https://europepmc.org/article/MED/23532336



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
`RoKAI` has been shown to produce more robust results when combined with any kinase activity inference method (see the publication by Yılmaz et al. below).
We use all 5 components of RoKAI's functional/structural neighbourhood network as information source (see Fig. 3 in the publication). 

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