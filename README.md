# GMMAD2_data
Parser for GMMDA2 database (microbes-metabolites, microbes-diseases, and metabolites-diseases)
***

## Microbe-Metabolite
Record example: <br>
```ruby
{
   "_id":"264203_associated_with_6287",
   "association":{
      "predicate":"biolink:associated_with",
      "infores":"wom",
      "sources":[
         "host",
         "microbiota",
         "food related",
         "drug related",
         "environment"
      ]
   },
   "object":{
      "id":"PUBCHEM.COMPOUND:6287",
      "name":"valine",
      "type":"biolink:SmallMolecule",
      "pubchem_cid":6287,
      "chemical_formula":"C5H11NO2",
      "smiles":"CC(C)C(C(=O)O)N",
      "xrefs":{
         "kegg":"C00183",
         "hmdb":"HMDB0000883"
      }
   },
   "subject":{
      "id":"taxid:264203",
      "name":"zymomonas mobilis strain zm4 (atcc 31821)",
      "type":"biolink:Bacterium",
      "taxid":264203,
      "scientific_name":"zymomonas mobilis subsp. mobilis zm4 = atcc 31821",
      "parent_taxid":120045,
      "lineage":[
         264203,
         120045,
         542,
         541,
         2844881,
         204457,
         28211,
         1224,
         2,
         131567,
         1
      ],
      "rank":"strain"
   }
}
```

## Microbe-Disease
Total record before excluding sample size is 0: 508141
<br>
Total record after excluding sample size is 0: 121285
<br>
Edge reduced by 386,856 <br>
Organism Type Count: <br>
Bacterium: 119590 <br>
OrganismTaxon: 830 <br>
Archaea: 865 <br>

Record example: <br>
```ruby
{
   "_id":"1755642_OrganismalEntityAsAModelOfDiseaseAssociation_D003093",
   "association":{
      "predicate":"OrganismalEntityAsAModelOfDiseaseAssociation",
      "control_name":"healthy control",
      "qualifier":"decrease",
      "qualifier_ratio":"-0.001222806",
      "disease_sample_size":"3",
      "disease_abundance_mean":"0.000544065",
      "disease_abundance_median":"0.000506614",
      "disease_abundance_sd":"0.000351942",
      "healthy_sample_size":"233",
      "healthy_abundance_mean":"0.002674384",
      "healthy_abundance_median":"0.00172942",
      "healthy_abundance_sd":"0.002712906"
   },
   "object":{
      "id":"MESH:D003093",
      "name":"colitis, ulcerative",
      "mesh":"D003093",
      "type":"biolink:Disease",
      "description":"Inflammation of the COLON that is predominantly confined to the MUCOSA. Its major symptoms include DIARRHEA, rectal BLEEDING, the passage of MUCUS, and ABDOMINAL PAIN."
   },
   "subject":{
      "id":"taxid:1755642",
      "taxid":1755642,
      "name":"lachnospiraceae bacterium mc_35",
      "type":"biolink:Bacterium",
      "scientific_name":"lachnospiraceae bacterium mc_35",
      "parent_taxid":186928,
      "lineage":[
         1755642,
         186928,
         186803,
         3085636,
         186801,
         1239,
         1783272,
         2,
         131567,
         1
      ],
      "rank":"species"
   }
}
```

Record shows healthy_sample_size is 0, which is excluded: <br>
```ruby
{
   "_id":"1872678_OrganismalEntityAsAModelOfDiseaseAssociation_D000855",
   "association":{
      "predicate":"OrganismalEntityAsAModelOfDiseaseAssociation",
      "control_name":"healthy control",
      "qualifier":"increase",
      "qualifier_ratio":"0.003266555",
      "disease_sample_size":"4",
      "disease_abundance_mean":"0.00459856",
      "disease_abundance_median":"0.003266555",
      "disease_abundance_sd":"0.002894081",
      "healthy_sample_size":"0",
      "healthy_abundance_mean":"0",
      "healthy_abundance_median":"0",
      "healthy_abundance_sd":"0"
   },
   "object":{
      "id":"MESH:D000855",
      "name":"anorexia",
      "mesh":"D000855",
      "type":"biolink:Disease",
      "description":"The lack or loss of APPETITE accompanied by an aversion to food and the inability to eat. It is the defining characteristic of the disorder ANOREXIA NERVOSA."
   },
   "subject":{
      "id":"taxid:1872678",
      "taxid":1872678,
      "name":"pseudoalteromonas fuliginea",
      "type":"biolink:Bacterium",
      "scientific_name":"pseudoalteromonas fuliginea",
      "parent_taxid":53246,
      "lineage":[
         1872678,
         53246,
         267888,
         135622,
         1236,
         1224,
         3379134,
         2,
         131567,
         1
      ],
      "rank":"species"
   }
}
```

## Metabolite-Gene
Record example: <br>
```ruby
{
   "_id":"6912_associated_with_7498",
   "association":{
      "predicate":"biolink:associated_with",
      "score":0.891,
      "sources":[
         "host",
         "microbiota",
         "food related",
         "drug related"
      ],
      "infores":[
         "stitch"
      ]
   },
   "object":{
      "id":"NCBIGene:7498",
      "symbol":"XDH",
      "type":"biolink:Gene",
      "entrezgene":7498,
      "protein_size":1333,
      "xrefs":{
         "ensemblgene":"ENSG00000158125",
         "hgnc":12805,
         "uniprotkb":"P47989"
      },
      "name":"xanthine dehydrogenase",
      "summary":"Xanthine dehydrogenase belongs to the group of molybdenum-containing hydroxylases involved in the oxidative metabolism of purines. The encoded protein has been identified as a moonlighting protein based on its ability to perform mechanistically distinct functions. Xanthine dehydrogenase can be converted to xanthine oxidase by reversible sulfhydryl oxidation or by irreversible proteolytic modification. Defects in xanthine dehydrogenase cause xanthinuria, may contribute to adult respiratory stress syndrome, and may potentiate influenza infection through an oxygen metabolite-dependent mechanism. [provided by RefSeq, Jan 2014]"
   },
   "subject":{
      "id":"PUBCHEM.COMPOUND:6912",
      "name":"xylitol",
      "type":"biolink:SmallMolecule",
      "pubchem_cid":6912,
      "drug_name":"xylitol",
      "chemical_formula":"C5H12O5",
      "smiles":"C(C(C(C(CO)O)O)O)O",
      "xrefs":{
         "kegg_compound":"C00379",
         "hmdb":"HMDB0002917",
         "drugbank":"DB11195"
      }
   }
}
```