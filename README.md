# GMMAD2_data
Parser for GMMDA2 database (microbes-metabolites, microbes-diseases, and metabolites-diseases)
***

microbe-metabolite record_example: <br>
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

microbe-disease record example: <br>
```ruby
{
   "_id":"1873459_OrganismalEntityAsAModelOfDiseaseAssociation_D003967",
   "association":{
      "predicate":"OrganismalEntityAsAModelOfDiseaseAssociation",
      "control_name":"healthy control",
      "qualifier":"decrease",
      "qualifier_ratio":"-0.00218686",
      "disease_sample_size":"0",
      "disease_abundance_mean":"0",
      "disease_abundance_median":"0",
      "disease_abundance_sd":"0",
      "healthy_sample_size":"2",
      "healthy_abundance_mean":"0.00218686",
      "healthy_abundance_median":"0.00218686",
      "healthy_abundance_sd":"0.001611793"
   },
   "object":{
      "id":"MESH:D003967",
      "name":"diarrhea",
      "mesh":"D003967",
      "type":"biolink:Disease",
      "description":"An increased liquidity or decreased consistency of FECES, such as running stool. Fecal consistency is related to the ratio of water-holding capacity of insoluble solids to total water, rather than the amount of water present. Diarrhea is not hyperdefecation or increased fecal weight."
   },
   "subject":{
      "id":"taxid:1873459",
      "taxid":1873459,
      "name":"blastococcus sp.",
      "type":"biolink:Bacterium",
      "scientific_name":"blastococcus sp.",
      "parent_taxid":2619396,
      "lineage":[
         1873459,
         2619396,
         38501,
         85030,
         1643682,
         1760,
         201174,
         1783272,
         2,
         131567,
         1
      ],
      "rank":"species"
   }
}
```