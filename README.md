# GMMAD2_data
Parser for GMMDA2 database (microbes-metabolites, microbes-diseases, and metabolites-diseases)
***
record_example: <br>
```ruby
{
   "_id":"6287_associated_with_264203",
   "association":{
      "predicate":"biolink:associated_with",
      "infores":"wom",
      "sources":[
         "Host",
         "Microbiota",
         "Food related",
         "Drug related",
         "Environment"
      ]
   },
   "object":{
      "id":"PUBCHEM.COMPOUND:6287",
      "name":"valine",
      "type":"biolink:ChemicalEntity",
      "pubchem_cid":6287,
      "kegg":"C00183",
      "hmdb":"HMDB0000883",
      "chemical_formula":"C5H11NO2",
      "smiles":"CC(C)C(C(=O)O)N"
   },
   "subject":{
      "id":"taxid:264203",
      "name":"zymomonas mobilis strain zm4 (atcc 31821)",
      "type":"biolink:OrganismalEntity",
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
```