# GMMAD2_data
Parser for GMMDA2 database (microbes-metabolites, microbes-diseases, and metabolites-diseases)
***
record_example: <br>
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