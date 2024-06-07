# GMMAD2_data
Parser for GMMDA2 database (microbes-metabolites, microbes-diseases, and metabolites-diseases)
***
record_example: <br>
```ruby{
   "_id":"23657835_associated_with_1089454",
   "association":{
      "predicate":"biolink:associated_with",
      "infores":"vmh",
      "source":[
         "Host",
         " Microbiota",
         " Food related"
      ]
   },
   "object":{
      "id":"PUBCHEM.COMPOUND:23657835",
      "name":"adenosylcob(iii)yrinic acid a,c-diamide",
      "type":"biolink:ChemicalEntity",
      "pubchem_cid":23657835,
      "kegg":"C06506",
      "hmdb":"HMDB0001083",
      "chemical_formula":"C55H73CoN11O15+",
      "smiles":"CC1=C2C(C(C([N-]2)C3(C(C(C(=N3)C(=C4C(C(C(=N4)C=C5C(C(C1=N5)CCC(=O)O)(C)C)CCC(=O)O)(C)CC(=O)N)C)CCC(=O)O)(C)CC(=O)N)C)CC(=O)O)(C)CCC(=O)O.[CH2-]C1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)O.[Co+3]"
   },
   "subject":{
      "id":"taxid:1089454",
      "name":"gordonia terrae nbrc 100016",
      "type":"biolink:OrganismalEntity",
      "taxid":1089454,
      "scientific_name":"gordonia terrae nbrc 100016",
      "parent_taxid":2055,
      "lineage":[
         1089454,
         2055,
         2053,
         85026,
         85007,
         1760,
         201174,
         1783272,
         2,
         131567,
         1
      ],
      "rank":"strain"
   }
}```