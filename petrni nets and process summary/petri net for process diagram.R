library(DOT)

dot_script <- "
digraph PetriNet {
    rankdir=LR;

    cif [shape=circle, label=\"AlphaFold .cif files\"];
    prep_cif [shape=box, label=\"Extract pLDDT scores\"];
    plddt [shape=circle, label=\"pLDDT scores\"];
    filter [shape=box, label=\"filter <50\"];
    filtered_cif [shape=circle, label=\"inferred IDPs\n from cif\nin UP acc.\"];

    uniprotkb [shape=box, label=\"UniProtKB proteome\"];
    aiupred [shape=circle, label=\"AIUPred\"];
    filter2 [shape=box, label=\"Statistic filter\"];
    filtered_up [shape=circle, label=\"inferred IDPs\n from UniProt\nin UP acc.\"];
    
    compare [shape=box, label=\"Compare, Intersect/Union\"];
    decided_idps [shape=circle, label=\"Decided IDPs\"]
    convert [shape=box, label=\"Convert to YEAST gene_id\"]
    product [shape=circle, label=\"IDPs with SGD acc.\"]



    cif -> prep_cif;
    prep_cif -> plddt;
    plddt -> filter;
    filter -> filtered_cif;
    filtered_cif;

    uniprotkb -> aiupred;
    aiupred -> filter2;
    filter2 -> filtered_up;
    filtered_up;

    filtered_cif -> compare;
    filtered_up -> compare;
    compare -> decided_idps;
    decided_idps -> convert;
    convert -> product

    
}
"

dot(dot_script, file = "petri net of pipeline.svg")

