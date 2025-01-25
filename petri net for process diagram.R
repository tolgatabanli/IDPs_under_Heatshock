library(DOT)

dot_script <- "
digraph PetriNet {
    rankdir=LR;

    cif [shape=circle, label=\"AlphaFold\\n.cif files\"];
    prep_cif [shape=box, label=\"Extract local\\npLDDT scores\"];
    plddt [shape=circle, label=\"local pLDDTs\"];
    filter [shape=box, label=\"Filter\"];
    filtered_cif [shape=circle, label=\"Inferred IDPs\\nfrom AlphaFold\"];

    uniprotkb [shape=circle, label=\"UniProtKB\\nproteome\"];
    aiupred [shape=box, label=\"AIUPred\"];
    disorder_scores [shape=circle, label=\"disorder_scores\"]
    filter2 [shape=box, label=\"Filter\"];
    filtered_up [shape=circle, label=\"Inferred IDPs\\n from AIUPred\"];
    
    compare [shape=box, label=\"Compare - Intersect\"];
    decided_idps [shape=circle, label=\"IDP Decision\"]
    convert [shape=box, label=\"Convert to\\nEnsembl Gene ID\"]
    product [shape=circle, label=\"IDPs with\\nEnsembl Acc.\"]



    cif -> prep_cif;
    prep_cif -> plddt;
    plddt -> filter;
    filter -> filtered_cif;
    filtered_cif;

    uniprotkb -> aiupred;
    aiupred -> disorder_scores;
    disorder_scores -> filter2
    filter2 -> filtered_up;
    filtered_up;

    filtered_cif -> compare;
    filtered_up -> compare;
    compare -> decided_idps;
    decided_idps -> convert;
    convert -> product

    
}
"

dot(dot_script, file = "Desktop/petri net of pipeline.svg")

# Deseq flow

dot_deseq <- "
digraph PetriNet {
    rankdir=LR;
    
    product [shape=circle, label=\"IDPs with\\nEnsembl Acc.\"]
    counts [shape=circle, label=<RNAseq Counts>]
    deseq [shape=box, label=<DESeq>]
    ds_results [shape=circle, label=<<b>List </b>of<br/>DESeq Results>]
    
    highlight [shape=box, label=<Highlight IDPs<br/>in Results>]
    sig_idps [shape=circle, label=<Significant IDPs>]
    reg [shape=circle, label=<Expression<br/>Change<br/>Comparisons>]
    enrichment [shape=box, label=<Gene Enrichment>]
    fxn [shape=circle, label=<Functions of<br/>Significant IDPs>]
    
    counts -> deseq
    deseq -> ds_results
    product -> highlight
    ds_results -> highlight
    
    highlight -> sig_idps
    highlight -> reg
    sig_idps -> enrichment
    enrichment -> fxn
    
}

"

dot(dot_deseq, file = "../petri net of deseq.svg")


