library(DOT)

dot_script <- "
digraph PetriNet {
    rankdir=LR;
    beautify=true
    
    input_idps [shape=circle, label=\"IDPs (SGD acc.)\"];
    input_hs [shape=circle, label=\"HS expression data\"];

    join [shape=box, label=\"Join/Filter\"];
    joined [shape=circle, label=\"Expr. Data for IDPs\"]
    
    deseq [shape=box, label=\"DESeq\"];
    deseq_product [shape=circle, label=\"Differential Expression\"]
    pw [shape=circle, label=\"Pathway DB\"]
    
    enrich [shape=box, label=\"Enrichment Analysis\"];
    
    input_idps -> join;
    input_hs -> join;
    join -> joined;
    joined -> deseq;
    deseq -> deseq_product;
    
    deseq_product -> enrich;
    pw -> enrich;
}

"
dot(dot_script)
