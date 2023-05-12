const kozak_colors = {
    Strong: "#E69F00",
    Moderate: "#56B4E9",
    Weak: "#999999",
    None: "#999999"
};

const smorf_sources = {
    'sorfDB' : 'sorfDB | Olexiouk, Volodimir, et al. "sORFs.org: a repository of small ORFs identified by ribosome profiling." Nucleic acids research 44.D1 (2016): D324-D329.', 
    'ribotaper' : 'Ribotaper | Calviello, Lorenzo, et al. "Detecting actively translated open reading frames in ribosome profiling data." Nature methods 13.2 (2016): 165-170.', 
    'ribotish' : 'Ribo-TISH | Zhang, Peng, et al. "Genome-wide identification and differential analysis of translational initiation." Nature communications 8.1 (2017): 1749.',
    'PRICE' : 'PRICE | Erhard, Florian, et al. "Improved Ribo-seq enables identification of cryptic translation events." Nature methods 15.5 (2018): 363-366.'
}

const pathogenicity_colors = {
    "Pathogenic": "#D55E00",
    "Likely pathogenic": "#D55E00",
    "Pathogenic/Likely pathogenic": "#D55E00",

    "Benign": "#009E73",
    "Likely benign": "#009E73",
    "Benign/Likely benign": "#009E73",

    "Conflicting interpretations of pathogenicity": "#CC79A7",
    "Conflicting interpretations": "#CC79A7",
    "Uncertain significance": "#CC79A7",
};

const detail_mapping = {
    // gnomAD mappings
    "alt": "ALT",
    "clinical_significance": "Clinical Significance",
    "clinvar_variation_id": "ClinVar Variation ID",
    "gold_stars": "Gold Stars",
    "hgvsc": "HGVSC",
    "in_gnomad": "Is variant in gnomAD?",
    "major_consequence": "Major Consequence (VEP)",
    "pos": "Position",
    "ref": "REF",
    "review_status": "Review Status",
    "tpos": "Transcript Position",
    "transcript_id": "Ensembl Transcript ID",
    "genome.af": "gnomAD v3 AF",
    "genome.ac": "gnomAD v3 AC",
    "genome.an": "gnomAD v3 AN",
    "efficiency": "Translational Efficiency",
    "lower_bound": "Translational Efficiency (Lower bound)",
    "upper_bound": "Translational Efficiency (Upper bound)",

    /// ORF Details
    "ensembl_transcript_id": "Ensembl Transcript ID",
    "orf_start_codon": "Transcript pos. Start Codon",
    "orf_seq": "ORF Sequence",
    "orf_stop_codon": "Transcript pos. Stop codon",
    "orf_start_codon_genome": "Genomic pos. Start Codon",
    "orf_stop_codon_genome": "Genomic pos. Stop Codon",
    "orf_type": "ORF Type",
    "frame": "ORF Frame w.r.t. CDS",
    "kozak_context": "7bp context",
    "context": "11 bp context",
    "kozak_consensus_strength": "Kozak Consensus Strength",
    "orf_id": "ORF ID",

    /// UTR Annotator details

    // uAUG gained mappings
    "variant_id": "Variant ID",


    "uAUG_gained_CapDistanceToStart": "uAUG-gained Distance from Cap to start",
    "uAUG_gained_DistanceToCDS": "uAUG-gained Distance to CDS",
    "uAUG_gained_DistanceToStop": "uAUG-gained Distance to Stop codon",
    "uAUG_gained_KozakContext": "uAUG-gained Kozak Context",
    "uAUG_gained_KozakStrength": "uAUG-gained Kozak Strength",
    "uAUG_gained_type": "uAUG-gained Type",

    // uAUG Lost
    "uAUG_lost_type": "uAUG-Lost Type",
    "uAUG_lost_KozakContext": "uAUG-Lost Kozak Context",
    "uAUG_lost_KozakStrength": "uAUG-Lost Kozak Strength",
    "uAUG_lost_CapDistanceToStart": "uAUG-Lost Distance from Cap to start",
    "uAUG_lost_DistanceToCDS": "uAUG-Lost Distance to CDS",
    "uAUG_lost_DistanceToSTOP": "uAUG-Lost Distance to Stop",

    // uStop Lost mapping
    "uSTOP_lost_AltStop": "uSTOP-Lost ALT Stop",
    "uSTOP_lost_AltStopDistanceToCDS": "uSTOP-Lost ALT Stop Distance to CDS",
    "uSTOP_lost_KozakContext": "uSTOP-Lost Kozak Context",
    "uSTOP_lost_KozakStrength": "uSTOP-Lost Kozak Strength",
    "uSTOP_lost_FrameWithCDS": "uSTOP-Lost Frame w.r.t. CDS",

    // uSTOP gained mappings
    "uSTOP_gained_ref_StartDistanceToCDS": "uSTOP-gained Distance to CDS",
    "uSTOP_gained_ref_type": "uSTOP-gained REF Type",
    "uSTOP_gained_KozakContext": "uSTOP-gained Kozak Context",
    "uSTOP_gained_KozakStrength": "uSTOP-gained Kozak Strength",
    "uSTOP_gained_newSTOPDistanceToCDS": "uSTOP-gained Distance from Stop to CDS",

    // uFrameshift mapping
    "uFrameshift_ref_type": "uFrameshift Ref Type",
    "uFrameshift_ref_type_length": "uFrameshift Ref Type Length",
    "uFrameshift_StartDistanceToCDS": "uFrameshift Start Distance to CDS",
    "uFrameshift_alt_type": "uFrameshift ALT Type",
    "uFrameshift_alt_type_length": "uFrameshift ALT Type Length",
    "uFrameshift_KozakContext": "uFrameshift Kozak Context",
    "uFrameshift_KozakStrength": "uFrameshift Kozak Strength", 
};

const possible_utr_annotations = Object.keys(detail_mapping).filter(e=> e[0]=="u" && e != "upper_bound")
