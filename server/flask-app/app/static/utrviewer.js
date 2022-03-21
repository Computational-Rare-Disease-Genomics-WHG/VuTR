var kozak_colors = {
    Strong: "#E69F00",
    Moderate: "#56B4E9",
    Weak: "#009E73",
    None: "#009E73"
};
var pathogenicity_colors = {
    "Pathogenic": "#D55E00",
    "Likely pathogenic": "#D55E00",
    "Benign": "#0072B2",
    "Likely benign": "#0072B2",
    "Conflicting interpretations": "#CC79A7",
    "Uncertain significance": "#CC79A7",
}
var detail_mapping = {

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


    /// ORF Details
    "ensembl_transcript_id": "Ensembl Transcript ID",
    "orf_start_codon": "ORF Start Codon",
    "orf_seq": "ORF Sequence",
    "orf_stop_codon": "ORF Stop codon",
    "orf_type": "ORF Type",
    "frame": "ORF Frame w.r.t. CDS",
    "kozak_context": "7bp context",
    "context": "11 bp context",
    "kozak_consensus_strength": "Kozak Consensus Strength",
    "orf_id": "ORF ID",

    /// UTR Annotator details
    "variant_id": "Variant ID",
    // uAUG gained mappings
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
    "uSTOP_gained_newSTOPDistanceToCDS": "uSTOP-gained Distance from Stop to CDS"

    /*
    To be defined when we have frameshift variants with indels

    "uFrameshift_ref_type":,
    "uFrameshift_ref_type_length":,
    "uFrameshift_StartDistanceToCDS":,
    "uFrameshift_alt_type":,
    "uFrameshift_alt_type_length":,
    "uFrameshift_KozakContext":,
    "uFrameshift_KozakStrength": "Kozak Strength", */
}

var strand_corrected_interval = function(
    start,
    end,
    start_site,
    buffer,
    strand
) {
    if (strand == '+') {
        return ({
            'start': start,
            'end': end
        })
    } else {
        return ({
            "start": (start_site + buffer) - end,
            "end": (start_site + buffer) - start,
        })
    }
}

/* On click event handlers*/
var open_modal = function(data) {
    // Clear previous
    $("#modal-container").empty();


    var custom = `
		<div class="modal fade" id="feature-modal" tabindex="-1" role="dialog" aria-labelledby="exampleModalLongTitle" aria-hidden="true">
		<div class="modal-dialog " role="document">
		  <div class="modal-content">
			<div class="modal-header">
			  <h5 class="modal-title" id="exampleModalLongTitle">Detail</h5>
			  <button type="button" class="close" data-dismiss="modal" aria-label="Close">
				<span aria-hidden="true">&times;</span>
			  </button>
			</div>
			<div class="modal-body">
			  <ul id="feature-modal-data">
			  </ul>
			</div>

		  </div>
		</div>
	  </div>
		`;
    document.getElementById('modal-container')
        .insertAdjacentHTML('beforeend',
            custom);

    // Filter data
    var ul = document.getElementById('feature-modal-data');
    for (const [key, value] of Object.entries(data)) {
        var li = document.createElement('li');
        ul.appendChild(li);
        // Map the VEP consequence terms to formatted strings
        li.innerHTML += `<b>${detail_mapping[key]}</b>: ${value}`;
    }
    $('#feature-modal').modal();
}

var search_obj = function(data, id, id_var) {
    /**
     * Quick wrapper to function to search through an
     * array of objects
     */
    return (data.filter(e => {
        return (e[id_var] == id)
    })[0])
}


var reverse = function(x) {
    const length = x.length;
    for (let i = 0; i < Math.floor(length / 2); i++) {
        const temp = x[i]
        x = x.substring(0, i) +
            x.charAt(length - i - 1) +
            x.substring(i + 1, length - i - 1) +
            temp +
            x.substring(length - i)
    }
    return x
};




var create_transcript_viewer = function(tr_obj,
    div,
    start_site,
    strand,
    buffer,
    gnomad_utr_impact,
    clinvar_utr_impact) {


    // Subset to the first 100 bases following the CDS
    var sequence = strand == "+" ? tr_obj["full_seq"]
        .substring(0, start_site + buffer) :
        reverse(tr_obj["full_seq"].substring(0, start_site + buffer));


    // Create the feature viewer
    var ft2 = new FeatureViewer.createFeature(sequence, div, {
        showAxis: false,
        showSequence: true,
        brushActive: true,
        toolbar: false,
        bubbleHelp: true,
        zoomMax: 10
    })

    //********** Gene Structure *********/
    // Plot where the coding sequence and 5' utrs
    if (start_site != 0) {
        ft2.addFeature({
            data: [{
                    x: strand_corrected_interval(start_site + 1, start_site + buffer, start_site, buffer, strand)['start'],
                    y: strand_corrected_interval(start_site + 1, start_site + buffer, start_site, buffer, strand)['end'],
                    color: '#909590',
                    id: 'cds_rect',
                },
                {
                    x: strand_corrected_interval(1, start_site, start_site, buffer, strand)['start'],
                    y: strand_corrected_interval(1, start_site, start_site, buffer, strand)['end'],
                    color: '#2C302E',
                    id: 'utr_rect'
                }
            ],
            name: "Gene Structure",
            className: "gene_struct",
            type: "rect",
            color: "#00000"
        });
        document.getElementById("fcds_rect").setAttribute("height", "30");
        document.getElementById("fcds_rect").setAttribute("y", "-8");
    }


    //********** ORFS *********/

    // Plot each separate ORF on a separate track based on frame
    var orf_groups = [{
            "grouping_name": "uORFs",
            "orf_type": "uORF",
            "frame": ["Inframe", "Out-of-Frame (2bp)", "Out-of-Frame (1bp)"]
        },
        {
            "grouping_name": "Inframe (oORF)",
            "orf_type": "oORF",
            "frame": ["Inframe"]
        },
        {
            "grouping_name": "Out-of-Frame (oORF)",
            "orf_type": "oORF",
            "frame": ["Out-of-Frame (2bp)", "Out-of-Frame (1bp)"]
        }
    ]
    var uorfs = tr_obj["orfs"]
    orf_groups.forEach(group => {
        // Filter each ofs based on the grouping criteria
        // defined above
        var curr_orf_type = uorfs.filter(function(obj) {
            return (group.frame.includes(obj.frame) &&
                obj.orf_type == group.orf_type)
        });
        curr_orf_frame_dat = [];
        curr_orf_type.forEach(e => curr_orf_frame_dat.push({
            x: strand_corrected_interval(e.orf_start_codon, e.orf_stop_codon, start_site, buffer, strand)['start'],
            y: strand_corrected_interval(e.orf_start_codon, e.orf_stop_codon, start_site, buffer, strand)['end'],
            color: kozak_colors[e.kozak_consensus_strength],
            id: e.orf_id,
        }))
        if (curr_orf_frame_dat.length > 0) {
            ft2.addFeature({
                data: curr_orf_frame_dat,
                type: "rect",
                className: "uorf_rect",
                name: group.grouping_name,
                color: '#474A48'
            })
        }

    });

    ft2.onFeatureSelected(function(d) {
        open_modal(search_obj(uorfs, d.detail.id, 'orf_id'));
    });
    //********** Variants *********/

    var gnomad_variant_ft = new FeatureViewer.createFeature(sequence,
        '#gnomad_tracks', {
            showAxis: false,
            showSequence: false,
            brushActive: true,
            toolbar: false,
            bubbleHelp: false,
            zoomMax: 10
        })

    var clinvar_variant_ft = new FeatureViewer.createFeature(sequence,
        '#clinvar_tracks', {
            showAxis: false,
            showSequence: false,
            brushActive: true,
            toolbar: false,
            bubbleHelp: false,
            zoomMax: 10
        })
    var clinvar_variants = gnomad_data['clinvar_variants'];
    var clinvar_var_feat_dat = [];
    clinvar_variants.forEach(element => {
        clinvar_var_feat_dat.push({
            x: strand_corrected_interval(element['tpos'], element['tpos'], start_site, buffer, strand)['start'],
            y: strand_corrected_interval(element['tpos'], element['tpos'], start_site, buffer, strand)['end'],
            color: pathogenicity_colors[element.clinical_significance]
        });
    });
    if (clinvar_var_feat_dat.length > 0) {
        clinvar_variant_ft.addFeature({
            data: clinvar_var_feat_dat,
            type: "rect",
            className: "clinvar_var",
            name: "ClinVar Variants",
            color: "#000000"
        });
    }

    var pop_variants = gnomad_data['variants']; // gnomAD variants from the API
    var pop_var_feat_dat = []; // Visualization intervals to pass to feature viewer
    var pop_var_tpos = []; // tmp for storing the transcript positions of the variants
    pop_variants.forEach(element => {
        tpos = element['tpos']
        // Only append to the track if didn't exist
        if (!pop_var_tpos.includes(tpos)) {
            pop_var_tpos.push(tpos);
            strand_corrected_tpos = strand_corrected_interval(element['tpos'], element['tpos'], start_site, buffer, strand)
            pop_var_feat_dat.push({
                x: strand_corrected_tpos['start'],
                y: strand_corrected_tpos['end'],
            });
        }


    });
    gnomad_variant_ft.addFeature({
        data: pop_var_feat_dat,
        type: "rect",
        className: "gnomAD_var",
        name: "gnomAD Variants",
        color: "#424242"
    });

    clinvar_utr_impact.forEach(
        element => {
            clinvar_variant_ft.addFeature({
                data: [{
                    x: strand_corrected_interval(element['start'], element['end'], start_site, buffer, strand)['start'],
                    y: strand_corrected_interval(element['start'], element['end'], start_site, buffer, strand)['end'],
                    id: element.variant_id
                }],
                type: "rect",
                className: "clinvar_high_impact_variant",
                name: element.variant_id,
                color: kozak_colors[element.kozak_strength]
            })
        }
    )

    gnomad_utr_impact.forEach(
        element => {
            gnomad_variant_ft.addFeature({
                data: [{
                    x: strand_corrected_interval(element['start'], element['end'], start_site, buffer, strand)['start'],
                    y: strand_corrected_interval(element['start'], element['end'], start_site, buffer, strand)['end'],
                }],
                type: "rect",
                className: "gnomAD_high_impact_variant",
                name: element.variant_id,
                color: kozak_colors[element.kozak_strength]
            })
        }
    )

    clinvar_variant_ft.onFeatureSelected(function(d) {
        open_modal(search_obj(clinvar_variants, d.detail.id, 'variant_id'))
    });
    gnomad_variant_ft.onFeatureSelected(function(d) {
        open_modal(search_obj(pop_variants, d.detail.id, 'variant_id'))

    });
}

var initialize_user_viewer = function(div, tr_obj, start_site, strand, buffer) {
    var sequence = strand == "+" ? tr_obj["full_seq"].substring(0, start_site + buffer) : reverse(tr_obj["full_seq"].substring(0, start_site + buffer));
    user_viewer = new FeatureViewer.createFeature(sequence, div, {
        showAxis: false,
        showSequence: false,
        brushActive: true,
        toolbar: false,
        bubbleHelp: true,
        zoomMax: 10
    })
    return user_viewer
}

var add_user_supplied_feature = function(
    user_viewer, payload, var_name,
    start_site, buffer, strand) {
    dat = payload["data"]["intervals"]
    intervals = [{
        x: strand_corrected_interval(dat['start'], dat['end'], start_site, buffer, strand)['start'],
        y: strand_corrected_interval(dat['start'], dat['end'], start_site, buffer, strand)['end'],
        id: dat['variant_id']
    }]

    user_viewer.addFeature({
        data: intervals,
        type: "rect",
        className: "user_var" + var_name,
        name: var_name,
        color: kozak_colors[dat.kozak_strength]
    });


}
