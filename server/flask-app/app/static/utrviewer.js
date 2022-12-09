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
};

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
    "orf_start_genome": "Genomic pos. Start Codon",
    "orf_stop_genome": "Genomic pos. Stop Codon",
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
};


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


var get_exon_structure = function(genomic_features, buffer, start_site,
    strand) {

    /* Filter to five prime UTR*/
    var genomic_features = genomic_features.filter(e => {
        return (e.type === 'exon');
    });

    /* Sort by exon number */
    genomic_features.sort((a, b) => (a.exon_number > b.exon_number) ? 1 : -
        1);

    var new_x = 1;
    var output = [];
    /* Loop over genomic features */
    genomic_features.forEach((element, index) => {
        /* Find the width of the exon*/

        width_exon = element['end'] - element['start'];
        exon_start = strand_corrected_interval(new_x, new_x +
            width_exon, start_site, buffer, strand)['start']
        exon_end = strand_corrected_interval(new_x, new_x +
            width_exon, start_site, buffer, strand)['end']

        /* Exclude exons that start before the CDS */
        if (exon_end > 0) {
            /* Trim exon that goes beyond buffer a little bit */
            if (exon_start < 0) {
                output.push({
                    'x': 0,
                    'y': Math.min(exon_end, start_site +
                        buffer),
                    color: '#A4AAAC',
                    description: 'Exon ' + (index + 1),

                });

            } else {

                output.push({
                    'x': exon_start,
                    'y': Math.min(exon_end, start_site +
                        buffer),
                    color: '#A4AAAC',
                    description: 'Exon ' + (index + 1),

                });
            }

        }
        /* Add exon length*/
        new_x = new_x + width_exon + 1;

    });

    return (output);
}
var flattenObj = (ob) => {

    // The object which contains the
    // final result
    let result = {};

    // loop through the object "ob"
    for (const i in ob) {

        // We check the type of the i using
        // typeof() function and recursively
        // call the function again
        if ((typeof ob[i]) === 'object' && !Array.isArray(ob[i])) {
            const temp = flattenObj(ob[i]);
            for (const j in temp) {

                // Store temp in result
                result[i + '.' + j] = temp[j];
            }
        }

        // Else store ob[i] in result directly
        else {
            result[i] = ob[i];
        }
    }
    return result;
};



/* On click event handlers*/
var open_modal = function(data, type) {

    /// Type is the data_type of the detail presented

    // Clear previous modal
    $("#modal-container").empty();


    if (type === 'gnomad') {
        delete data['transcript_consequence'];
        data = flattenObj(data);
    }

    if (type === 'gnomad' || type === 'clinvar') {
        delete data['type'];
        delete data['kozak_strength'];
        delete data['start'];
        delete data['end'];
        delete data['viz_color'];
        delete data['viz_type'];
    }


    var custom = `
		<div class="modal fade" id="feature-modal" tabindex="-1" role="dialog" aria-labelledby="exampleModalLongTitle" aria-hidden="true">
		<div class="modal-dialog modal-lg" role="document">
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
    if (type === 'clinvar') {
        var li = document.createElement('li');
        ul.appendChild(li);
        li.innerHTML +=
            `<b>View variant in ClinVar</b>: <a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/${data['clinvar_variation_id']}'>${data['clinvar_variation_id']}</a>`;
    }

    if (type === 'gnomad') {
        var li = document.createElement('li');
        ul.appendChild(li);
        li.innerHTML +=
            `<b>View variant in gnomad</b>: <a href='https://gnomad.broadinstitute.org/variant/${data['variant_id']}?dataset=gnomad_r3'>${data['variant_id']}</a>`;
    }

    // Add TE
    var dat = data;
    if ('efficiency' in data) {
        var te_chartdiv = document.createElement('li');
        ul.appendChild(te_chartdiv);
        te_chartdiv.innerHTML +=
            `<b>Translational Efficiency Distribution</b><br/><canvas id = 'te-plot'></canvas`;
        const ctx = document.getElementById('te-plot').getContext('2d');
        // Histogram data from R
        const labels = ["10-20", "20-30,", "30-40", "40-50", "50-60",
            "60-70", "70-80", "80-90", "90-100", "100-110", "110-120",
            "120-130", "130-140", "140-150"
        ]
        const counts = [10, 326, 795, 1230, 1628, 2182, 3563, 3510, 2600,
            1687, 811, 239, 48, 7
        ]
        var bg_color = 
        ['#134DF1',
        '#134DF1',
        '#134DF1',
        '#134DF1',
        '#134DF1',
        '#134DF1',
        '#134DF1',
        '#134DF1',
        '#134DF1',
        '#134DF1',
        '#134DF1',
        '#134DF1',
        '#134DF1',
        '#134DF1',
        '#134DF1'];

        // Find index where efficiency falls under 
        for (var i = 0; i < labels.length; ++i) {
            const r = labels[i].split('-');
            
            if (parseInt(r[0]) < dat['efficiency'] &&
                dat['efficiency'] < parseInt(r[1])) {
                bg_color[i] = '#E83A5A';
            }
        }

        const data = {
            labels: labels,
            datasets: [{
                label: 'Translational Efficiency of all uORFs',
                data: counts,
                borderWidth: 0,
                backgroundColor: bg_color,
                //borderColor: 'rgb(75, 192, 192)',
            }]
        };
        // Create configurationT
        const config = {
            type: 'bar',
            data: data,
            options: {
                scales: {
                    y: {
                        beginAtZero: true
                    }
                },
                legend: {
                    display: false
                },
                tooltip: {
                    enabled: false // <-- this option disables tooltips
                  }            
        

            },
        };
        new Chart(ctx, config);
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

var handleZoom = function(ft, d) {
    var start = d.detail.start;
    var end = d.detail.end;
    ft.zoom(start, end);

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




var create_transcript_viewer = function(
    tr_obj,
    div,
    start_site,
    strand,
    buffer,
    seq,
    gnomad_utr_impact,
    clinvar_utr_impact,
    genomic_features
) {

    // Subset to the first 100 bases following the CDS
    var sequence = strand == "+" ? seq
        .substring(0, start_site + buffer) :
        reverse(seq.substring(0, start_site + buffer));


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
                    x: strand_corrected_interval(start_site + 1,
                        start_site + buffer, start_site,
                        buffer, strand)['start'],
                    y: strand_corrected_interval(start_site + 1,
                        start_site + buffer, start_site,
                        buffer, strand)['end'],
                    color: '#58565F',
                    description: "\t\tCDS",
                    id: 'cds_rect',
                },
                {
                    x: strand_corrected_interval(1, start_site,
                        start_site, buffer, strand)[
                        'start'],
                    y: strand_corrected_interval(1, start_site,
                        start_site, buffer, strand)['end'],
                    color: '#A4AAAC',
                    description: "\t\t5' UTR",
                    id: 'utr_rect'
                }
            ],
            name: "Gene Structure",
            className: "gene_struct",
            type: "rect",
            color: "#A4AAAC"
        });
        document.getElementById("fcds_rect").setAttribute("height", "30");
        document.getElementById("fcds_rect").setAttribute("y", "-8");

        ft2.addFeature({
            data: get_exon_structure(genomic_features, buffer,
                start_site, strand),
            color: '#2C302E',
            name: 'Exon Structure',
            className: 'exon_struct',
            type: 'rect',
            color: '#000000'
        });
    }


    //********** ORFS *********/

    // Plot each separate ORF on a separate track based on frame
    var orf_groups = [{
            "grouping_name": "uORFs",
            "orf_type": "uORF",
            "frame": ["Inframe", "Out-of-Frame (2bp)",
                "Out-of-Frame (1bp)"
            ]
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
    var uorfs = tr_obj;
    orf_groups.forEach(group => {
        // Filter each ofs based on the grouping criteria
        // defined above
        var curr_orf_type = uorfs.filter(function(obj) {
            return (group.frame.includes(obj.frame) &&
                obj.orf_type == group.orf_type)
        });
        curr_orf_frame_dat = [];
        curr_orf_type.forEach(e => curr_orf_frame_dat.push({
            x: strand_corrected_interval(e
                .orf_start_codon, e.orf_stop_codon,
                start_site, buffer, strand)[
                'start'],
            y: strand_corrected_interval(e
                .orf_start_codon, e.orf_stop_codon,
                start_site, buffer, strand)['end'],
            color: kozak_colors[e
                .kozak_consensus_strength],
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
        open_modal(search_obj(uorfs, d.detail.id, 'orf_id'), 'orf');
    });
    //********** ClinVar *********/



    /* gnomAD Variant Track */
    var gnomad_variant_ft = new FeatureViewer.createFeature(sequence,
        '#gnomad_tracks', {
            showAxis: false,
            showSequence: false,
            brushActive: false,
            toolbar: false,
            bubbleHelp: true,
            zoomMax: 10
        })

    var pop_variants = gnomad_data[
        'variants']; // gnomAD variants from the API
    var
        pop_var_feat_dat = []; // Visualization intervals to pass to feature viewer
    var
        pop_var_tpos = []; // tmp for storing the transcript positions of the variants
    pop_variants.forEach(element => {
        tpos = element['tpos']
        // Only append to the track if didn't exist
        if (!pop_var_tpos.includes(tpos)) {
            pop_var_tpos.push(tpos);
            strand_corrected_tpos = strand_corrected_interval(
                element['tpos'], element['tpos'], start_site,
                buffer, strand)
            pop_var_feat_dat.push({
                x: strand_corrected_tpos['start'],
                y: strand_corrected_tpos['end'],
                id: element['variant_id'],
                className: 'no_impact_gnomad'
            });
        }


    });
    if (pop_var_feat_dat.length > 0) {
        gnomad_variant_ft.addFeature({
            data: pop_var_feat_dat,
            type: "rect",
            className: "gnomAD_var",
            name: "gnomAD Variants",
            color: "#424242"
        });
    }


    gnomad_utr_impact.forEach(
        element => {
            gnomad_variant_ft.addFeature({
                data: [{
                    x: strand_corrected_interval(
                        element['start'],
                        element['end'],
                        start_site,
                        buffer,
                        strand)['start'],
                    y: strand_corrected_interval(
                        element['start'], element[
                            'end'], start_site,
                        buffer, strand)['end'],
                    id: element["variant_id"]
                }],
                type: "rect",
                className: "gnomAD_high_impact_variant",
                name: element['variant_id'],
                color: kozak_colors[element.kozak_strength]
            })
        }
    )
    gnomad_variant_ft.onFeatureSelected(function(m) {
        gnomad_utr_conseq = search_obj(gnomad_utr_impact, m.detail
            .id, 'variant_id')
        open_modal(Object.assign(search_obj(pop_variants, m.detail
                .id, 'variant_id'), gnomad_utr_conseq),
            'gnomad');



    });
    var clinvar_variant_ft = new FeatureViewer.createFeature(sequence,
        '#clinvar_tracks', {
            showAxis: false,
            showSequence: false,
            brushActive: false,
            toolbar: false,
            bubbleHelp: true,
            zoomMax: 10
        })
    var clinvar_variants = gnomad_data['clinvar_variants'];
    var clinvar_var_feat_dat = [];
    clinvar_variants.forEach(element => {
        clinvar_var_feat_dat.push({
            x: strand_corrected_interval(element['tpos'],
                element['tpos'], start_site, buffer,
                strand)['start'],
            y: strand_corrected_interval(element['tpos'],
                element['tpos'], start_site, buffer,
                strand)['end'],
            color: pathogenicity_colors[element
                .clinical_significance],
            id: element['variant_id'],
            className: 'no_impact_clinvar'

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

    clinvar_utr_impact.forEach(
        element => {
            clinvar_variant_ft.addFeature({
                data: [{
                    x: strand_corrected_interval(
                        element['start'], element[
                            'end'], start_site,
                        buffer, strand)['start'],
                    y: strand_corrected_interval(
                        element['start'], element[
                            'end'], start_site,
                        buffer, strand)['end'],
                    id: element.variant_id
                }],
                type: "rect",
                className: "clinvar_high_impact_variant",
                name: element['variant_id'],
                color: kozak_colors[element.kozak_strength]
            });
        }
    )

    /*Event handler to view the Clinvar variant detail page*/
    clinvar_variant_ft.onFeatureSelected(function(d) {
        clinvar_utr_conseq = search_obj(clinvar_utr_impact, d.detail
            .id, 'variant_id')
        open_modal(Object.assign(
                search_obj(clinvar_variants, d.detail.id,
                    'variant_id'),
                clinvar_utr_conseq),
            'clinvar');
    });

    var viewers = {
        'clinvar_ft': clinvar_variant_ft,
        'gnomad_ft': gnomad_variant_ft,
        'arch_ft': ft2
    };

    return viewers;

}


var initialize_user_viewer = function(div,
    seq,
    start_site,
    strand,
    buffer) {
    var sequence = strand == "+" ? seq.substring(0, start_site + buffer) :
        reverse(seq.substring(0, start_site + buffer));
    user_viewer = new FeatureViewer.createFeature(sequence, div, {
        showAxis: false,
        showSequence: false,
        brushActive: false,
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
        x: strand_corrected_interval(dat['start'], dat['end'],
            start_site, buffer, strand)['start'],
        y: strand_corrected_interval(dat['start'], dat['end'],
            start_site, buffer, strand)['end'],
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