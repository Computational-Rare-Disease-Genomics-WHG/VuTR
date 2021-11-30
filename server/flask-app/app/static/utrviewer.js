
var strand_corrected_interval = function (start, end, start_site, buffer, strand){
    if (strand == '+'){
    return ({
        'start' : start, 
        'end' : end
    })
    }
    else {
        return (
            {
            "start": (start_site+buffer)-end, 
            "end" : (start_site+buffer)-start, 
            }
        )
    }
}

var reverse = function (x) {
    const length = x.length;
    for(let i = 0; i < Math.floor(length / 2); i++) {
      const temp = x[i]
      x = x.substring(0, i) + 
        x.charAt(length - i - 1) +  
        x.substring(i + 1, length - i - 1) + 
        temp + 
        x.substring(length - i)
    }
    return x
  };
  


var create_transcript_viewer = function (tr_obj,
	div,
	start_site,
    strand,
	buffer,
	gnomad_utr_impact,
	clinvar_utr_impact) {


	// Subset to the first 100 bases following the CDS
	var sequence = strand=="+" ? tr_obj["full_seq"].substring(0, start_site + buffer) : reverse(tr_obj["full_seq"].substring(0, start_site + buffer));

	var kozak_colors = {
		Strong: "#E69F00",
		Moderate: "#56B4E9",
		Weak: "#009E73",
		None: "#009E73"
	};
    var pathogenicity_colors = {
        "Pathogenic" : "#D55E00",
        "Likely pathogenic" : "#D55E00",
        "Benign" : "#0072B2",
        "Likely benign" : "#0072B2",
        "Conflicting interpretations" : "#CC79A7",
        "Uncertain significance" : "#000000",
    }

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
					id: 'cds_rect'
				},
				{
					x: strand_corrected_interval( 1, start_site, start_site, buffer, strand)['start'],
					y: strand_corrected_interval( 1, start_site, start_site, buffer, strand)['end'],
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
	var orf_groups  = [
    {
            "grouping_name" : "uORFs", 
            "orf_type" : "uORF", 
            "frame" : ["Inframe", "Out-of-Frame (2bp)", "Out-of-Frame (1bp)"]
    }, 
    {
            "grouping_name" : "Inframe (oORF)", 
            "orf_type" : "oORF", 
            "frame" : ["Inframe"]
    }, 
    {
            "grouping_name" : "Out-of-Frame (oORF)",
            "orf_type" : "oORF", 
            "frame" : ["Out-of-Frame (2bp)", "Out-of-Frame (1bp)"]
    }]
	var uorfs = tr_obj["orfs"]
	orf_groups.forEach(group => {
        
        // Filter based on groups
		var curr_orf_type = uorfs.filter(function (obj) {
			return (group.frame.includes(obj.frame) &&
                    obj.orf_type == group.orf_type)
		});
		curr_orf_frame_dat = [];
		curr_orf_type.forEach(e => curr_orf_frame_dat.push({
			x: strand_corrected_interval(e.orf_start_codon, e.orf_stop_codon , start_site, buffer, strand)['start'], 
			y: strand_corrected_interval(e.orf_start_codon, e.orf_stop_codon , start_site, buffer, strand)['end'],
			color: kozak_colors[e.kozak_consensus_strength]
		}))
        if (curr_orf_frame_dat.length > 0){
		ft2.addFeature({
			data: curr_orf_frame_dat,
			type: "rect",
			className: "uorf_rect",
			name: group.grouping_name,
			color: '#474A48'
		})
    }

	});

    //********** Variants *********/

	var variant_ft= new FeatureViewer.createFeature(sequence, 
        '#variant_viewer', 
        {
		showAxis: false,
		showSequence: false,
		brushActive: true,
		toolbar: false,
		bubbleHelp: false,
		zoomMax: 10
	})

	var clinvar_variants = gnomad_data['clinvar_variants'];
    console.log(clinvar_variants)
	var clinvar_var_feat_dat = [];
	clinvar_variants.forEach(element => {
		clinvar_var_feat_dat.push({
			x: strand_corrected_interval(element['tpos'],element['tpos'],start_site, buffer, strand)['start'],  
			y: strand_corrected_interval(element['tpos'],element['tpos'],start_site, buffer, strand)['end'], 
            color : pathogenicity_colors[element.clinical_significance]
		});
	});
    if (clinvar_var_feat_dat.length>0){
        variant_ft.addFeature({
            data: clinvar_var_feat_dat,
            type: "rect",
            className: "clinvar_var",
            name: "ClinVar Variants",
            color: "#000000"
        });
}

	var pop_variants = gnomad_data['variants'];
	var pop_var_feat_dat = [];
	pop_variants.forEach(element => {
		pop_var_feat_dat.push({
			x: strand_corrected_interval(element['tpos'],element['tpos'],start_site, buffer, strand)['start'],
			y: strand_corrected_interval(element['tpos'],element['tpos'],start_site, buffer, strand)['end'],
		});
	});
	variant_ft.addFeature({
		data: pop_var_feat_dat,
		type: "rect",
		className: "gnomAD_var",
		name: "gnomAD Variants",
		color: "#424242"
	});


	clinvar_utr_impact.forEach(
		element => {
			variant_ft.addFeature({
				data: [{
					x: strand_corrected_interval(element['start'],element['end'],start_site, buffer, strand)['start'],
					y: strand_corrected_interval(element['start'],element['end'],start_site, buffer, strand)['end'],
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
			variant_ft.addFeature({
				data: [{
					x: strand_corrected_interval(element['start'],element['end'],start_site, buffer, strand)['start'],
					y: strand_corrected_interval(element['start'],element['end'],start_site, buffer, strand)['end'],
				}],
				type: "rect",
				className: "gnomAD_high_impact_variant",
				name: element.variant_id,
				color: kozak_colors[element.kozak_strength]
			})
		}
	)

}