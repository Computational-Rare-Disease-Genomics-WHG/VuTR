var create_transcript_viewer = function (tr_obj,
	div,
	start_site,
	buffer,
	gnomad_utr_impact,
	clinvar_utr_impact) {


	// Subset to the first 100 bases following the CDS
	var sequence = tr_obj["full_seq"].substring(0, start_site + buffer);

	var kozak_colors = {
		Strong: "#D55E00",
		Moderate: "#0072B2",
		Weak: "#009E73",
		None: "#009E73"
	};

    // Create the feature viewer
	var ft2 = new FeatureViewer.createFeature(sequence, div, {
		showAxis: false,
		showSequence: true,
		brushActive: true,
		toolbar: true,
		bubbleHelp: true,
		zoomMax: 10
	})

	// Plot where the coding sequence and 5' utrs 
	if (start_site != 0) {
		ft2.addFeature({
			data: [{
					x: start_site + 1,
					y: start_site + buffer,
					color: '#909590',
					id: 'cds_rect'
				},
				{
					x: 1,
					y: start_site,
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
			x: e.orf_start_codon,
			y: e.orf_stop_codon,
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



	var clinvar_variants = gnomad_data['clinvar_variants'];
	var clinvar_var_feat_dat = [];
	clinvar_variants.forEach(element => {
		clinvar_var_feat_dat.push({
			x: element['tpos'],
			y: element['tpos']
		});
	});
    if (clinvar_var_feat_dat.length>0){
        ft2.addFeature({
            data: clinvar_var_feat_dat,
            type: "rect",
            className: "clinvar_var",
            name: "ClinVar Variants",
            color: "#FFA69E"
        });
}

	var pop_variants = gnomad_data['variants'];
	var pop_var_feat_dat = [];
	pop_variants.forEach(element => {
		pop_var_feat_dat.push({
			x: element['tpos'],
			y: element['tpos']
		});
	});
	ft2.addFeature({
		data: pop_var_feat_dat,
		type: "rect",
		className: "gnomAD_var",
		name: "gnomAD Variants",
		color: "#424242"
	});


	clinvar_utr_impact.forEach(
		element => {
			ft2.addFeature({
				data: [{
					x: element.start,
					y: element.end
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
			ft2.addFeature({
				data: [{
					x: element.start,
					y: element.end
				}],
				type: "rect",
				className: "gnomAD_high_impact_variant",
				name: element.variant_id,
				color: kozak_colors[element.kozak_strength]
			})
		}
	)
}