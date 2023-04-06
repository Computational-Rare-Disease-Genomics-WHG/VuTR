/** 
utrviewer.js 

The following functions allow for the creation of the transcript viewer 
in the vutr.rarediseasegenomics.org/viewer/<ENST> page.

Namely, utilising the features from CalipoSIB's Feature Viewer.
See their Github page for more details.

Following functions add additional functional such as a modal detail upon 
clicking one of the elements, and correcting for the fact that 
FeatureViewer isn't strand aware, whereas VuTR requires separate 
views for reverse strand genes
*/



/** Reverses a string
@param {str}
@returns {str} The string reversed
*/
function reverse(str) {
    return str.split('').reverse().join('');
  }

/**
Processes a sequences and pads introns with X's
TODO: Change X's to intron sequences

@param {pos}
@param {exonBoundaries}
@returns {num}
*/
function processSequence(cdnaSequence, strand, startSite, buffer, exonBoundaries) {
    let seqWithIntronBuffer = '';
    for (let i = 0; i < exonBoundaries.length; i++) {
      const e = exonBoundaries[i];
      // Check if the exon is after the start site + buffer
      if (e.y > startSite + buffer) {
        seqWithIntronBuffer += cdnaSequence.substring(e.x - 1, startSite + buffer);
        break;
      }
      // Pad introns with blank spaces
      seqWithIntronBuffer += cdnaSequence.substring(e.x - 1, e.y) + 'X'.repeat(45);
    }
  
    const sequence = strand === '+' ? seqWithIntronBuffer.toUpperCase() : reverse(seqWithIntronBuffer.toUpperCase());
    return sequence;
  }
  
/**
@param {pos}
@param {exonBoundaries}
@returns {num}
*/
var findExonNumber = function(
    pos, 
    exonBoundaries
){
    var foundExon = exonBoundaries.filter(
        e => {return ((e.x <= pos) && (pos<= e.y))}
    );
    if (foundExon.length==1){
        return (foundExon[0].exon_number)
    }
    return (-1)
}

/**
@param {genomic_features}
@returns {List[Obj]}
*/
var getExonBoundaries = function (
    genomic_features
){
    /* Filter to five prime UTR*/
    var genomic_features = genomic_features.filter(e => {
        return (e.type === 'exon');
    });
    /* Sort by exon number */
    genomic_features.sort((a, b) => (a.exon_number > b.exon_number) ? 1 : -
        1);
    
    output = [];
    x = 1;
    genomic_features.forEach((element, index) => {
        width_exon = element['end'] - element['start'];
        output.push({
            'exon_number' : index+1,
            'x' : x,
            'y' : x+width_exon

        }
        );
        x += width_exon+1; 
    });
    return (output);
}

/**
This function converts a transcript [start, end] to an appropriate coordinate 
for Feature viewer depending on which strand the gene is located on. 

Also ensure that the interval doesn't go out of bounds for Feature Viewer

@param {number} start - The transcript start of the feature
@param {number} end - The transcript end of the feature
@param {number} start_site - The start site of the CDS
@param {number} buffer - The amount of bps of buffer from the CDS to the end of the viz.
@param {string} strand - Which strand is this gene on? "-" or "+"
@returns {Obj} {start : The strand aware start, end: The strand aware end}
*/
var scInterval = function(
    start,
    end,
    start_site,
    buffer,
    strand, 
    exonBoundaries
) {

    var intronSize = 45;
    var startExonNumber = findExonNumber(start, exonBoundaries);
    var endExonNumber = findExonNumber(end, exonBoundaries);
    var startIntronDelta = intronSize*(startExonNumber-1)
    var endIntronDelta = intronSize*(endExonNumber-1)

    if (strand == '+') {
        return ({
            'start': start+startIntronDelta,
            'end': Math.min(end - 1, start_site+buffer+1)+endIntronDelta
        });
    } else {
        return ({
            "start": Math.max((start_site + buffer + 2) - end, 1),
            "end": (start_site + buffer + 1) - start,
        });
    }
}


/*
Determines the exon structure, essential a list of start, and ends (aware of strand)
depending on which is 

@param {[Objs...]} genomic_features - A list of Objs with each obj an exon from the MANE gff 
@param {number} start_site - The start site of the CDS
@param {number} buffer - The amount of bps of buffer from the CDS to the end of the viz.
@param {string} strand - Which strand is this gene on? "-" or "+"
@returns {Obj} [{start : The strand aware start, end: The strand aware end}] for the exons
*/
var getExonStructure = function(
    genomic_features,
    buffer,
    start_site,
    strand, 
    exonBoundaries) {

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

        exon_start = scInterval(new_x, new_x +
            width_exon, start_site, buffer, strand, exonBoundaries)['start']
        exon_end = scInterval(new_x, new_x +
            width_exon, start_site, buffer, strand, exonBoundaries)['end'];

        /* Exclude exons that start before the CDS */
        if (exon_end > 0) {
            /* Trim exon that goes beyond buffer a little bit */
            if (exon_start < 0) {
                if (strand == '-') {
                    output.push({
                        'x': Math.max(buffer + 0.0001,
                            exon_start),
                        'y': Math.min(exon_end + 1,
                            start_site +
                            buffer),
                        color: '#A4AAAC',
                        description: 'Exon ' + (index + 1),
                    });
                } else {
                    output.push({
                        'x': Math.max(0, exon_start),
                        'y': Math.min(exon_end, start_site +
                            buffer),
                        color: '#A4AAAC',
                        description: 'Exon ' + (index + 1)
                    });
                }

            } else {
                if (strand == '-') {
                    output.push({
                        'x': Math.max(buffer + 0.001,
                            exon_start),
                        'y': Math.min(exon_end, start_site +
                            buffer),
                        color: '#A4AAAC',
                        description: 'Exon ' + (index + 1),
                    });
                } else {
                    output.push({
                        'x': Math.max(0, exon_start),
                        'y': Math.min(exon_end, start_site +  
                            0.999),
                        color: '#A4AAAC',
                        description: 'Exon ' + (index + 1),

                    });
                }
            }

        }
        /* Add exon length*/
        new_x = new_x + width_exon + 1;

    });

    return (output);
}



/*
Flattens nested objects
@param {Obj} obj - Where there might be nesting.
@returns {Obj} - Returns an object that is exactly 1 level deep
*/

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


/** 
Creates the <li> elements for the UTR annotation
@param {Obj} data - The obj with output from UTR annotator
@returns {string} - HTML string to be used in the modal
*/

var createUtrAnnotationList = function(data) {

    var ann_obj = Object.entries(data).filter(([key, value]) =>
        possible_utr_annotations.includes(key));
    var annotation = '';
    if (ann_obj.length > 0) {
        var annotation = `<h5>5'UTR Annotation</h5><ul>`;
        // Filter data to only the UTR annotations and then make them into <li> element
        for (const [key, value] of ann_obj) {
            // Map the VEP consequence terms to formatted strings
            annotation = annotation.concat(
                `<li><b>${detail_mapping[key]}</b>: ${value}</li>`);
        }
        annotation = annotation.concat('</ul><hr>');
    }
    return annotation;
}


/* 
On click event handlers to open the modal details
@param {obj} data - The obj with data to be displayed
@param {string} type - The type either 'gnomad', 'clinvar', 'uorf'
@returns None

*/
var openModal = function(data, type) {

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


    var custom;

    // gnomAD 
    if (type == 'gnomad') {
        custom = `
    <div class="modal fade" 
    id="feature-modal"
    tabindex="-1" 
    role="dialog" 
    aria-labelledby="exampleModalLongTitle" 
    aria-hidden="true">

<div class="modal-dialog modal-lg" role="document">
  <div class="modal-content">
    <div class="modal-header">
      <h5 class="modal-title" id="exampleModalLongTitle">ORF Details</h5>
      <button type="button" class="close" data-dismiss="modal" aria-label="Close">
        <span aria-hidden="true">&times;</span>
      </button>
    </div>


    <div class="modal-body">
    <h5>Variant</h5> 
    <ul>
      <li><b>REF</b> : ${data['ref']}</li>
      <li><b>ALT</b> : ${data['alt']} </li>
      <li><b>Genome Position</b> : ${data['pos']} </li>
      <li><b>Transcript Position</b> : ${data['tpos']} </li>
      <li><b>Variant ID</b> : ${data['variant_id']} </li>
      <li><b>HGVSC</b> : ${data['hgvsc']} </li>
    </ul>
    <hr>
    ${createUtrAnnotationList(data)}
      <h5>gnomAD Frequency</h5> 
      <ul>
        <li><b>Allele Count</b> : ${data['genome.ac']}</li>
        <li><b>Allele Number</b> : ${data['genome.an']} </li>
        <li><b>Allele Frequency (all pop)</b> : ${data['genome.af']} </li>
      </ul>
      <hr>
      <b>View variant in gnomAD</b>: <a href='https://gnomad.broadinstitute.org/variant/${data['variant_id']}?dataset=gnomad_r3'>${data['variant_id']}</a>    
    </div>

  </div>
</div>
</div>`;
    }

    if (type=='user-supplied'){
        custom = `
        <div class="modal fade" 
        id="feature-modal"
        tabindex="-1" 
        role="dialog" 
        aria-labelledby="exampleModalLongTitle" 
        aria-hidden="true">
    
    <div class="modal-dialog modal-lg" role="document">
      <div class="modal-content">
        <div class="modal-header">
          <h5 class="modal-title" id="exampleModalLongTitle">Variant ${data['variant_id']} </h5>
          <button type="button" class="close" data-dismiss="modal" aria-label="Close">
            <span aria-hidden="true">&times;</span>
          </button>
        </div>
    
    
        <div class="modal-body">

        ${createUtrAnnotationList(data)}
        </div>
    
      </div>
    </div>
    </div>`;
        

    }

    if (type == 'clinvar') {
        custom = `
    <div class="modal fade" 
    id="feature-modal"
    tabindex="-1" 
    role="dialog" 
    aria-labelledby="exampleModalLongTitle" 
    aria-hidden="true">

<div class="modal-dialog modal-lg" role="document">
  <div class="modal-content">
    <div class="modal-header">
      <h5 class="modal-title" id="exampleModalLongTitle">ORF Details</h5>
      <button type="button" class="close" data-dismiss="modal" aria-label="Close">
        <span aria-hidden="true">&times;</span>
      </button>
    </div>


    <div class="modal-body">
    <h5>Variant</h5> 
    <ul>
      <li><b>REF</b> : ${data['ref']}</li>
      <li><b>ALT</b> : ${data['alt']} </li>
      <li><b>Genome Position</b> : ${data['pos']} </li>
      <li><b>Transcript Position</b> : ${data['tpos']} </li>
      <li><b>Variant ID</b> : ${data['variant_id']} </li>
      <li><b>HGVSC</b> : ${data['hgvsc']} </li>
    </ul>

    <hr>

    <h5>ClinVar variant details</h5> 
    <ul>
    <li><b>Clinical Significance</b> : ${data['clinical_significance']}</li>
    <li><b>Review Status</b> : ${data['review_status']}</li>
    <li><b>Gold Stars</b> : ${data['gold_stars']} </li>
    <li><b>ClinVar Variation ID</b> : ${data['clinvar_variation_id']} </li>
    <li><b>Variant in gnomAD?</b> : ${data['in_gnomad']} </li>
    </ul>
  <hr>
    ${createUtrAnnotationList(data)}
    <b>View variant in ClinVar</b>: <a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/${data['clinvar_variation_id']}'>${data['clinvar_variation_id']}</a>    </div>

  </div>
</div>
</div>

    
    `;

    }

    if (type == 'orf') {
        custom = `
		<div class="modal fade" 
            id="feature-modal"
            tabindex="-1" 
            role="dialog" 
            aria-labelledby="exampleModalLongTitle" 
            aria-hidden="true">

		<div class="modal-dialog modal-lg" role="document">
		  <div class="modal-content">
			<div class="modal-header">
			  <h5 class="modal-title" id="exampleModalLongTitle">ORF Details</h5>
			  <button type="button" class="close" data-dismiss="modal" aria-label="Close">
				<span aria-hidden="true">&times;</span>
			  </button>
			</div>


			<div class="modal-body">
              <h5>Coordinates</h5>
              <ul>
              <li><b>Transcript</b>
                <ul>
                <li><b>Start codon</b> : ${data['orf_start_codon']} </li>
                <li><b>Stop codon</b> : ${data['orf_stop_codon']}</li>
                </ul>      
              </li>
              <li><b>Genome</b>
              <ul>
              <li><b>Start codon</b> : ${data['orf_start_codon_genome']}</li>
              <li><b>Stop codon</b> : ${data['orf_stop_codon_genome']}</li>
              </ul>      

              </li> 
              </ul>
              <hr>

              <h5>ORF details</h5> 
              <ul>
                <li><b>ORF Sequence</b> : <p style='word-wrap: break-word'>${data['orf_seq']}</p></li>
                <li><b>ORF Type</b> : ${data['orf_type']} </li>
                <li><b>ORF Frame</b> : ${data['frame']} </li>
              </ul>
              <hr>


              <h5>Translation details</h5>
              <ul>
                <li><b>Translational Efficiency</b> : ${data['efficiency']} (${data['lower_bound']}-${data['upper_bound']})</li>
                <li><b>Kozak Consensus Strength</b> : ${data['kozak_consensus_strength']} </li>
                <li><b>Kozak Consensus Context</b> :  ${data['kozak_context']}</li>
              </ul>
              <hr>
              <h5>Translation Efficiency Distribution</h5>
              
              <ul id="feature-modal-data">
			  </ul>
			</div>

		  </div>
		</div>
	  </div>
		`;


    }

    if (type == 'smorf') {
        custom = `
		<div class="modal fade" 
            id="feature-modal"
            tabindex="-1" 
            role="dialog" 
            aria-labelledby="exampleModalLongTitle" 
            aria-hidden="true">

		<div class="modal-dialog modal-lg" role="document">
		  <div class="modal-content">
			<div class="modal-header">
			  <h5 class="modal-title" id="exampleModalLongTitle">smORF with evidence</h5>
			  <button type="button" class="close" data-dismiss="modal" aria-label="Close">
				<span aria-hidden="true">&times;</span>
			  </button>
			</div>


			<div class="modal-body">
                <ul>
                <li><b>iORF ID</b> : ${data['smorf_iorf_id']} </li>
                <li><b>Source</b> : Chothani, Sonia P., et al. <a href="https://doi.org/10.1016/j.molcel.2022.06.023" > A high-resolution map of human RNA translation.</a> Molecular Cell 82.15 (2022): 2885-2899. </li>
                <li><b>Length</b> : ${data['len']} bps</li>
                <li><b>Start Codon</b> : ${data['starts']} </li>
			</div>
		  </div>
		</div>
	  </div>
		`;


    }
    document.getElementById('modal-container')
        .insertAdjacentHTML('beforeend',
            custom);

    // Filter data
    var ul = document.getElementById('feature-modal-data');




    // Add TE
    var dat = data;
    if ('efficiency' in data) {
        var te_chartdiv = document.createElement('div');
        ul.appendChild(te_chartdiv);
        te_chartdiv.innerHTML +=
            `<canvas id = 'te-plot'></canvas`;
        const ctx = document.getElementById('te-plot').getContext('2d');
        // Histogram data from R
        const labels = ["10-20", "20-30,", "30-40", "40-50", "50-60",
            "60-70", "70-80", "80-90", "90-100", "100-110", "110-120",
            "120-130", "130-140", "140-150"
        ]
        const counts = [10, 326, 795, 1230, 1628, 2182, 3563, 3510, 2600,
            1687, 811, 239, 48, 7
        ]
        var bg_color = ['#134DF1',
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
            '#134DF1'
        ];

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


/**
 * Quick wrapper to function to search through an
 * array of objects for a certain object 
 * @param {[Obj]} data - The list of objects
 * @param {string} id - The identifier to search
 * @param {string} id_var - The key to search for
 * @returns {Obj or [Obj]}  The object (s) found
 */
var searchObj = function(data, id, id_var) {
    found = data.filter(e => {
        return (e[id_var] == id)
    })
    return (found[0])
}

/**
 * This handles the zooming from the featureViewer object
 * @param {FeatureViewer} ft - The feature viewer to be manipulated
 * @param {FeatureViewer.ClickEvent} - the event with the [start, end] of the zoom.
 * @returns {None}
 */
var handleZoom = function(ft, d) {
    var start = d.detail.start;
    var end = d.detail.end;
    ft.zoom(start, end);
}

/**
 * Quick wrapper to function to reverse complement a string
 * @param {string} x - a string to reverse
 * @returns {string} x - reverse complemented
 */
var reverse = function(x) {
    const length = x.length;
    // Reverse the string 
    for (let i = 0; i < Math.floor(length / 2); i++) {
        const temp = x[i]
        x = x.substring(0, i) +
            x.charAt(length - i - 1) +
            x.substring(i + 1, length - i - 1) +
            temp +
            x.substring(length - i)
    }
    // Complement String 
    var mapObj = {
        a: "t",
        t: "a",
        c: "g",
        g: "c"
    };
    var re = new RegExp(Object.keys(mapObj).join("|"), "gi");
    x_complement = x.replace(re, function(matched) {
        return mapObj[matched];
    });
    return x_complement;
};


/**
 * Creates the main transcript viewer consists of 
 * uORFS, gnomAD and clinvar using Feature Viewer.
 * @param {[Obj]} tr_obj - The uORFs 
 * @param {str} div - id of the <div> to add
 * @param {number} start_site - The start site of the CDS
 * @param {string} strand - Which strand is this gene on? "-" or "+"
 * @param {number} buffer - The amount of bps of buffer from the CDS to the end of the viz.
 * @param {str} seq - The sequence of the cDNA 
 * @param {[Obj]} smorf - The smORF with evidence dataset
 * @param {[Obj]} gnomad_utr_impact - gnomAD dataset for v3.1.2 variants
 * @param {[Obj]} clinvar_utr_impact - The object from gnomAD API
 * @param {[Obj]} genomic_features - The genomic features (exons) from MANE gff
 * @returns {FeatureViewer}
 */
var createTranscriptViewer = function(
    tr_obj,
    div,
    start_site,
    strand,
    buffer,
    seq,
    smorf,
    gnomad_utr_impact,
    clinvar_utr_impact,
    genomic_features, 
    exonBoundaries
) {
    // Process Sequence
    var sequence = processSequence(seq, strand, start_site, buffer, exonBoundaries)

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
                x: scInterval(start_site + 1,
                    start_site + buffer + 2, start_site,
                    buffer, strand, exonBoundaries)['start'],
                y: scInterval(start_site + 1,
                    start_site + buffer + 2, start_site,
                    buffer, strand, exonBoundaries)['end'],
                color: '#58565F',
                description: "CDS",
                id: 'cds_rect',
            }, ...getExonStructure(genomic_features, buffer,
                start_site, strand, exonBoundaries)],
            color: '#2C302E',
            name: 'Gene Structure',
            className: 'exon_struct',
            type: 'rect',
            color: '#000000'
        });
        document.getElementById("fcds_rect").setAttribute("height", "30");
        document.getElementById("fcds_rect").setAttribute("y", "-8");

        if (smorf.length > 0) {
            smorf_data = [];
            smorf.forEach(e => smorf_data.push({
                x: scInterval(
                    e.transcript_start, e.transcript_end +
                    1,
                    start_site, buffer, strand, exonBoundaries
                )['start'],
                y: scInterval(
                    e.transcript_start, e.transcript_end +
                    1,
                    start_site, buffer, strand, exonBoundaries
                )['end'],
                id: e.smorf_iorf_id
            }));

            ft2.addFeature({
                data: smorf_data,
                type: "rect",
                className: "uorf_rect",
                name: "ORFs w. evidence",
                color: '#474A48'
            })
        }
    }


    /* ORFS Variant Track */
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

    // Add ORFs to feature viewer
    orf_groups.forEach(group => {
        // Filter each ofs based on the grouping criteria
        // defined above
        var curr_orf_type = uorfs.filter(function(obj) {
            return (group.frame.includes(obj.frame) &&
                obj.orf_type == group.orf_type)
        });
        curr_orf_frame_dat = [];
        curr_orf_type.forEach(e => curr_orf_frame_dat.push({
            x: scInterval(e
                .orf_start_codon, e.orf_stop_codon,
                start_site, buffer, strand, exonBoundaries)[
                'start'],
            y: scInterval(e
                .orf_start_codon, e.orf_stop_codon,
                start_site, buffer, strand, exonBoundaries)['end'],
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

    // Add event handlers
    ft2.onFeatureSelected(function(d) {
        if (searchObj(smorf, d.detail.id, 'smorf_iorf_id')) {
            openModal(searchObj(smorf, d.detail.id,
                'smorf_iorf_id'), 'smorf');
        } else {
            openModal(searchObj(uorfs, d.detail.id, 'orf_id'),
                'orf');
        }
    });



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
            strand_corrected_tpos = scInterval(
                element['tpos']-1.25, 
                element['tpos']+1,
                start_site,
                buffer, strand, exonBoundaries)
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
                    x: scInterval(
                        element['start'],
                        element['end'],
                        start_site,
                        buffer,
                        strand, exonBoundaries)['start'],
                    y: scInterval(
                        element['start'], element[
                            'end'], start_site,
                        buffer, strand, exonBoundaries)['end'],
                    id: element["annotation_id"]
                }],
                type: "rect",
                className: "gnomAD_high_impact_variant",
                name: element['variant_id'],
                color: kozak_colors[element.kozak_strength]
            })
        }
    )
    gnomad_variant_ft.onFeatureSelected(function(m) {
        let id_sel = m.detail.id;
        if (id_sel.startsWith('annotation')){
            gnomad_utr_conseq = searchObj(gnomad_utr_impact, id_sel, 'annotation_id');
            openModal(Object.assign(searchObj(pop_variants, gnomad_utr_conseq['variant_id'], 'variant_id'), gnomad_utr_conseq),
                'gnomad');
        }
        else{
        gnomad_utr_conseq = searchObj(gnomad_utr_impact, id_sel, 'variant_id')
        openModal(Object.assign(searchObj(pop_variants, id_sel, 'variant_id'), gnomad_utr_conseq),
            'gnomad');
        }


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
            x: scInterval(element['tpos']-1.25,
                element['tpos'] + 1.1, start_site, buffer,
                strand, exonBoundaries)['start'],
            y: scInterval(element['tpos']-1.25,
                element['tpos'] + 1.1, start_site, buffer,
                strand, exonBoundaries)['end'],
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
                    x: scInterval(
                        element['start']-1.25, element[
                            'end']+1.1, start_site,
                        buffer, strand, exonBoundaries)['start'],
                    y: scInterval(
                        element['start']-1.25, element[
                            'end']+1.1, start_site,
                        buffer, strand, exonBoundaries)['end'],
                    id: element.annotation_id
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
        let id_sel = d.detail.id;
        if (id_sel.startsWith('annotation')){
            clinvar_utr_conseq = searchObj(clinvar_utr_impact,
                id_sel, 'annotation_id')    
            openModal(Object.assign(
                    searchObj(clinvar_variants, clinvar_utr_conseq['variant_id'],
                        'variant_id'),
                    clinvar_utr_conseq),
                'clinvar');
        }
        else{
            clinvar_utr_conseq = searchObj(clinvar_utr_impact,
                id_sel, 'variant_id')
            openModal(Object.assign(
                    searchObj(clinvar_variants, id_sel,
                        'variant_id'),
                    clinvar_utr_conseq),
                'clinvar');
        }

    });

    var viewers = {
        'clinvar_ft': clinvar_variant_ft,
        'gnomad_ft': gnomad_variant_ft,
        'arch_ft': ft2
    };

    return viewers;

}

/** 
 * Creates the viewer for user supplied variants
 * @param {str} div - id of the <div> to add
 * @param {str} seq - The sequence of the cDNA (Used for the purposes of length rather than actually plotting)
 * @param {number} start_site - The start site of the CDS
 * @param {number} buffer - The amount of bps of buffer from the CDS to the end of the viz.
 * @param {string} strand - Which strand is this gene on? "-" or "+"
 * @returns {FeatureViewer} The feature viewer object
 */

var initialiseUserViewer = function(
    div,
    seq,
    start_site,
    strand,
    buffer) {
    var sequence = processSequence(seq, strand, start_site, buffer, exonBoundaries)

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


/** 
 *  Adds a track for each user supplied variant 
 * @param {FeatureViewer} user_viewer - The object to add to
 * @param {Obj} payload - The object from the API with the start, end
 * @param {string} var_name - The name of the variant (variant_id)
 * @param {number} start_site - The start site of the CDS
 * @param {number} buffer - The amount of bps of buffer from the CDS to the end of the viz.
 * @param {string} strand - Which strand is this gene on? "-" or "+"
 * @returns {None}
 */

var addUserSuppliedFeature = function(
    user_viewer,
    payload,
    var_name,
    start_site,
    buffer,
    strand, 
    exonBoundaries
) {
    dat = payload["data"]["intervals"]
    intervals = [{
        x: scInterval(dat['start']-1.25, dat['end']+1.1,
            start_site, buffer, strand, exonBoundaries)['start'],
        y: scInterval(dat['start']+1.25, dat['end']+1.1,
            start_site, buffer, strand, exonBoundaries)['end'],
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