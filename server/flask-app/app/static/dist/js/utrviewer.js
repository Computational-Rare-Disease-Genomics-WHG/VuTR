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
    strand
) {
    if (strand == '+') {
        return ({
            'start': start,
            'end': Math.min(end - 1, start_site + buffer + 1) + 1
        });
    } else {
        return ({
            "start": Math.max((start_site + buffer + 1) - end, 1),
            "end": (start_site + buffer + 1) - start,
        });
    }
}


/**
 * Mutates a sequence with the variant and returns the mutated sequence
 * @function mutateSeq
 * @param {string} seq - The sequence to be mutated
 * @param {number} pos - The position of the mutation
 * @param {string} ref - The reference allele
 * @param {string} alt - The alternative allele
 * @returns {string} - The mutated sequence
 */
function mutateSeq(seq, pos, ref, alt) {
    return seq.slice(0, pos - 1) + alt + seq.slice(pos + ref.length - 1);
}










/**
 * 
 * This function converts translation effeciency to a color based on a linear scale given 
 * by c0 and c1 and a fraction f
 * @function interpolateColor
 * @param {str} c0 beginning color in hex
 * @param {str} c1 end color in hex
 * @param {float} f - the fraction of the color to interpolate
 * @returns Hex color without the #
 */
function interpolateColor(c0, c1, f){
    c0 = c0.match(/.{1,2}/g).map((oct)=>parseInt(oct, 16) * (1-f))
    c1 = c1.match(/.{1,2}/g).map((oct)=>parseInt(oct, 16) * f)
    let ci = [0,1,2].map(i => Math.min(Math.round(c0[i]+c1[i]), 255))
    return ci.reduce((a,v) => ((a << 8) + v), 0).toString(16).padStart(6, "0")
}


/**
Converts an int to a color in hex based on a linear scale from 0 to 150
@function intToColor
@param {number} value - The value to convert to color
@returns {string} - The color in hex
*/

function intToColor(value) {
    value = Math.min(Math.max(value, 0), 150)/150;
    return '#' + interpolateColor("ffffff", "007bff", value)
}


/** 
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
        exon_start = scInterval(new_x, new_x +
            width_exon, start_site, buffer, strand)['start']
        exon_end = scInterval(new_x, new_x +
            width_exon, start_site, buffer, strand)['end'];

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

    if (data.length == 1){
        data = data[0]
    }

    $('#modal-container').empty();

    /// Type is the data_type of the detail presented

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

    var title = {
        'gnomad': 'gnomAD Variant',
        'gnomad_multiple_variants': 'gnomAD Variant',
        'clinvar': 'ClinVar Variants',
        'clinvar_multiple_variants': 'ClinVar Variants',
        'orf': 'ORF detail',
        'smorf': 'smORF',
        'user-supplied': 'Variant'
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
      <h5 class="modal-title" id="exampleModalLongTitle">${title[type]}</h5>
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
    if (type == 'gnomad_multiple_variants') {
        let tabsHtml = '';
        let tabContentHtml = '';


        /* Create Tabs and Modals */
        data.forEach((variant, index) => {
            const isActive = index === 0 ? 'active' : '';
            const tabId = `variant-tab-${index}`;
            const isShown = index === 0 ? 'show' : '';
            const tabContentId = `variant-tab-content-${index}`;

            // Generate the tab headers
            tabsHtml += `
              <a class="nav-link ${isActive}" id="${tabId}-tab" data-toggle="pill" href="#${tabContentId}" role="tab" aria-controls="${tabContentId}" aria-selected="${index === 0}">${variant['variant_id']}</a>
            `;

            // Generate the tab content
            tabContentHtml += `
              <div class="tab-pane fade ${isActive} ${isShown}" id="${tabContentId}" role="tabpanel" aria-labelledby="${tabId}-tab">
                <!-- Insert your variant details here using variant.data -->
                <ul>
                  <li><b>REF</b> : ${variant.ref}</li>
                  <li><b>ALT</b> : ${variant.alt}</li>
                  <li><b>Genome Position</b> : ${variant.pos}</li>
                    <li><b>Transcript Position</b> : ${variant.tpos}</li>
                    <li><b>Variant ID</b> : ${variant.variant_id}</li>
                    <li><b>HGVSC</b> : ${variant.hgvsc}</li>
                </ul>
                <hr>
                ${createUtrAnnotationList(variant)}
                <h5>gnomAD Frequency</h5>
                <ul>
                    <li><b>Allele Count</b> : ${variant.genome.ac}</li>
                    <li><b>Allele Number</b> : ${variant.genome.an}</li>
                    <li><b>Allele Frequency (all pop)</b> : ${variant.genome.af}</li>
                </ul>
                <hr>
                <b>View variant in gnomAD</b>: <a href='https://gnomad.broadinstitute.org/variant/${variant['variant_id']}?dataset=gnomad_r3'>${variant['variant_id']}</a>
              </div>
            `;
        });

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
                    <h5 class="modal-title" id="exampleModalLongTitle">${title[type]}</h5>
                    <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                    <span aria-hidden="true" >&times;</span>
                    </button>
                </div>
                <div class="modal-body">
                    <!-- Add your tab navigation here -->
                    <h6>There are multiple variants at this position</h6>
                    <div class="row">
                    <div class="col-5 col-sm-3">
                        <div class="nav flex-column nav-tabs h-100" id="variant-tabs" role="tablist" aria-orientation="vertical">
                        ${tabsHtml}
                        </div>
                    </div>
                    <div class="col-7 col-sm-9">
                        <!-- Add your tab content here -->
                        <div class="tab-content" id="variant-tab-content">
                        ${tabContentHtml}
                    </div>
                    </div>
                    </div>
                </div>
            </div>`;

    }

    if (type == 'user-supplied') {
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

    if (type == 'clinvar_multiple_variants') {
        let tabsHtml = '';
        let tabContentHtml = '';

        /* Create Tabs and Modals */
        data.forEach((variant, index) => {
            const isActive = index === 0 ? 'active' : '';
            const tabId = `variant-tab-${index}`;
            const isShown = index === 0 ? 'show' : '';
            const tabContentId = `variant-tab-content-${index}`;

            // Generate the tab headers
            tabsHtml += `
                <a class="nav-link ${isActive}" id="${tabId}-tab" data-toggle="pill" href="#${tabContentId}" role="tab" aria-controls="${tabContentId}" aria-selected="${index === 0}">${variant['variant_id']}</a>
            `;
            // Generate the tab content
            tabContentHtml += `
                <div class="tab-pane fade ${isActive} ${isShown}" id="${tabContentId}" role="tabpanel" aria-labelledby="${tabId}-tab">
                    <!-- Insert your variant details here using variant.data -->
                    <ul>
                        <li><b>REF</b> : ${variant.ref}</li>
                        <li><b>ALT</b> : ${variant.alt}</li>
                        <li><b>Genome Position</b> : ${variant.pos}</li>
                        <li><b>Transcript Position</b> : ${variant.tpos}</li>
                        <li><b>Variant ID</b> : ${variant.variant_id}</li>
                        <li><b>HGVSC</b> : ${variant.hgvsc}</li>
                    </ul>
                    <hr>
                    <h5>ClinVar variant details</h5>
                    <ul>
                        <li><b>Clinical Significance</b> : ${variant.clinical_significance}</li>
                        <li><b>Review Status</b> : ${variant.review_status}</li>
                        <li><b>Gold Stars</b> : ${variant.gold_stars}</li>
                        <li><b>ClinVar Variation ID</b> : ${variant.clinvar_variation_id}</li>
                        <li><b>Variant in gnomAD?</b> : ${variant.in_gnomad}</li>
                    </ul>
                    <hr>
                    ${createUtrAnnotationList(variant)}
                    <b>View variant in ClinVar</b>: <a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/${variant['clinvar_variation_id']}'>${variant['clinvar_variation_id']}</a>
                </div>
            `;
        });
            
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
                    <h5 class="modal-title" id="exampleModalLongTitle">${title[type]}</h5>
                    <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                    <span aria-hidden="true" >&times;</span>
                    </button>
                </div>
                <div class="modal-body">
                    <!-- Add your tab navigation here -->
                    <h6>There are multiple variants at this position</h6>
                    <div class="row">
                    <div class="col-5 col-sm-3">
                        <div class="nav flex-column nav-tabs h-100" id="variant-tabs" role="tablist" aria-orientation="vertical">
                        ${tabsHtml}
                        </div>
                    </div>
                    <div class="col-7 col-sm-9">
                        <!-- Add your tab content here -->
                        <div class="tab-content" id="variant-tab-content">
                        ${tabContentHtml}
                    </div>
                    </div>
                    </div>
                </div>
            </div>`;
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
                <li><b>ORF Sequence</b> : <p style='word-wrap: break-word'>${data['orf_seq'].substring(0,data['orf_seq'].length-3)}</p></li>
                <li><b>ORF Stop Codon</b> : ${data['orf_seq'].slice(-3)} </li>
                <li><b>ORF Type</b> : ${data['orf_type']} </li>
                <li><b>ORF Frame</b> : ${data['frame']} </li>
              </ul>
              <hr>


              <h5>Translation details</h5>
              <ul>
                <li><b>Translational Efficiency</b> : ${data['efficiency']}</li>
                <li><b>Kozak Consensus Strength</b> : ${data['kozak_consensus_strength']} </li>
                <li><b>Kozak Consensus Context</b> :  ${data['kozak_context']}</li>
              </ul>
              <hr>
              <h5>Translational efficiency of MANE CDS Start Sites</h5>
              
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
                <li><b>Source</b> : ${data['smorf_datasets']}</li>
                <li><b>Type</b> : ${data['type']}</li>
                <li><b>Length</b> : ${data['smorf_length']} aa</li>
                <li><b>Start Codon</b> : ${data['start_codon']} </li>
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
    if ('efficiency' in data && type != 'smorf') {
        var te_chartdiv = document.createElement('div');
        ul.appendChild(te_chartdiv);
        te_chartdiv.innerHTML +=
            `<canvas id = 'te-plot'></canvas`;
        const ctx = document.getElementById('te-plot').getContext('2d');
        // Histogram data from R
        // b <- fread("MANE_transcripts_v1.0.tsv") 
        // fread("../../translational_efficiency.txt")
        // cdt <- b[, .(context=substr(seq, five_prime_utr_length-5,five_prime_utr_length+5))]
        // cdt %<>% .[nchar(context) == 11]
        // hist_obj <- hist(te[cdt]$efficiency, breaks=20)
        // hist_data <- data.table(breaks = hist_obj$breaks[-length(hist_obj$breaks)], counts = hist_obj$counts)
        // hist_data[, label:= paste0(breaks, "-", breaks+5)]
        // hist_data$label %>% paste0(., collapse='","')
        // hist_data$counts %>% paste0(., collapse=",")

        const labels = ["15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-80", "80-85", "85-90", "90-95", "95-100", "100-105", "105-110", "110-115", "115-120", "120-125", "125-130", "130-135", "135-140", "140-145", "145-150"];
        const midpoints = [17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 92.5, 97.5, 102.5, 107.5, 112.5, 117.5, 122.5, 127.5, 132.5, 137.5, 142.5, 147.5];
        const counts = [2, 8, 16, 38, 75, 82, 106, 154, 219, 411, 866, 1539, 2132, 1972, 1794, 1776, 1745, 1694, 1227, 1320, 831, 431, 198, 82, 21, 11, 1];
        var bg_color = midpoints.map((x) => '#111111');

        // Find index where efficiency falls under 
        for (var i = 0; i < labels.length; ++i) {
            const r = labels[i].split('-');
            if (parseInt(r[0]) <= dat['efficiency'] && dat['efficiency'] <= parseInt(r[1])) {
                bg_color[i] = '#E83A5A';
                break;
            }
        }

        const data = {
            labels: labels,
            datasets: [{
                label: 'Translational effienciency of all MANE CDS starts',
                data: counts,
                borderWidth: 0,
                backgroundColor: bg_color,
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


/***
 * Create the first track of the Feature Viewer for gnomAD variants
 * and clinvar variants
 * @param {[Obj]} variants - The list of variants
 * @param {string} variant_track - The type of variant track (gnomad or clinvar)
 * @returns {[Obj]} - The list of variants with the track_id and trackDescription
 */
var createTrackVariants = function(variants, variant_track) {
    const track_vars = variants.reduce((accumulator, variant) => {
        const trackKey = `${variant.tpos}-${variant.ref}`;
        const existingVariant = accumulator.find(item => item.track_id === trackKey);

        if (existingVariant) {
            // If a variant with the same tpos / ref pair exists, set trackDescription to 'M'
            existingVariant.trackDescription = existingVariant.trackDescription + ',' + variant.alt.toUpperCase();

            if (variant_track == 'clinvar') {
                pathogenicity_order = Object.keys(pathogenicity_colors);
                // Store the more pathogenic annotation
                if (pathogenicity_order.indexOf(variant['clinical_significance']) < pathogenicity_order.indexOf(existingVariant['pathogenicity'])) {
                    existingVariant['pathogenicity'] = variant['clinical_significance'];
                }
            }

        } else {
            // Otherwise, add a new entry to the accumulator

            tvar = {
                track_id: trackKey,
                ref: variant.ref,
                tpos: variant.tpos,
                trackDescription: variant.alt.toUpperCase(),
            };

            if (variant_track == 'clinvar') {
                tvar['pathogenicity'] = variant['clinical_significance'];
            }

            accumulator.push(tvar);
        }

        return accumulator;
    }, []);
    return track_vars;
}


/**
 * Quick wrapper to function to search through an
 * array of objects for a certain object 
 * @param {[Obj]} data - The list of objects
 * @param {string} id - The identifier to search
 * @param {string} id_var - The key to search for
 * @returns {Obj or [Obj]}  The object (s) found
 * TODO: 
 */
var searchObj = function(data, id, id_var) {
    // If its a variant 
    if (id_var == 'track_id' || id_var == 'variant_id') {
        found = data.filter(e => {
            return (e[id_var] == id)
        })
        return (found)
    }

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
    conservation_data
) {

    // Subset to the first 100 bases following the CDS
    var sequence = strand == "+" ? seq
        .substring(0, start_site + buffer).toUpperCase() :
        reverse(seq.substring(0, start_site + buffer)).toUpperCase();


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
                    buffer, strand)['start'],
                y: scInterval(start_site + 1,
                    start_site + buffer + 2, start_site,
                    buffer, strand)['end'],
                color: '#58565F',
                description: "CDS",
                id: 'cds_rect',
            }, ...getExonStructure(genomic_features, buffer,
                start_site, strand)],
            color: '#2C302E',
            name: 'Gene Structure',
            className: 'exon_struct',
            type: 'rect',
        });
        document.getElementById("fcds_rect").setAttribute("height", "30");
        document.getElementById("fcds_rect").setAttribute("y", "-8");

        if (smorf.length > 0) {
            smorf_data = [];
            smorf.forEach(e => smorf_data.push({
                x: scInterval(
                    e.transcript_start, e.transcript_end,
                    start_site, buffer, strand
                )['start'],
                y: scInterval(
                    e.transcript_start, e.transcript_end,
                    start_site, buffer, strand
                )['end'],
                id: e.smorf_id,
                color: kozak_colors[e.kozak_consensus_strength]

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
                .orf_start_codon, e.orf_stop_codon-1,
                start_site, buffer, strand)[
                'start'],
            y: scInterval(e.orf_start_codon, e.orf_stop_codon-1,
                start_site, buffer, strand)['end'],
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

    // Add event handlers
    ft2.onFeatureSelected(function(d) {
        console.log(d.detail.id);
        
        if (searchObj(smorf, d.detail.id, 'smorf_id')) {
            openModal(searchObj(smorf, d.detail.id, 'smorf_id'), 'smorf');
        } else {
            openModal(searchObj(uorfs, d.detail.id, 'orf_id'), 'orf');
        }
    });

    /*
    Conservation Track
    */

    var phylop = [];
    var cadd = [];

    // Populate conservation tracks
    conservation.forEach(e => {
        const pos = scInterval(e['tpos'], e['tpos'], start_site, buffer, strand)['start'];

        phylop.push({
            x: pos,
            y: e['phylop']
        });

        cadd.push({
            x: pos,
            y: e['phred_cadd']
        })

    });

    ft2.addFeature({
        data: phylop,
        type: "line",
        className: "phylop_line",
        name: "PhyloP Score",
        color: '#99999',
        height: "2"
    });


    ft2.addFeature({
        data: cadd,
        type: "line",
        className: "cadd_line",
        name: "CADD Score",
        color: '#99999',
        height: "2"
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

    var pop_variants = gnomad_data['variants']; // gnomAD variants from the API
    var pop_variants = pop_variants.map(e => ({
        ...e,
        track_id: `${e.tpos}-${e.ref}`
    }));
    var gnomad_track_vars = createTrackVariants(pop_variants, 'gnomad'); // Create the no impact track
    var pop_var_feat_data = []; // Visualization intervals to pass to feature viewer
    var pop_var_tpos = []; // tmp for storing the transcript positions of the variants
    gnomad_track_vars.forEach(element => {
        tpos = element['tpos']
        // Only append to the track if didn't exist
        if (!pop_var_tpos.includes(tpos)) {
            pop_var_tpos.push(tpos);
            strand_corrected_tpos = scInterval(
                element['tpos'],
                element['tpos'],
                start_site,
                buffer, strand)
            pop_var_feat_data.push({
                x: strand_corrected_tpos['start'],
                y: strand_corrected_tpos['end'] + 0.8,
                id: element['track_id'],
                description: element['trackDescription'],
                className: 'gnomad_no_impact_track',
            });
        }


    });
    if (pop_var_feat_data.length > 0) {
        gnomad_variant_ft.addFeature({
            data: pop_var_feat_data,
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
                        strand)['start'],
                    y: scInterval(
                        element['start'], element[
                            'end'], start_site,
                        buffer, strand)['end'],
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
        if (id_sel.startsWith('annotation')) {
            gnomad_utr_conseq = searchObj(gnomad_utr_impact, id_sel, 'annotation_id');
            openModal(Object.assign(searchObj(pop_variants, gnomad_utr_conseq['variant_id'], 'variant_id'), gnomad_utr_conseq),
                'gnomad');
        } else {
            consequences = searchObj(pop_variants, id_sel, 'track_id')
            if (Array.isArray(consequences)) {
                if (consequences.length > 1) {
                    openModal(Object.assign(searchObj(pop_variants, id_sel, 'track_id'), consequences), 'gnomad_multiple_variants');
                } else if (consequences.length == 1) {
                    openModal(Object.assign(searchObj(pop_variants, id_sel, 'track_id'), consequences[0]), 'gnomad');
                }
            }

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
    var clinvar_variants = clinvar_variants.map(e => ({
        ...e,
        track_id: `${e.tpos}-${e.ref}`
    }));
    var clinvar_track_variants = createTrackVariants(clinvar_variants, 'clinvar');
    var clinvar_var_feat_dat = [];
    clinvar_track_variants.forEach(element => {
        clinvar_var_feat_dat.push({
            x: scInterval(element['tpos'], element['tpos'], start_site, buffer, strand)['start'],
            y: scInterval(element['tpos'], element['tpos'], start_site, buffer, strand)['end'] + 0.8,
            color: pathogenicity_colors[element['pathogenicity']],
            id: element['track_id'],
            description: element['trackDescription'],
            className: 'no_impact_clinvar',
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
                        element['start'] - 1.25, element[
                            'end'] + 1.1, start_site,
                        buffer, strand)['start'],
                    y: scInterval(
                        element['start'] - 1.25, element[
                            'end'] + 1.1, start_site,
                        buffer, strand)['end'],
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
        if (id_sel.startsWith('annotation')) {
            clinvar_utr_conseq = searchObj(clinvar_utr_impact,
                id_sel, 'annotation_id')
            openModal(Object.assign(
                    searchObj(clinvar_variants, clinvar_utr_conseq['variant_id'],
                        'variant_id'),
                    clinvar_utr_conseq),
                'clinvar');
        } else {
            clinvar_utr_conseq = searchObj(clinvar_variants,
                id_sel, 'track_id');
            
            if (Array.isArray(clinvar_utr_conseq)) {
                if (clinvar_utr_conseq.length > 1) {
                    openModal(Object.assign(searchObj(clinvar_variants, id_sel, 'track_id'), clinvar_utr_conseq), 'clinvar_multiple_variants');
                } else if (clinvar_utr_conseq.length == 1) {
                openModal(Object.assign(
                        searchObj(clinvar_variants, id_sel,
                            'track_id'),
                        clinvar_utr_conseq[0]),
                    'clinvar');
        }
    

    }}});

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
    strand
) {
    dat = payload["data"]["intervals"]
    intervals = [{
        x: scInterval(dat['start'] - 1.25, dat['end'] + 1.1,
            start_site, buffer, strand)['start'],
        y: scInterval(dat['start'] + 1.25, dat['end'] + 1.1,
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