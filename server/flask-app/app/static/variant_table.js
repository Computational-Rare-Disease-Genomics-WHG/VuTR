create_variant_table = function (div, data){
    $(div).jsGrid({
        height: "100%",
        width: "100%",
        sorting: true,
        pageSize: 6,
        paging: true,
        data: data,
        fields: [{
                 name: "clinvar_variation_id",
                 type: "text",
                 width: 150,
                 title: "ClinVar Variation ID"
           },
           {
                 name: "ref",
                 type: "text",
                 title: "REF"
           },
           {
                 name: "alt",
                 type: "text",
                 title: "ALT"
           },
           {
                 name: "pos",
                 type: "text",
                 title: "POS"
           },
           {
                 name: "major_consequence",
                 type: "text",
                 title: "Major Consequence"
           },
           {
                 name: "in_gnomad",
                 type: "checkbox",
                 title: "In gnomAD?"
           },
           {
                 name: "clinical_significance",
                 type: "text",
                 title: "Clincal Significance"
           },
           {
                 name: "hgvsc",
                 type: "text",
                 title: "HGVSC"
           },
        ]
     });
}
