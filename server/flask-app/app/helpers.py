"""
A set of sqlite3 helper functions
"""
import json
import requests

# import the datasets
from . import variant_db
from . import features_db


def get_all_te_values():
    """
    Gets the translational efficiency values for all orfs
    """
    cursor = features_db.get_db()
    query = cursor.execute('SELECT efficiency from orf_features')
    result = [te['efficiency'] for te in query.fetchall() if te is not None]
    features_db.close_db()
    return result


def parse_five_prime_utr_variant_consequence(conseq_str):
    """
    Parses the consequence str into a keyed dictionary as per
    https://github.com/ImperialCardioGenetics/UTRannotator#the-detailed-annotation-for-each-consequence

    """
    return {
        annotation.split(':')[0]: annotation.split(':')[1]
        for annotation in conseq_str.split(',')
    }


def search_enst_by_transcript_id(variant_id):
    """
    Find all of the transcripts associated with a gpos
    """
    gpos = variant_id.split('-')[1]
    cursor = features_db.get_db()

    # Query and search for results
    query = cursor.execute(
        'SELECT ensembl_transcript_id FROM genome_to_transcript_coordinates WHERE genomic_pos=? ',  # noqa: E501 # pylint: disable=C0301
        [gpos],
    )
    result = query.fetchall()
    features_db.close_db()
    return result


def parse_values(val, start_site, buffer_length):
    """
    Converts the val into an int
    """
    if val == 'NA':  # pylint: disable=R1705
        return start_site + buffer_length
    else:
        return int(val)


def convert_uploaded_variation_to_variant_id(uploaded_variation):
    """
    Replaces the uploaded variation in VEP to a gnomad-esq variant id
    """
    return uploaded_variation.replace('_', '-').replace('/', '-')


def get_utr_annotation_for_list_variants(
    list_variants, possible_variants_dict, start_site, buffer_length
):
    """
    Get the utr annotation for a list of variants
    """

    high_impact_utr_variants = list(
        set(  # pylint: disable=R1718
            [
                v['variant_id']
                for v in possible_variants_dict
                if v['variant_id'] in list_variants
            ]
        )
    )

    if len(high_impact_utr_variants) > 0:
        return [
            find_intervals_for_utr_consequence(
                var_id=v['variant_id'],
                conseq_type=v['five_prime_UTR_variant_consequence'],
                conseq_dict=v['five_prime_UTR_variant_annotation'],
                cdna_pos=v['cDNA_position'],
                start_site=start_site,
                buffer_length=buffer_length,
            )
            for v in possible_variants_dict
            if v['variant_id'] in high_impact_utr_variants
        ]
    return []


def find_intervals_for_utr_consequence(
    var_id, conseq_type, conseq_dict, cdna_pos, start_site, buffer_length
):
    """
    Parses the output of UTR annotator to a dictionary of
    intervals [start, end] for the visualization
    """
    intervals = {}
    intervals['variant_id'] = var_id
    conseq_dict = parse_five_prime_utr_variant_consequence(conseq_dict)
    if conseq_type == 'uAUG_gained':
        # Done
        intervals['start'] = cdna_pos.split('-')[0]
        intervals['end'] = int(cdna_pos.split('-')[0]) + parse_values(
            conseq_dict['uAUG_gained_DistanceToStop'], start_site, buffer_length
        )
        intervals['viz_type'] = 'New Feature'
        intervals['viz_color'] = 'main'
        intervals['type'] = 'uAUG_gained'
        intervals['kozak_strength'] = conseq_dict['uAUG_gained_KozakStrength']

    elif conseq_type == 'uAUG_lost':
        # Done
        intervals['start'] = int(conseq_dict['uAUG_lost_CapDistanceToStart'])

        intervals['end'] = int(start_site) - parse_values(
            conseq_dict['uAUG_lost_DistanceToCDS'], start_site, buffer_length
        )
        intervals['viz_type'] = 'New Feature'
        intervals['viz_color'] = 'null'
        intervals['type'] = 'uAUG_lost'
        intervals['kozak_strength'] = conseq_dict['uAUG_lost_KozakStrength']

    elif conseq_type == 'uSTOP_lost':
        intervals['start'] = cdna_pos
        intervals['end'] = cdna_pos
        intervals['viz_type'] = 'New Feature'
        intervals['viz_color'] = 'main'
        intervals['type'] = 'uSTOP_lost'
        intervals['kozak_strength'] = conseq_dict['uSTOP_lost_KozakStrength']

    elif conseq_type == 'uSTOP_gained':
        # Get the cdna position of the start site
        intervals['start'] = start_site - int(
            conseq_dict['uSTOP_gained_ref_StartDistanceToCDS']
        )
        # cDNA position of the new stop gained
        intervals['end'] = start_site - parse_values(
            conseq_dict['uSTOP_gained_newSTOPDistanceToCDS'], start_site, buffer_length
        )
        intervals['viz_type'] = 'New Feature'
        intervals['viz_color'] = 'main'
        intervals['type'] = 'uSTOP_gained'
        intervals['kozak_strength'] = conseq_dict['uSTOP_gained_KozakStrength']

    # Once we have indels as well
    elif conseq_type == 'uFrameShift':
        frame_shift_pos = int(start_site) - parse_values(
            int(conseq_dict['uFrameShift_ref_StartDistanceToCDS']),
            start_site, buffer_length
        )

        intervals['start'] = frame_shift_pos
        intervals['end'] = frame_shift_pos + parse_values(
            conseq_dict['uFrameShift_alt_type_length'], start_site, buffer_length
        )
        print("SAJIDSODJASIDASJODISAJASOIDJSAIODSHELLLLLLLLOOOOOOOO_____________")
        print(intervals['start'])
        print(intervals['end'])
        intervals['viz_type'] = 'New Feature'
        intervals['viz_color'] = 'main'
        intervals['type'] = 'uFrameshift'
        intervals['kozak_strength'] = conseq_dict['uFrameShift_KozakStrength']


    intervals.update(conseq_dict)

    return intervals


def convert_between_ids(from_id, from_entity, to_entity):
    """
    Converts between different entity
    """
    list_possible_cols = [
        'ensembl_transcript_id',
        'ensembl_gene_id',
        'ensembl_protein_id',
        'ncbi_gene_id',
        'refseq_transcript_id',
        'refseq_protein_id',
        'mane_status',
        'name',
        'hgnc_symbol',
        'hgnc_id',
    ]
    if from_entity in list_possible_cols and to_entity in list_possible_cols:
        rows = features_db.get_db()
        query = f"SELECT {to_entity} FROM mane_summary WHERE {from_entity}='{from_id}'"
        results = rows.execute(query).fetchone()
        features_db.close_db()
        if results is not None:
            return results[to_entity]
    return None


def get_genomic_features(ensg):
    """
    Gets the genomic features tab
    """
    cursor = features_db.get_db()

    # Query and search for results
    query = cursor.execute(
        'SELECT * FROM mane_genomic_features WHERE ensembl_gene_id=? ', [ensg]
    )
    result = query.fetchall()
    features_db.close_db()
    return result


def find_all_high_impact_utr_variants(ensembl_transcript_id):
    """
    Finds all possible UTR variants for a
     given transcript id from the database
    """
    db = variant_db.get_db()
    cursor = db.execute(
        'SELECT variant_id FROM variant_annotations WHERE ensembl_transcript_id=?',
        [ensembl_transcript_id],
    )
    rows = cursor.fetchall()
    variant_db.close_db()
    return [i[0] for i in rows]


def get_transcript_position(ensembl_transcript_id, gpos):
    """
    Gets the transcript position for the transcript / gpos combo
    """
    db = features_db.get_db()
    cursor = db.execute(
        'SELECT transcript_pos FROM genome_to_transcript_coordinates WHERE ensembl_transcript_id=? AND genomic_pos=?',  # noqa: E501 # pylint: disable=C0301
        [ensembl_transcript_id, int(gpos)],
    )
    result = cursor.fetchone()
    features_db.close_db()
    return result['transcript_pos']


def get_possible_variants(ensembl_transcript_id):
    """
    Searches the database for variants
    """
    var_db = variant_db.get_db()
    cursor = var_db.execute(
        'SELECT annotations FROM variant_annotations WHERE ensembl_transcript_id =?',
        [ensembl_transcript_id],
    )
    rows = cursor.fetchall()
    variants = [json.loads(row[0]) for row in rows]
    variant_db.close_db()

    for v in variants:
        v.update(
            {
                'variant_id': convert_uploaded_variation_to_variant_id(
                    v['#Uploaded_variation']
                )
            }
        )

    return variants


def process_gnomad_data(gnomad_data, ensembl_transcript_id):
    """
    Get the gnomAD data and find their transcript coordinates
    and filter to 5' UTR variants
    """
    # Filtering to SNVs for now and those with a 5' UTR consequence
    gnomad_data['clinvar_variants'] = [
        clinvar
        for clinvar in gnomad_data['clinvar_variants']
        if clinvar['major_consequence'] == '5_prime_UTR_variant'
        and len(clinvar['ref']) == 1
        and len(clinvar['alt']) == 1
    ]
    gnomad_data['variants'] = [
        var
        for var in gnomad_data['variants']
        if var['transcript_consequence']['major_consequence'] == '5_prime_UTR_variant'
        and len(var['ref']) == 1
        and len(var['alt']) == 1
    ]

    # Add the transcript relative positions for both
    for clinvar in gnomad_data['clinvar_variants']:
        clinvar.update(
            {'tpos': get_transcript_position(ensembl_transcript_id, clinvar['pos'])}
        )
    for var in gnomad_data['variants']:
        var.update({'tpos': get_transcript_position(ensembl_transcript_id, var['pos'])})

    # Get the variant ids
    gnomad_variants_list = [var['variant_id'] for var in gnomad_data['variants']]

    # Get the list of clinvar variants
    clinvar_variants_list = [
        var['variant_id'] for var in gnomad_data['clinvar_variants']
    ]

    return gnomad_data, gnomad_variants_list, clinvar_variants_list


def get_transcript_features(ensembl_transcript_id):
    """
    Get transcript features
    """
    db = features_db.get_db()
    cursor = db.execute(
        'SELECT * FROM mane_transcript_features WHERE ensembl_transcript_id=?',
        [ensembl_transcript_id],
    )
    rows = cursor.fetchone()
    features_db.close_db()
    return rows


def get_genome_to_transcript_intervals(ensembl_transcript_id, tpos):
    """
    Retrieves all of the features of the uorfs / uorfs
    for the native architechure of the gene
    """
    print(tpos)
    db = features_db.get_db()
    cursor = db.execute(
        'SELECT genomic_pos FROM genome_to_transcript_coordinates WHERE ensembl_transcript_id=? AND transcript_pos=?',  # noqa: E501 # pylint: disable=C0301
        [ensembl_transcript_id, tpos],
    )
    result = cursor.fetchone()
    features_db.close_db()
    return result['genomic_pos']



def get_all_orfs_features(ensembl_transcript_id):
    """
    Retrieves all of the features of the uorfs / uorfs
    for the native architechure of the gene
    """
    db = features_db.get_db()
    cursor = db.execute(
        'SELECT * FROM orf_features WHERE ensembl_transcript_id=?',
        [ensembl_transcript_id],
    )
    rows = cursor.fetchall()
    features_db.close_db()

    return rows


def get_constraint_score(ensembl_gene_id):
    """
    Get constraint score from the features db
    @param ensembl_gene_id
    @returns constraint score (double)
    """
    db = features_db.get_db()
    cursor = db.execute(
        'SELECT loeuf FROM loeuf_constraint WHERE ensembl_gene_id=?',
        [ensembl_gene_id],
    )
    result = cursor.fetchone()
    features_db.close_db()
    return result['loeuf']


def get_gnomad_variants_in_utr_regions(utr_regions):
    """
    gnomAD search in utr regions
    """
    searches = [
        gnomad_api_search_by_region(
            chrom=ur['chr'][3:], start=ur['start'], stop=ur['end']
        )['region']
        for ur in utr_regions
    ]
    data = {}
    # Append
    data['variants'] = sum([s['variants'] for s in searches], [])
    data['clinvar_variants'] = sum([s['clinvar_variants'] for s in searches], [])
    return data


def gnomad_api_search_by_region(chrom, start, stop):
    """
    For prototyping purposes
    """
    region_variant_query = """
    query get_data ($chrom : String!,
                    $start : Int!,
                    $stop : Int!){
    region(chrom: $chrom, start:$start, stop:$stop, reference_genome:GRCh38){
    clinvar_variants {
            transcript_id
            ref
            pos
            alt
            in_gnomad
            clinvar_variation_id
            gold_stars
            variant_id
            review_status
            hgvsc
            clinical_significance
            major_consequence
        }
        variants(dataset: gnomad_r3) {
            ref
            pos
            alt
            hgvsc
            variant_id
            genome {
            af
            an
            ac
            }
            transcript_consequence {
            is_mane_select
            major_consequence
            sift_prediction
            polyphen_prediction
            is_mane_select_version
            }
        }
    }
  }
  """
    response = requests.post(
        'https://gnomad.broadinstitute.org/api',
        data=json.dumps(
            {
                'query': region_variant_query,
                'variables': {'start': start, 'stop': stop, 'chrom': chrom},
            }
        ),
        headers={
            'Content-Type': 'application/json',
        },
    ).json()

    return response['data']
