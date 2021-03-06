"""
download_sorf_query.py
Download Sorfs ribo-profiling data through the SOAP API
"""
import io
import pandas as pd
from suds.client import Client

client = Client('http://biomart.biobix.be/martsoap?wsdl')
QUERY = """<!DOCTYPE Query>
<Query client="true" processor="TSV" limit="-1" header="1">
    <Dataset name="BioMart" config="Human">
        <Filter name="human__annotation_104" value="5UTR" filter_list=""/>
            <Attribute name="human__sorf_id_104"/>
            <Attribute name="human__chr_104"/>
            <Attribute name="human__sorf_end_104"/>
            <Attribute name="human__strand_104"/>
            <Attribute name="human__sorf_begin_104"/>
            <Attribute name="human__upstream_gene_distance_104"/>
            <Attribute name="human__downstream_gene_distance_104"/>
            <Attribute name="human__tr_seq_104"/>
            <Attribute name="human__id_104"/>
            <Attribute name="human__start_codon_104"/>
    </Dataset>
</Query>
"""
print('Querying biomart.biobix.be')
output = client.service.getResults(QUERY)
print('Reading into memory')
df = pd.read_csv(io.StringIO(output), sep='\t')
print('Writing to file')
df.to_csv('../../data/pipeline/SORFS/sorfs.tsv', sep='\t')
print('Completed SORFS.org download')
