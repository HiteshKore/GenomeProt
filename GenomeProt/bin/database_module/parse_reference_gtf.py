#############GTF FILE PARSER###########################################
# Script By: Hitesh Kore
#######################################################################

class GTFParser():
    def __init__(self, raw_line):
        line = raw_line.strip()
        cols = [col.strip() for col in line.split('\t')]

        self.valid = not line.startswith('#') and (len(cols) >= 9)
        if not self.valid:
            return

        chrom = cols[0]
        self.feature = cols[2]
        start = int(cols[3])
        end = int(cols[4])
        self.strand = cols[6]
        self.coordinates = f"{chrom}:{start}-{end}"

        meta_info_vals = cols[8].split(';')
        def concat_meta_attributes(attribute_name):
            return ','.join([meta_info_val.replace(attribute_name, '').replace('"', '').strip() for meta_info_val in meta_info_vals if attribute_name in meta_info_val])

        self.gene_id = concat_meta_attributes("gene_id")
        self.gene_name = concat_meta_attributes("gene_name")
        self.gene_type = concat_meta_attributes("gene_type")
        self.transcript_id = concat_meta_attributes("transcript_id")
        self.transcript_type = concat_meta_attributes("transcript_type")

    def featureExists(self, feature):
        return feature in self.feature