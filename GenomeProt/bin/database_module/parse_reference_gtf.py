#############GTF FILE PARSER###########################################
# Script By: Hitesh Kore
#######################################################################

class GTFParser():
    def __init__(self, raw_line):
        line = raw_line.strip()
        cols = [col.strip() for col in line.split('\t')]
        if len(cols) <= 2:
            return

        self.chr = cols[0]
        self.feature = cols[2]
        self.start = int(cols[3])
        self.end = int(cols[4])
        self.strand = cols[6]
        self.coord_range = f"{self.start}-{self.end}"
        self.coordinates = f"{self.chr}:{self.coord_range}"
        self.length = self.end - self.start + 1

        meta = cols[8].split(';')
        self.ensgene = ','.join([i.replace("gene_id", '').replace('"', '') for i in meta if "gene_id" in i]).strip()
        self.transcriptid = ','.join([i.replace("transcript_id", '').replace('"', '') for i in meta if "transcript_id" in i]).strip()
        self.genename = ','.join([i.replace("gene_name", '').replace('"', '').strip() for i in meta if "gene_name" in i]).strip()
        self.genetype = ','.join([i.replace("gene_type", '').replace('"', '').strip() for i in meta if "gene_type" in i]).strip()
        self.transcript_type = ','.join([i.replace("transcript_type", '').replace('"', '').strip() for i in meta if "transcript_type" in i]).strip()

    def featureExists(self, feature):
        return feature in self.feature