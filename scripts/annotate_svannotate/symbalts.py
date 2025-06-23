import sys
import re
from pysam import VariantFile, FastaFile
def main():
    if len(sys.argv)!=2:
        sys.stderr.write("Usage: program.py vcf\n")
        sys.exit(1)
    vcf_in = VariantFile(sys.argv[1])
    vcf_out = VariantFile('-', 'w', header=vcf_in.header)

    vcf_in.header.add_line('##INFO=<ID=BND_ALT,Number=1,Type=String,Description="BND info from ALT field">')
    vcf_out.header.add_line('##INFO=<ID=BND_ALT,Number=1,Type=String,Description="BND info from ALT field">')

    # this prevents a cryptic error about the INFO tag id being too large
    if 'END' not in vcf_in.header.info:
        vcf_in.header.add_line('##INFO=<ID=END,Number=1,Type=Integer>')
        vcf_out.header.add_line('##INFO=<ID=END,Number=1,Type=Integer>')

    for rec in vcf_in.fetch():
        # for BNDs, preserve ALT info in a tag (BND_ALT)
        if rec.info['SVTYPE'] == "BND":
            rec.info["BND_ALT"] = rec.alts[0]

        rec.alts = ("<%s>" % rec.info['SVTYPE'],)
        rec.ref = 'N'

        # setting symbolic alts breaks some END tags, this is a fix:
        if rec.info['SVTYPE'] in ['DEL','DUP','INV']:
            if isinstance(rec.info['SVLEN'], tuple):
                svlen = abs(rec.info['SVLEN'][0])
            else:
                svlen = abs(rec.info['SVLEN'])
            rec.stop = rec.pos+svlen

        vcf_out.write(rec) # write every variant
    
if __name__=="__main__":main()
