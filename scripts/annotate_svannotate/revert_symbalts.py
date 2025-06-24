import sys
from pysam import VariantFile

def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: program.py <annotated_vcf> <original_vcf>\n")
        sys.exit(1)

    annotated_vcf_path = sys.argv[1]
    original_vcf_path = sys.argv[2]

    with VariantFile(annotated_vcf_path) as annotated_vcf, \
         VariantFile(original_vcf_path) as original_vcf:

        # Create a dictionary for quick lookups of original records
        original_records = {}
        for record in original_vcf.fetch():
            original_records[record.id] = record

        # Setup output VCF
        vcf_out = VariantFile('-', 'w', header=annotated_vcf.header)

        for annotated_record in annotated_vcf.fetch():
            if annotated_record.id in original_records:
                original_record = original_records[annotated_record.id]
                
                # Create a new record to avoid modifying the original
                new_record = annotated_record.copy()

                # Revert REF, ALT and SVLEN
                new_record.ref = original_record.ref
                new_record.alts = original_record.alts
                if 'SVLEN' in original_record.info:
                    new_record.info['SVLEN'] = original_record.info['SVLEN']
                
                # For BNDs, if the BND_ALT tag was used, clear it
                if 'BND_ALT' in new_record.info:
                    del new_record.info['BND_ALT']
                
                vcf_out.write(new_record)
            else:
                # If for some reason the variant is not in the original VCF, write it as is
                vcf_out.write(annotated_record)

if __name__ == "__main__":
    main() 