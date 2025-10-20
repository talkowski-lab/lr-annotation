import sys
import numpy

def main():
    if len(sys.argv)!=5:
        sys.stderr.write("Usage: program.py PALMER_calls type sample ref_fai\n")
        sys.exit(1)

    callfile = open(sys.argv[1],'r')
    MEtype = sys.argv[2]
    sample = sys.argv[3]
    ref_fai = open(sys.argv[4],'r')

    min_conf = 1

    # write VCF header
    print('##fileformat=VCFv4.2')
    for line in ref_fai:
        line = line.strip().split('\t')
        contig, length = line[0], int(line[1])
        print('##contig=<ID=%s,length=%d>' % (contig, length))
    ref_fai.close()

    print('##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print('##INFO=<ID=CONF_READS,Number=.,Type=Integer,Description="Number of confident supporting reads">')
    print('##INFO=<ID=ORI,Number=.,Type=String,Description="Orientation of insertion">')
    print('##INFO=<ID=ME_TYPE,Number=.,Type=String,Description="Type of mobile element">')
    print('##INFO=<ID=SVTYPE,Number=.,Type=String,Description="Type of SV">')
    print('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length of SV">')
    print('##INFO=<ID=POLYA_LEN,Number=.,Type=Integer,Description="Length of polyA tail">')
    print('##INFO=<ID=TSD_5PRIME_LEN,Number=.,Type=Integer,Description="Length of 5prime TSD">')
    print('##INFO=<ID=TSD_3PRIME_LEN,Number=.,Type=Integer,Description="Length of 3prime TSD">')
    print('##INFO=<ID=TRANSD_LEN,Number=.,Type=Integer,Description="Length of transduction">')
    print('##INFO=<ID=INVERSION_5PRIME,Number=0,Type=Flag,Description="Whether the MEI has a 5 prime inversion">')
    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % sample)

    # Process calls
    calls = {}
    for line in callfile:
        if line.startswith('cluster_id\t'): # header
            continue
        line = line.strip().split('\t')
        cur = {}
        cid, cur["chrom"], cur["conf"], cur["ori"], cur["polyAlen"], cur["5TSDlen"], cur["3TSDlen"], cur["transDlen"], cur['5inv'] = line[0], line[1], line[10], line[13], line[14], line[15], line[16], line[17], int(line[18])

        # Take the median of all possible positions
        cur["pos"] = round(numpy.median([int(line[2]), int(line[3]), int(line[4]), int(line[5])]))

        # Calculate MEI length
        tmp_end = round(numpy.median([int(line[8]), int(line[9])]))
        tmp_start = round(numpy.median([int(line[6]), int(line[7])]))
        if cur['5inv']: # Calculate length for 5' inversion differently
            tmp_5invend = int(line[19])
            tmp_5invstart = int(line[20])
            cur["length"] = tmp_end - tmp_start + int(cur["transDlen"]) + tmp_5invend - tmp_5invstart

        else:
            cur["length"] = tmp_end - tmp_start + int(cur["transDlen"])

        if int(cur["conf"]) < min_conf:
            continue
        
        if cid not in calls:
            calls[cid] = cur
        else:
            calls[cid]["conf"]      = ','.join(list(set([cur["conf"], calls[cid]["conf"]])))
            calls[cid]["ori"]       = ','.join(list(set([cur["ori"], calls[cid]["ori"]])))
            calls[cid]["polyAlen"]  = ','.join(list(set([cur["polyAlen"], calls[cid]["polyAlen"]])))
            calls[cid]["5TSDlen"]   = ','.join(list(set([cur["5TSDlen"], calls[cid]["5TSDlen"]])))
            calls[cid]["3TSDlen"]   = ','.join(list(set([cur["3TSDlen"], calls[cid]["3TSDlen"]])))
            calls[cid]["transDlen"] = ','.join(list(set([cur["transDlen"], calls[cid]["transDlen"]])))
            calls[cid]["5inv"] = cur['5inv']

    for cid in calls:
        cur = calls[cid]
        if cur['5inv']:
            print('%s\t%s\t%s\tN\tN\t60\t.\tSVTYPE=INS;ME_TYPE=%s;SVLEN=%d;CONF_READS=%s;ORI=%s;POLYA_LEN=%s;TSD_5PRIME_LEN=%s;TSD_3PRIME_LEN=%s;TRANSD_LEN=%s;INVERSION_5PRIME\tGT\t0/1' % (cur["chrom"], cur["pos"], cid,  MEtype, cur["length"], cur["conf"], cur["ori"], cur["polyAlen"], cur["5TSDlen"],   cur["3TSDlen"], cur["transDlen"]))
        else:
            print('%s\t%s\t%s\tN\tN\t60\t.\tSVTYPE=INS;ME_TYPE=%s;SVLEN=%d;CONF_READS=%s;ORI=%s;POLYA_LEN=%s;TSD_5PRIME_LEN=%s;TSD_3PRIME_LEN=%s;TRANSD_LEN=%s\tGT\t0/1' % (cur["chrom"], cur["pos"], cid, MEtype, cur["length"], cur["conf"], cur["ori"], cur["polyAlen"], cur["5TSDlen"], cur["3TSDlen"], cur["transDlen"]))


if __name__=="__main__":
    main()
