import sys


def main():
    if len(sys.argv) != 6:
        sys.stderr.write(
            "Usage: program.py PALMER_calls PALMER_read_seqs type sample ref_fai\n"
        )
        sys.exit(1)

    callfile = open(sys.argv[1], "r")
    seqfile = open(sys.argv[2], "r")
    MEtype = sys.argv[3]
    sample = sys.argv[4]
    ref_fai = open(sys.argv[5], "r")

    min_conf = 1

    # process SEQ info
    seqs_all = {}
    for line in seqfile:
        line = line.strip().split("\t")
        ID, seq = line[0], line[6]
        if ID not in seqs_all:
            seqs_all[ID] = [seq]
        else:
            seqs_all[ID].append(seq)
    seqfile.close()

    # get median-length seq for each cluster
    seq = {}
    for ID in seqs_all:
        seqs_all[ID].sort(key=lambda s: len(s))  # sort by length of seq
        seq[ID] = seqs_all[ID][
            int(len(seqs_all[ID]) / 2)
        ]  # this retrieves the median value from the sorted list

    # write VCF header
    print("##fileformat=VCFv4.2")
    for line in ref_fai:
        line = line.strip().split("\t")
        contig, length = line[0], int(line[1])
        print("##contig=<ID=%s,length=%d>" % (contig, length))
    ref_fai.close()

    print(
        '##ALT=<ID=INS,Description="Insertion of novel sequence '
        'relative to the reference">'
    )
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print(
        "##INFO=<ID=CONF_READS,Number=.,Type=Integer,"
        'Description="Number of confident supporting reads">'
    )
    print('##INFO=<ID=ORI,Number=.,Type=String,Description="Orientation of insertion">')
    print(
        '##INFO=<ID=ME_TYPE,Number=.,Type=String,Description="Type of mobile element">'
    )
    print('##INFO=<ID=SVTYPE,Number=.,Type=String,Description="Type of SV">')
    print('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length of SV">')
    print(
        '##INFO=<ID=POLYA_LEN,Number=.,Type=Integer,Description="Length of polyA tail">'
    )
    print(
        "##INFO=<ID=TSD_5PRIME_LEN,Number=.,Type=Integer,"
        'Description="Length of 5prime TSD">'
    )
    print(
        "##INFO=<ID=TSD_3PRIME_LEN,Number=.,Type=Integer,"
        'Description="Length of 3prime TSD">'
    )
    print(
        "##INFO=<ID=TRANSD_LEN,Number=.,Type=Integer,"
        'Description="Length of transduction">'
    )
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % sample)

    # process PALMER calls
    calls = {}
    for line in callfile:
        if line.startswith("cluster_id\t"):  # header
            continue
        line = line.strip().split("\t")
        cur = {}
        (
            cid,
            cur["chrom"],
            cur["start"],
            cur["conf"],
            cur["ori"],
            cur["polyAlen"],
            cur["5TSDlen"],
            cur["3TSDlen"],
            cur["transDlen"],
        ) = (
            line[0],
            line[1],
            line[2],
            line[10],
            line[13],
            line[14],
            line[15],
            line[16],
            line[17],
        )
        if int(cur["conf"]) < min_conf:
            continue
        if cid not in calls:
            calls[cid] = cur
        else:
            calls[cid]["conf"] = ",".join(list(set([cur["conf"], calls[cid]["conf"]])))
            calls[cid]["ori"] = ",".join(list(set([cur["ori"], calls[cid]["ori"]])))
            calls[cid]["polyAlen"] = ",".join(
                list(set([cur["polyAlen"], calls[cid]["polyAlen"]]))
            )
            calls[cid]["5TSDlen"] = ",".join(
                list(set([cur["5TSDlen"], calls[cid]["5TSDlen"]]))
            )
            calls[cid]["3TSDlen"] = ",".join(
                list(set([cur["3TSDlen"], calls[cid]["3TSDlen"]]))
            )
            calls[cid]["transDlen"] = ",".join(
                list(set([cur["transDlen"], calls[cid]["transDlen"]]))
            )

        # add the sequence
        if cid not in seq:  # this shouldn't happen if using the matching input files
            sys.stderr.write("Error: cluster '%s' not found in read file\n" % cid)
            sys.exit()
        calls[cid]["seq"] = seq[cid]

    for cid in calls:
        cur = calls[cid]
        print(
            "%s\t%s\t%s\tN\t%s\t60\t."
            "SVTYPE=INS;ME_TYPE=%s;SVLEN=%d;CONF_READS=%s;ORI=%s;"
            "POLYA_LEN=%s;TSD_5PRIME_LEN=%s;TSD_3PRIME_LEN=%s;"
            "TRANSD_LEN=%s\tGT\t0/1"
            % (
                cur["chrom"],
                cur["start"],
                cid,
                cur["seq"],
                MEtype,
                len(cur["seq"]),
                cur["conf"],
                cur["ori"],
                cur["polyAlen"],
                cur["5TSDlen"],
                cur["3TSDlen"],
                cur["transDlen"],
            )
        )


if __name__ == "__main__":
    main()
