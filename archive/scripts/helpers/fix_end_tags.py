import sys
# pysam has a bug where it doesn't write END consistently for each SVTYPE, so we need to do it manually

def main():
    if len(sys.argv)!=2:
        sys.stderr.write("Usage: program.py vcf\n")
        sys.exit(1)
    vcf_in = open(sys.argv[1],'r')

    for line in vcf_in:
        if line.startswith('#'):
            sys.stdout.write(line)
            continue
        line = line.strip().split('\t')

        #parse INFO field
        infos = line[7].split(';')
        has_end = False
        SVLEN, SVTYPE = None, None
        for info in infos:
            if '=' not in info:
                continue
            tag, val = info.split('=')
            if tag=="END":
                has_end = True
            elif tag=="SVLEN":
                SVLEN = abs(int(val))
            elif tag=="SVTYPE":
                SVTYPE = val

        if not has_end:
            pos = int(line[1])
            if SVTYPE in ['INS','BND']:
                end=pos
                infos.append("END=%d" % end )
            elif SVTYPE in ['DUP','DEL','INV']:
                if not SVLEN:
                    sys.stderr.write("SVLEN missing in line:\n%s\n" % '\t'.join(line))
                    sys.exit()
                end = pos+SVLEN
                infos.append("END=%d" % end )

        # write up to sample data
        sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (line[0], line[1], line[2], line[3], line[4], line[5], line[6], ';'.join(infos), line[8]))
        # add sample data
        for i in range(9,len(line)):
            sys.stdout.write("\t%s" % line[i])
        sys.stdout.write("\n")

if __name__=="__main__":main()
