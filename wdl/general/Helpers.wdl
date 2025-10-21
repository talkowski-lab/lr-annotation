version 1.0

import "Structs.wdl"

task AddGenotypes {
    input {
        File annot_vcf_file
        File vcf_file
        String hail_docker
        String genome_build
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([annot_vcf_file, vcf_file], "GB") 
    Float base_disk_gb = 10.0

    RuntimeAttr default_attr = object {
                                      mem_gb: 4,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    String filename = basename(annot_vcf_file)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(annot_vcf_file, ".vcf.gz") else basename(annot_vcf_file, ".vcf.bgz")
    String combined_vcf_name = "~{prefix}.GT.vcf.bgz"

    command <<<
    cat <<EOF > merge_vcfs.py
    from pyspark.sql import SparkSession
    import hail as hl
    import numpy as np
    import pandas as pd
    import sys
    import ast
    import os
    import json
    import argparse

    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('-i', dest='vcf_file', help='Input VCF file before annotation')
    parser.add_argument('-a', dest='annot_vcf_file', help='Input VCF file after annotation')
    parser.add_argument('-o', dest='combined_vcf_name', help='Output filename')
    parser.add_argument('--cores', dest='cores', help='CPU cores')
    parser.add_argument('--mem', dest='mem', help='Memory')
    parser.add_argument('--build', dest='build', help='Genome build')

    args = parser.parse_args()

    vcf_file = args.vcf_file
    annot_vcf_file = args.annot_vcf_file
    combined_vcf_name = args.combined_vcf_name
    cores = args.cores  # string
    mem = int(np.floor(float(args.mem)))
    build = args.build

    hl.init(min_block_size=128, 
            local=f"local[*]", 
            spark_conf={
                        "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                        "spark.speculation": 'true'
                        }, 
            tmp_dir="tmp", local_tmpdir="tmp",
                        )

    mt = hl.import_vcf(vcf_file, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome=build)
    try:
        # for haploid (e.g. chrY)
        mt = mt.annotate_entries(
            GT = hl.if_else(
                    mt.GT.ploidy == 1, 
                    hl.call(mt.GT[0], mt.GT[0]),
                    mt.GT)
        )
    except:
        pass

    header = hl.get_vcf_metadata(annot_vcf_file) 
    annot_mt = hl.import_vcf(annot_vcf_file, force_bgz=True, array_elements_required=False, call_fields=[], reference_genome=build)

    annot_mt = annot_mt.union_cols(mt)

    hl.export_vcf(dataset=annot_mt, output=combined_vcf_name, metadata=header, tabix=True)
    EOF

    python3 merge_vcfs.py -i ~{vcf_file} -a ~{annot_vcf_file} -o ~{combined_vcf_name} --cores ~{select_first([runtime_override.cpu_cores, default_attr.cpu_cores])} --mem ~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} \
        --build ~{genome_build}
    >>>

    output {
        File combined_vcf_file = combined_vcf_name
        File combined_vcf_idx = combined_vcf_name + ".tbi"
    }

}

task SplitFileWithHeader {
    input {
        File file
        Int shards_per_chunk
        String cohort_prefix
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    command <<<
        cat <<EOF > split_file.py
        import pandas as pd
        import numpy as np
        import os
        import sys

        file = sys.argv[1]
        shards_per_chunk = int(sys.argv[2])

        file_ext = file.split('.')[-1]
        base_filename = os.path.basename(file).split(f".{file_ext}")[0]
        pd.Series([base_filename]).to_csv('base_filename.txt', index=False, header=None)
        df = pd.read_csv(file, sep='\t')
        n_chunks = np.ceil(df.shape[0]/shards_per_chunk)
        i=0
        while i<n_chunks:
            filename = f"{base_filename}.shard_{i}.{file_ext}"
            df.iloc[i*shards_per_chunk:(i+1)*shards_per_chunk, :].to_csv(filename, sep='\t', index=False)
            i+=1
        
        EOF
        python3 split_file.py ~{file} ~{shards_per_chunk} > stdout
    >>>

    String base_filename = read_lines('base_filename.txt')
    output {
        Array[File] chunks = glob("~{base_filename}.shard_*")
    }

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }
}

task SplitSamples {
    input {
        File vcf_file
        Int samples_per_chunk
        String cohort_prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command <<<
        set -eou pipefail
        bcftools query -l ~{vcf_file} > ~{cohort_prefix}_samples.txt

        cat <<EOF > split_samples.py 
        import os
        import sys
        import pandas as pd
        import numpy as np

        cohort_prefix = sys.argv[1]
        samples_per_chunk = int(sys.argv[2])

        samples = sorted(pd.read_csv(f"{cohort_prefix}_samples.txt", header=None)[0].tolist())
        n = samples_per_chunk  # number of samples in each chunk
        chunks = [samples[i * n:(i + 1) * n] for i in range((len(samples) + n - 1) // n )]  
        
        shard_samples = []
        for i, chunk1 in enumerate(chunks):
            for chunk2 in chunks[i+1:]:
                shard_samples.append(chunk1 + chunk2)

        for i, shard in enumerate(shard_samples):
            pd.Series(shard).to_csv(f"{cohort_prefix}_shard_{i}.txt", index=False, header=None)
        EOF

        python3 split_samples.py ~{cohort_prefix} ~{samples_per_chunk}
    >>>

    output {
        Array[File] sample_shard_files = glob("~{cohort_prefix}_shard_*.txt")
    }
}

task SplitFamilies {
    input {
        File ped_uri
        Int families_per_chunk
        String cohort_prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(ped_uri, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command <<<
        set -eou pipefail

        cat <<EOF > split_samples.py 
        import os
        import sys
        import pandas as pd
        import numpy as np

        cohort_prefix = sys.argv[1]
        families_per_chunk = int(sys.argv[2])  # approximation for trios
        ped_uri = sys.argv[3]

        ped = pd.read_csv(ped_uri, sep='\t').iloc[:, :6]
        ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']

        n = families_per_chunk
        families = ped.family_id.unique().tolist()

        chunks = [families[i * n:(i + 1) * n] for i in range((len(families) + n - 1) // n )]  
        sample_chunks = [ped[ped.family_id.isin(fams)].sample_id.to_list() for fams in chunks]

        for i, shard in enumerate(sample_chunks):
            pd.Series(shard).to_csv(f"{cohort_prefix}_shard_{i}.txt", index=False, header=None)
        EOF

        python3 split_samples.py ~{cohort_prefix} ~{families_per_chunk} ~{ped_uri}
    >>>

    output {
        Array[File] family_shard_files = glob("~{cohort_prefix}_shard_*.txt")
    }
}

task MergeResultsPython {
     input {
        Array[String] tsvs
        String hail_docker
        String merged_filename
        Float input_size
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command <<<
        cat <<EOF > merge_tsvs.py
        import pandas as pd
        import numpy as np
        import sys

        tsvs = pd.read_csv(sys.argv[1], header=None)[0].tolist()
        merged_filename = sys.argv[2]

        dfs = []
        tot = len(tsvs)
        merged_df = pd.DataFrame()
        for i, tsv in enumerate(tsvs):
            if (i+1)%100==0:
                print(f"Loading tsv {i+1}/{tot}...")
            df = pd.read_csv(tsv, sep='\t')
            merged_df = pd.concat([merged_df, df])
        merged_df.to_csv(merged_filename, sep='\t', index=False)
        EOF

        python3 merge_tsvs.py ~{write_lines(tsvs)} ~{merged_filename} > stdout
    >>>

    output {
        File merged_tsv = merged_filename 
    }   
}

task MergeResults {
    input {
        Array[File] tsvs
        String hail_docker
        String merged_filename
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(tsvs, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0
    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command <<<
        head -n 1 ~{tsvs[0]} > ~{merged_filename}; 
        tail -n +2 -q ~{sep=' ' tsvs} >> ~{merged_filename}
    >>>

    output {
        File merged_tsv = merged_filename
    }
}

task GetHailMTSize {
    input {
        String mt_uri
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command <<<
        tot_size=$(gsutil -m du -sh ~{mt_uri} | awk -F '    ' '{ print $1 }')

        cat <<EOF > convert_to_gb.py
        import sys
        size = sys.argv[1]
        unit = sys.argv[2]
        def convert_to_gib(size, unit):
            size_dict = {"KiB": 2**10, "MiB": 2**20, "GiB": 2**30, "TiB": 2**40}
            return float(size) * size_dict[unit] / size_dict["GiB"]
        size_in_gib = convert_to_gib(size, unit)
        print(size_in_gib)        
        EOF

        python3 convert_to_gb.py $tot_size > mt_size.txt
    >>>

    output {
        Float mt_size = read_lines('mt_size.txt')[0]
    }
}

task GetHailMTSizes {
    input {
        Array[String] mt_uris
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command <<<
        tot_size=$(gsutil -m du -shc ~{sep=' ' mt_uris} | awk -F '    ' '{ print $1 }' | tail -n 1)

        cat <<EOF > convert_to_gb.py
        import sys
        size = sys.argv[1]
        unit = sys.argv[2]
        def convert_to_gib(size, unit):
            size_dict = {"KiB": 2**10, "MiB": 2**20, "GiB": 2**30}
            return float(size) * size_dict[unit] / size_dict["GiB"]
        size_in_gib = convert_to_gib(size, unit)
        print(size_in_gib)        
        EOF

        python3 convert_to_gb.py $tot_size > tot_size.txt
    >>>

    output {
        Float mt_size = read_lines('tot_size.txt')[0]
    }
}

task FilterIntervalsToMT {
    input {
        File bed_file
        Float input_size
        String mt_uri
        String bucket_id        
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command <<<
    cat <<EOF > filter_intervals.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    mt_uri = sys.argv[1]
    bed_file = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))
    bucket_id = sys.argv[5]

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.read_matrix_table(mt_uri)
    intervals = hl.import_bed(bed_file, reference_genome='GRCh38')
    mt_filt = mt.filter_rows(hl.is_defined(intervals[mt.locus]))

    filename = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{os.path.basename(mt_uri).split('.mt')[0]}_{os.path.basename(bed_file).split('.bed')[0]}.mt"
    mt_filt.write(filename)
    pd.Series([filename]).to_csv('mt_uri.txt', index=False, header=None)
    EOF
    
    python3 filter_intervals.py ~{mt_uri} ~{bed_file} ~{select_first([runtime_override.cpu_cores, default_attr.cpu_cores])} ~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} ~{bucket_id}
    >>>

    output {
        String mt_filt = read_lines('mt_uri.txt')[0]
    }
}

task FilterIntervalsToVCF {
    input {
        File bed_file
        Float input_size
        String mt_uri
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command <<<
    cat <<EOF > filter_intervals.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    mt_uri = sys.argv[1]
    bed_file = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    if mt_uri.split('.')[-1]=='mt':
        mt = hl.read_matrix_table(mt_uri)
    else:
        mt = hl.import_vcf(mt_uri, reference_genome='GRCh38', force_bgz=True, array_elements_required=False, call_fields=[])
    
    bed_df = pd.read_csv(bed_file, sep='\t')
    try:
        float(bed_df.columns[0]) 
        bed_df = pd.read_csv(bed_file, sep='\t', header=None)
    except:  # has header
        pass
    new_bed_file = f"{os.path.basename(bed_file).split('.bed')[0]}_no_header.bed"
    bed_df.iloc[:,:3].to_csv(new_bed_file, sep='\t', header=None, index=False)
    intervals = hl.import_bed(new_bed_file, reference_genome='GRCh38')
    mt_filt = mt.filter_rows(hl.is_defined(intervals[mt.locus]))

    filename = f"{os.path.basename(mt_uri).split('.mt')[0]}_{os.path.basename(bed_file).split('.bed')[0]}.vcf.bgz"
    hl.export_vcf(mt_filt, filename)
    EOF
    
    python3 filter_intervals.py ~{mt_uri} ~{bed_file} ~{select_first([runtime_override.cpu_cores, default_attr.cpu_cores])} ~{select_first([runtime_override.mem_gb, default_attr.mem_gb])}
    >>>

    output {
        File vcf_filt = "~{basename(mt_uri, '.mt')}_~{basename(bed_file, '.bed')}.vcf.bgz"
    }
}

task SubsetVCFSamplesHail {
    input {
        File vcf_file
        File samples_file  # .txt extension  
        String hail_docker
        String genome_build='GRCh38'
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command <<<
    cat <<EOF > subset_samples.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    vcf_file = sys.argv[1]
    samples_file = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))
    genome_build = sys.argv[5]

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.import_vcf(vcf_file, reference_genome = genome_build, array_elements_required=False, force_bgz=True)
    header = hl.get_vcf_metadata(vcf_file)

    samples = pd.read_csv(samples_file, header=None)[0].astype(str).tolist()
    try:
        # for haploid (e.g. chrY)
        mt = mt.annotate_entries(
            GT = hl.if_else(
                    mt.GT.ploidy == 1, 
                    hl.call(mt.GT[0], mt.GT[0]),
                    mt.GT)
        )
    except:
        pass

    mt_filt = mt.filter_cols(hl.array(samples).contains(mt.s))
    mt_filt = hl.variant_qc(mt_filt)
    mt_filt = mt_filt.filter_rows(mt_filt.variant_qc.AC[1]>0)
    hl.export_vcf(mt_filt, os.path.basename(samples_file).split('.txt')[0]+'.vcf.bgz', metadata=header)
    EOF

    python3 subset_samples.py ~{vcf_file} ~{samples_file} ~{select_first([runtime_override.cpu_cores, default_attr.cpu_cores])} ~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} ~{genome_build}
    >>>

    output {
        File vcf_subset = basename(samples_file, '.txt') + '.vcf.bgz'
    }
}

task MergeMTs {
    input {
        Array[String] mt_uris
        String cohort_prefix
        String bucket_id
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command <<<
    cat <<EOF > merge_mts.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    mt_uris = sys.argv[1].split(',')
    merged_filename = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))
    bucket_id = sys.argv[5]

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    tot_mt = len(mt_uris)
    for i, mt_uri in enumerate(mt_uris):
        if (((i+1)%10)==0):
            print(f"Merging MT {i+1}/{tot_mt}...")
        if i==0:
            mt = hl.read_matrix_table(mt_uri)
        else:
            mt2 = hl.read_matrix_table(mt_uri)
            mt = mt.union_rows(mt2)
    filename = f"{bucket_id}/hail/merged_mt/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{merged_filename}.mt"
    mt.write(filename, overwrite=True)
    pd.Series([filename]).to_csv('mt_uri.txt', index=False, header=None)
    EOF

    python3 merge_mts.py ~{sep=',' mt_uris} ~{cohort_prefix}_merged ~{select_first([runtime_override.cpu_cores, default_attr.cpu_cores])} ~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} ~{bucket_id}
    >>>

    output {
        String merged_mt = read_lines('mt_uri.txt')[0]
    }
}

task MergeHTs {
    input {
        Array[String] ht_uris
        String merged_filename
        String hail_docker
        Float input_size
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command <<<
    cat <<EOF > merge_hts.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    ht_uris = sys.argv[1].split(',')
    merged_filename = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    tot_ht = len(ht_uris)
    merged_df = pd.DataFrame()

    for i, uri in enumerate(ht_uris):
        print(f"{i+1}/{tot_ht}")
        ht = hl.read_table(uri)
        merged_df = pd.concat([merged_df, ht.to_pandas()])

    merged_df.to_csv(f"{merged_filename}.tsv", sep='\t', index=False)
    EOF

    python3 merge_hts.py ~{sep=',' ht_uris} ~{merged_filename} ~{select_first([runtime_override.cpu_cores, default_attr.cpu_cores])} ~{select_first([runtime_override.mem_gb, default_attr.mem_gb])}
    >>>

    output {
        File merged_tsv = merged_filename + '.tsv'
    }
}

task SubsetVCFs {
    input {
        File vcf_uri
        File vcf_idx
        File bed_file
        String output_name
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float relatedness_size = size(vcf_uri, "GB") 
    Float base_disk_gb = 10.0
    RuntimeAttr default_attr = object {
                                      mem_gb: 4,
                                      disk_gb: ceil(base_disk_gb + (relatedness_size) * 5.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command {
        bcftools view -R ~{bed_file} ~{vcf_uri} -o ~{output_name}
        bcftools index -t ~{output_name}
    }

    output {
        File subset_vcf = output_name
        File subset_vcf_idx = output_name + '.tbi'
    }
}

task SplitVcfIntoShards {
  input {
    File input_vcf
    File input_vcf_index
    Int variants_per_shard
    String output_prefix
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -e

    mkdir chunks

    # Extract header
    bcftools view -h ~{input_vcf} > chunks/header.vcf

    # Output body lines (no header), then split
    bcftools view -H ~{input_vcf} | split -l ~{variants_per_shard} - chunks/body_

    # Reconstruct chunked VCFs with header
    for body in chunks/body_*; do
      chunk_name=chunks/~{output_prefix}_$(basename "$body")
      cat chunks/header.vcf "$body" | bgzip -c > "${chunk_name}.vcf.gz"
      tabix -p vcf -f "${chunk_name}.vcf.gz"
    done
 
   >>>

  output {
    Array[File] split_vcfs = glob("chunks/*.vcf.gz")
    Array[File] split_vcf_indexes = glob("chunks/*.vcf.gz.tbi")
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75 + ceil(size(input_vcf,"GiB")*2.5),
    disk_gb: 10 + ceil(size(input_vcf,"GiB")*2.5),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File]? vcfs_idx
    Boolean merge_sort = false
    String outfile_prefix = "concat"
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String outfile_name = outfile_prefix + ".vcf.gz"
  String merge_flag = if merge_sort then "--allow-overlaps" else ""

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(vcfs, "GB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr default_attr = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
    docker: docker_image
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    VCFS_FILE="~{write_lines(vcfs)}"
    if ~{!defined(vcfs_idx)}; then
      cat ${VCFS_FILE} | xargs -n1 tabix
    fi
    bcftools concat -a ~{merge_flag} --output-type z --file-list ${VCFS_FILE} --output "~{outfile_name}"
    tabix -p vcf -f "~{outfile_name}"
  >>>

  output {
    File concat_vcf = outfile_name
    File concat_vcf_idx = outfile_name + ".tbi"
  }
}

task SubsetVcfToContig {
    input {
        File vcf
        File vcf_index
        String contig
        String? args_string
        String prefix
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view \
            ~{vcf} \
            --regions ~{contig} \
            ~{if defined(args_string) then args_string else ""} \
            -Oz -o ~{prefix}.~{contig}.vcf.gz
        tabix -p vcf -f ~{prefix}.~{contig}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.~{contig}.vcf.gz"
        File subset_vcf_index = "~{prefix}.~{contig}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConcatFiles {
    input {
        Array[File] files
        String outfile_name
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat ~{sep=' ' files} > ~{outfile_name}
    >>>

    output {
        File concatenated_file = "~{outfile_name}"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(files, "GB") * 1.2) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SplitQueryVcf {
    input {
        File vcf
        File? vcf_idx
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        svtk vcf2bed -i SVTYPE -i SVLEN ~{vcf} tmp.bed
        cut -f1-4,7-8 tmp.bed > ~{prefix}.bed
        set +o pipefail
        head -1 ~{prefix}.bed > header
        set -o pipefail
        cat header <(awk '{if ($5=="DEL") print}' ~{prefix}.bed )> ~{prefix}.DEL.bed
        cat header <(awk '{if ($5=="DUP") print}' ~{prefix}.bed )> ~{prefix}.DUP.bed
        cat header <(awk '{if ($5=="INS" || $5=="INS:ME" || $5=="INS:ME:ALU" || $5=="INS:ME:LINE1" || $5=="INS:ME:SVA" || $5=="ALU" || $5=="LINE1" || $5=="SVA" || $5=="HERVK" ) print}' ~{prefix}.bed )> ~{prefix}.INS.bed
        cat header <(awk '{if ($5=="INV" || $5=="CPX") print}' ~{prefix}.bed )> ~{prefix}.INV_CPX.bed
        cat header <(awk '{if ($5=="BND" || $5=="CTX") print}' ~{prefix}.bed )> ~{prefix}.BND_CTX.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
        File del_bed = "~{prefix}.DEL.bed"
        File dup_bed = "~{prefix}.DUP.bed"
        File ins_bed = "~{prefix}.INS.bed"
        File inv_bed = "~{prefix}.INV_CPX.bed"
        File bnd_bed = "~{prefix}.BND_CTX.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task BedtoolsClosest {
    input {
        File bed_a
        File bed_b
        String svtype
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eu
        paste <(head -1 ~{bed_a}) <(head -1 ~{bed_b}) | sed -e "s/#//g" > ~{svtype}.bed
        set -o pipefail
        bedtools closest -wo -a <(sort -k1,1 -k2,2n ~{bed_a}) -b <(sort -k1,1 -k2,2n ~{bed_b}) >> ~{svtype}.bed
    >>>

    output {
        File output_bed = "~{svtype}.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SelectMatchedSVs {
    input {
        File input_bed
        String svtype
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(input_bed, ".bed")

    command <<<
        set -euo pipefail
        Rscript /opt/gnomad-lr/scripts/benchmark/R_scripts/R1.bedtools_closest_CNV.R \
            -i ~{input_bed} \
            -o ~{prefix}.comparison
    >>>

    output {
        File output_comp = "~{prefix}.comparison"
    }    

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SelectMatchedINSs {
    input {
        File input_bed
        String svtype
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(input_bed, ".bed")

    command <<<
        Rscript /opt/gnomad-lr/scripts/benchmark/R_scripts/R2.bedtools_closest_INS.R \
            -i ~{input_bed} \
            -o ~{prefix}.comparison
    >>>

    output {
        File output_comp = "~{prefix}.comparison"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConvertToSymbolic {
    input {
        File vcf
        File vcf_idx
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python /opt/gnomad-lr/scripts/helpers/symbalts.py ~{vcf} | \
            python /opt/gnomad-lr/scripts/helpers/abs_svlen.py /dev/stdin | \
            bcftools view -Oz > ~{prefix}.vcf.gz

        tabix -f ~{prefix}.vcf.gz
    >>>

    output {
        File processed_vcf = "~{prefix}.vcf.gz"
        File processed_tbi = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(10 + size(vcf, "GB") * 2),
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task RenameVariantIds {
    input {
        File vcf
        File vcf_index
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        bcftools annotate --set-id '%CHROM-%POS-%REF-%ALT' ~{vcf} -Oz -o ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 4, 
        disk_gb: ceil(size(vcf, "GB")) * 2 + 5,
        boot_disk_gb: 10, 
        preemptible_tries: 2, 
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task FinalizeToFile {
    meta{
        description: "Copies the given file to the specified bucket."
    }

    parameter_meta {
        file: {
            description: "file to finalize",
            localization_optional: true
        }
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
        outdir: "directory to which files should be uploaded"
        name:   "name to set for uploaded file"
    }

    input {
        File file
        String outdir
        String? name

        File? keyfile

        RuntimeAttr? runtime_attr_override
    }



    String gcs_output_dir = sub(outdir, "/+$", "")
    String gcs_output_file = gcs_output_dir + "/" + select_first([name, basename(file)])

    command <<<
        set -euo pipefail

        gsutil -m cp "~{file}" "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FinalizeToDir {

    meta {
        description: "Copies the given file to the specified bucket."
    }

    parameter_meta {
        files: {
            description: "files to finalize",
            localization_optional: true
        }
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
        outdir: "directory to which files should be uploaded"
    }

    input {
        Array[File] files
        String outdir

        File? keyfile

        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_dir = sub(outdir, "/+$", "")

    command <<<
        set -euo pipefail

        cat ~{write_lines(files)} | gsutil -m cp -I "~{gcs_output_dir}"
    >>>

    output {
        String gcs_dir = gcs_output_dir
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FinalizeTarGzContents {
    meta {
        description : "Copies the contents of the given tar.gz file to the specified bucket."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        tar_gz_file : "Gzipped tar file whose contents we'll copy."
        outdir : "Google cloud path to the destination folder."
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
        runtime_attr_override : "[optional] Additional runtime parameters."
    }

    input {
        File tar_gz_file
        String outdir

        File? keyfile

        RuntimeAttr? runtime_attr_override
    }

    # This idiom ensures that we don't accidentally have double-slashes in our GCS paths
    String gcs_output_dir = sub(sub(outdir + "/", "/+", "/"), "gs:/", "gs://")

    command <<<
        set -euo pipefail

        mkdir tmp
        cd tmp
        tar -zxf ~{tar_gz_file}

        gsutil -m cp -r * ~{gcs_output_dir}
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task WriteCompletionFile {

    meta {
        description : "Write a file to the given directory indicating the run has completed."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        outdir : "Google cloud path to the destination folder."
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
    }

    input {
        String outdir
        File? keyfile
    }

    command <<<
        set -euo pipefail

        completion_file="COMPLETED_AT_$(date +%Y%m%dT%H%M%S).txt"
        touch $completion_file

        gsutil cp $completion_file ~{outdir}
    >>>

    #########################

    runtime {
        cpu:                    1
        memory:                 1 + " GiB"
        disks: "local-disk " +  10 + " HDD"
        bootDiskSizeGb:         10
        preemptible:            2
        maxRetries:             2
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
}

task CompressAndFinalize {

    meta {
        description: "Gzip a file and finalize"
    }

    parameter_meta {
        file : "File to compress and finalize."
        outdir : "Google cloud path to the destination folder."
        name : "[optional] Name of the file to write.  If not specified, the name of the input file will be used."
        runtime_attr_override : "[optional] Additional runtime parameters."
    }

    input {
        File file
        String outdir
        String? name

        RuntimeAttr? runtime_attr_override
    }

    String base = basename(file)
    String out = sub(select_first([name, base]), ".gz$", "") +  ".gz"
    # THIS IS ABSOLUTELY CRITICAL: DON'T CHANGE TYPE TO FILE, AS CROMWELL WILL TRY TO LOCALIZE THIS NON-EXISTENT FILE
    String gcs_output_file = sub(outdir, "/+$", "") + "/" + out

    Int disk_size = 2 * ceil(size(file, "GB"))

    command <<<
        set -euo pipefail

        gzip -vkc ~{file} > "~{base}.gz"
        gsutil cp "~{base}.gz" "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FinalizeAndCompress {
    meta {
        description: "Gzip a bunch of files and finalize to the same \'folder\'"
    }

    parameter_meta {
        files : "Files to compress and finalize."
        outdir : "Google cloud path to the destination folder."
        prefix : "[optional] Prefix to add to the output files."
        runtime_attr_override : "[optional] Additional runtime parameters."
    }

    input {
        Array[File] files
        String outdir

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_file = sub(outdir, "/+$", "") + "/" + prefix + "/"

    Int disk_size = 5 * ceil(size(files, "GB"))

    command <<<
        set -euo pipefail

        for ff in ~{sep=' ' files};
        do
            base="$(basename -- ${ff})"
            mv "${ff}" "${base}" && gzip -vk "${base}"
        done

        gsutil -m cp /cromwell_root/*.gz "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             7,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task WriteNamedFile {

    meta {
        description : "Write a file to the given directory with the given name."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        String name
        String outdir
        File? keyfile
    }

    parameter_meta {
        name : "Name of the file to write."
        outdir : "Google cloud path to the destination folder."
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
    }

    command <<<
        set -euo pipefail

        touch "~{name}"

        gsutil cp "~{name}" ~{outdir}
    >>>

    runtime {
        cpu:                    1
        memory:                 1 + " GiB"
        disks: "local-disk " +  10 + " HDD"
        bootDiskSizeGb:         10
        preemptible:            2
        maxRetries:             2
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
}

task ChunkManifest {
    meta {
        description: "Chunk a manifest file into smaller files"
    }

    parameter_meta {
        manifest: "The manifest file to chunk"
        manifest_lines_per_chunk: "The number of lines to include in each chunk"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File manifest
        Int manifest_lines_per_chunk

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        split -a 5 -d --additional-suffix=".txt" -l ~{manifest_lines_per_chunk} ~{manifest} chunk_
    >>>

    output {
        Array[File] manifest_chunks = glob("chunk_*")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

# TODO: investigate if samtools is better for this task
# Sort BAM file by coordinate order
# copied from "dsde_pipelines_tasks/BamProcessing.wdl", with
# customization on the runtime block, and "preemptible_tries" taken out
task SortSam {

    meta {
        description: "Sort a BAM file by coordinate order"
    }

    parameter_meta {
        input_bam: "The BAM file to sort"
        output_bam_basename: "The basename for the output BAM file"
        compression_level: "The compression level for the output BAM file"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File input_bam
        String output_bam_basename
        Int compression_level

        RuntimeAttr? runtime_attr_override
    }

    # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
    # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
    Float sort_sam_disk_multiplier = 3.25
    Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20

    command {
        java -Dsamjdk.compression_level=~{compression_level} -Xms4000m -jar /usr/gitc/picard.jar \
            SortSam \
            INPUT=~{input_bam} \
            OUTPUT=~{output_bam_basename}.bam \
            SORT_ORDER="coordinate" \
            CREATE_INDEX=true \
            CREATE_MD5_FILE=true \
            MAX_RECORDS_IN_RAM=300000 \
            VALIDATION_STRINGENCY=SILENT
    }

    output {
        File output_bam = "~{output_bam_basename}.bam"
        File output_bam_index = "~{output_bam_basename}.bai"
        File output_bam_md5 = "~{output_bam_basename}.bam.md5"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             5,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MakeChrIntervalList {

    meta {
        description: "Make a Picard-style list of intervals for each chromosome in the reference genome"
    }

    parameter_meta {
        ref_dict: "The reference dictionary"
        filter: "A list of strings to filter out of the reference dictionary"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File ref_dict
        Array[String] filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    command <<<
        set -euo pipefail

        grep '^@SQ' ~{ref_dict} | \
            awk '{ print $2 "\t" 1 "\t" $3 }' | \
            sed 's/[SL]N://g' | \
            grep -v -e '^@HD' ~{true='-e' false='' length(filter) > 0} ~{sep=" -e " filter} | \
            tee chrs.txt

        cat chrs.txt | awk '{printf("%s:%d-%d\n", $1,$2,$3)}' > intervalList.intervals

        # Now make another output - a set of individual contig interval list files:
        while read line ; do
            contig=$(echo "${line}" | awk '{print $1}')
            echo "${line}" | awk '{printf("%s:%d-%d\n", $1,$2,$3)}' > contig.${contig}.intervals
        done < chrs.txt
    >>>

    output {
        Array[Array[String]] chrs = read_tsv("chrs.txt")
        File interval_list = "intervalList.intervals"
        Array[String] contig_interval_strings = read_lines("intervalList.intervals")
        Array[File] contig_interval_list_files = glob("contig.*.intervals")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ExtractIntervalNamesFromIntervalOrBamFile {

    meta {
        description: "Pulls the contig names and regions out of an interval list or bed file."
    }

    parameter_meta {
        interval_file: "Interval list or bed file from which to extract contig names and regions."
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File interval_file

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    String interval_tsv_filename = "intervals.tsv"

    command <<<
        set -euo pipefail

        python3 <<CODE

        import re
        interval_re = re.compile('''(.*?):(\d+)-(\d+)''')

        def parse_bed_or_interval_file_for_names(interval_file):

            interval_names = []

            with open(interval_file, 'r') as f:
                if interval_file.endswith(".bed"):
                    for line in f:
                        fields = line.strip().split("\t")
                        contig = fields[0]
                        start = int(fields[1])
                        end = int(fields[2])

                        interval_names.append((contig, start, end))
                else:
                    for line in f:
                        match = interval_re.match(line.strip())
                        if not match:
                            raise RuntimeError(f"BAD NEWS!  DIDN'T MATCH OUR REGEX FOR INTERVAL LISTS!  ARE YOUR FILE EXTENSIONS CORRECT?  Line: {line}")
                        else:
                            contig = match.group(1)
                            start = int(match.group(2))
                            end = int(match.group(3))

                        interval_names.append((contig, start, end))
            return interval_names

        # Get our intervals:
        intervals_names = parse_bed_or_interval_file_for_names("~{interval_file}")

        # Print our interval names to stdout:
        print("Interval info:")
        for interval in intervals_names:
            print(f"{interval[0]}:{interval[1]}-{interval[2]}")
        print()

        # Write out the master interval list and individual interval files:
        with open("~{interval_tsv_filename}", 'w') as f:
            for interval in intervals_names:
                f.write(f"{interval[0]}\t{interval[1]}\t{interval[2]}\n")

        print(f"Wrote all interval info to: ~{interval_tsv_filename}")

        CODE
    >>>

    output {
        Array[Array[String]] interval_info = read_tsv(interval_tsv_filename)
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MakeIntervalListFromSequenceDictionary {

    meta {
        description: "Make a Picard-style list of intervals that covers the given reference genome dictionary, with intervals no larger than the given size limit."
    }

    parameter_meta {
        ref_dict: "The reference dictionary"
        ignore_contigs: "A list of strings to filter out of the reference dictionary"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File ref_dict
        Int max_interval_size = 10000
        Array[String] ignore_contigs = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    String out_master_interval_filename = "intervalList.intervals"
    String interval_tsv_filename = "intervals.tsv"

    command <<<
        set -euo pipefail

        python3 <<CODE

        import re

        def create_intervals(sequence_dictionary_file: str, max_interval_size_bp: int, ignore_contigs: list = []) -> list:

            # First read in our sequence dictionary:
            sq_re = re.compile("@SQ[\t ]+SN:(\w+)[\t ]+LN:(\d+)[\t ].*")
            seq_dict = dict()
            with open(sequence_dictionary_file, 'r') as f:
                for line in f:
                    m = sq_re.match(line)
                    if m and m.group(1) not in ignore_contigs:
                        seq_dict[m.group(1)] = int(m.group(2))

            # Now create a list of intervals with size constraints:
            intervals = list()
            for contig, length in seq_dict.items():
                if length <= max_interval_size_bp:
                    intervals.append((contig, 1, length))
                else:
                    # We need to split the contig into parts:
                    for i in range(1, length+1, max_interval_size_bp):
                        if i+max_interval_size_bp-1 > length:
                            intervals.append((contig, i, length))
                        else:
                            intervals.append((contig, i, i+max_interval_size_bp-1))

                    if intervals[-1][-1] != length:
                        intervals.append((contig, i+max_interval_size_bp, length))

            return intervals

        # Get our intervals:
        intervals = create_intervals("~{ref_dict}", ~{max_interval_size}, ["~{sep='","' ignore_contigs}"])

        # Print our intervals to stdout:
        print(f"Generated {len(intervals)} intervals:")
        for interval in intervals:
            print(f"{interval[0]}:{interval[1]}-{interval[2]}")
        print()

        # Write out the master interval list and individual interval files:
        with open("~{out_master_interval_filename}", 'w') as f:
            with open("~{interval_tsv_filename}", 'w') as f2:
                for interval in intervals:
                    this_interval_filename = f"{interval[0]}.{interval[1]}_{interval[2]}.intervals"

                    f.write(f"{interval[0]}:{interval[1]}-{interval[2]}\n")
                    f2.write(f"{interval[0]}\t{interval[1]}\t{interval[2]}\n")

        print(f"Wrote all intervals to: ~{out_master_interval_filename}")
        print(f"Wrote individual intervals to separate named interval files.")

        CODE
    >>>

    output {
        File interval_list = out_master_interval_filename
        Array[Array[String]] interval_info = read_tsv(interval_tsv_filename)
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CreateIntervalListFileFromIntervalInfo {

    meta {
        description: "Make a Picard-style interval list file from the given interval info."
    }

    parameter_meta {
        contig: "Contig for the interval."
        start: "Start position for the interval."
        end: "End position for the interval."
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        String contig
        String start
        String end
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    String out_interval_list = "~{contig}.~{start}_~{end}.intervals"

    command <<<
        set -euo pipefail

        echo "~{contig}:~{start}-~{end}" > ~{out_interval_list}
    >>>

    output {
        File interval_list = out_interval_list
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "ubuntu:22.04"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CountBamRecords {

    meta {
        description: "Count the number of records in a bam file"
    }

    parameter_meta {
        bam: {
                 localization_optional: true,
                 description: "The bam file"
             }
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 100

    command <<<
        set -eux
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -c ~{bam} > "count.txt" 2>"error.log"
        if [[ -f "error.log" ]]; then
            if [[ -s "error.log" ]]; then echo "samtools has warn/error" && cat "error.log" && exit 1; fi
        fi
    >>>

    output {
        File? samools_error = "error.log"
        Int num_records = read_int("count.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task DownsampleSam {

    meta {
        description : "Downsample the given bam / sam file using Picard/GATK's DownsampleSam tool."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        bam:   "BAM file to be filtered."
        probability : "[Optional] Probability that a read will be emitted (Default: 0.01)."
        strategy : "[Optional] Strategy to use to downsample the given bam file (Default: HighAccuracy)."
        prefix : "[Optional] Prefix string to name the output file (Default: downsampled_reads)."
        extra_args : "[Optional] Extra arguments to pass into DownsampleSam."
    }

    input {
        File bam

        Float probability = 0.01
        String strategy = "HighAccuracy"
        String prefix = "downsampled_reads"

        Int random_seed = 1

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(11 * size(bam, "GiB"))

    command <<<

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        gatk DownsampleSam --VALIDATION_STRINGENCY SILENT --RANDOM_SEED ~{random_seed} -I ~{bam} -O ~{prefix}.bam -S ~{strategy} -P ~{probability} ~{extra_args}
        samtools index -@$np ~{prefix}.bam
    >>>

    output {
        File output_bam = "~{prefix}.bam"
        File output_bam_index = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.2.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Sum {

    meta {
        description : "Sum a list of integers."
    }

    parameter_meta {
        ints: "List of integers to be summed."
        prefix: "[Optional] Prefix string to name the output file (Default: sum)."
    }

    input {
        Array[Int] ints
        String prefix = "sum"

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        python -c "print(~{sep="+" ints})" > ~{prefix}.txt
    >>>

    output {
        Int sum = read_int("~{prefix}.txt")
        File sum_file = "~{prefix}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Uniq {

    meta {
        description : "Find the unique elements in a list of strings."
    }

    parameter_meta {
        strings: "List of strings to be filtered."
    }

    input {
        Array[String] strings

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1

    command <<<
        set -euo pipefail

        sort ~{write_lines(strings)} | uniq > uniq.txt
    >>>

    output {
        Array[String] unique_strings = read_lines("uniq.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Timestamp {

    meta {
        description : "Get the current timestamp."
    }

    parameter_meta {
        dummy_dependencies: "List of dummy dependencies to force recomputation."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        Array[String] dummy_dependencies

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        date --iso-8601=ns > timestamp.txt
    >>>

    output {
        String timestamp = read_string("timestamp.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task BamToBed {

    meta {
        description : "Convert a BAM file to a bed file."
    }

    parameter_meta {
        bam: "BAM file to be converted."
        prefix: "Prefix for the output bed file."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euo pipefail
        bedtools genomecov -ibam ~{bam} -bg > ~{prefix}.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task BamToFastq {

    meta {
        description : "Convert a BAM file to a fastq file."
    }

    parameter_meta {
        bam: "BAM file to be converted."
        prefix: "Prefix for the output fastq file."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3*ceil(size(bam, "GB"))

    command <<<
        set -euo pipefail

        samtools fastq ~{bam} | gzip > ~{prefix}.fq.gz
    >>>

    output {
        File reads_fq = "~{prefix}.fq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MergeFastqs {

    meta {
        description : "Merge fastq files."
    }

    parameter_meta {
        fastqs: "Fastq files to be merged."
        prefix: "Prefix for the output fastq file."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        Array[File] fastqs

        String prefix = "merged"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3 * ceil(size(fastqs, "GB"))

    String disk_type = if disk_size < 375 then "LOCAL" else "HDD"

    Int memory = 8

    command <<<
        FILE="~{fastqs[0]}"
        if [[ "$FILE" =~ \.gz$ ]]; then
            cat ~{sep=' ' fastqs} > ~{prefix}.fq.gz
        else
            cat ~{sep=' ' fastqs} | gzip > ~{prefix}.fq.gz
        fi
    >>>

    output {
        File merged_fastq = "~{prefix}.fq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

# A utility to merge several input BAMs into a single BAM.
task MergeBams {

    meta {
        description : "Merge several input BAMs into a single BAM."
    }

    parameter_meta {
        bams: "Input array of BAMs to be merged."
        prefix: "Prefix for the output BAM."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        Array[File] bams
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(bams, "GB"))

    command <<<
        set -euo pipefail

        samtools merge \
            -p -c --no-PG \
            -@ 2 \
            --write-index \
            -o "~{prefix}.bam##idx##~{prefix}.bam.bai" \
            ~{sep=" " bams}
    >>>

    output {
        File merged_bam = "~{prefix}.bam"
        File merged_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Index {

    meta {
        description : "samtools index a BAM file."
    }

    parameter_meta {
        bam: "BAM file to be indexed."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(bam, "GB"))

    String prefix = basename(bam)

    command <<<
        set -euo pipefail

        mv ~{bam} ~{prefix}
        samtools index ~{basename(prefix)}
    >>>

    output {
        File bai = "~{prefix}.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

# A utility to subset a BAM to specifed loci
task SubsetBam {

    meta {
        description : "Subset a BAM file to a specified locus."
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        File bai
        String locus
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }



    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam
    >>>

    output {
        File subset_bam = "~{prefix}.bam"
        File subset_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ResilientSubsetBam {

    meta {
        description: "For subsetting a high-coverage BAM stored in GCS, without localizing (more resilient to auth. expiration)."
    }

    parameter_meta {
        bam: {
            localization_optional: true
        }
        interval_list_file:  "a Picard-style interval list file to subset reads with"
        interval_id:         "an ID string for representing the intervals in the interval list file"
        prefix: "prefix for output bam and bai file names"
    }

    input {
        File bam
        File bai

        File interval_list_file
        String interval_id
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Array[String] intervals = read_lines(interval_list_file)

    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    String subset_prefix = prefix + "." + interval_id

    command <<<

        # the way this works is the following:
        # 0) relying on the re-auth.sh script to export the credentials
        # 1) perform the remote sam-view subsetting in the background
        # 2) listen to the PID of the background process, while re-auth every 1200 seconds
        source /opt/re-auth.sh
        set -euo pipefail

        # see man page for what '-M' means
        samtools view \
            -bhX \
            -M \
            -@ 1 \
            --verbosity=8 \
            --write-index \
            -o "~{subset_prefix}.bam##idx##~{subset_prefix}.bam.bai" \
            ~{bam} ~{bai} \
            ~{sep=" " intervals} && exit 0 || { echo "samtools seem to have failed"; exit 77; } &
        pid=$!

        set +e
        count=0
        while true; do
            sleep 1200 && date && source /opt/re-auth.sh
            count=$(( count+1 ))
            if [[ ${count} -gt 6 ]]; then exit 0; fi
            if ! pgrep -x -P $pid; then exit 0; fi
        done
    >>>

    output {
        File subset_bam = "~{subset_prefix}.bam"
        File subset_bai = "~{subset_prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Bamtools {

    meta {
        description: "Runs a given bamtools command on a bam file"
    }

    parameter_meta {
        bamfile:    "bam file to run bamtools on"
        cmd:        "bamtools command to run"
        args:       "arguments to pass to bamtools"
    }

    input {
        File bamfile
        String cmd
        String args

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 1 + ceil(2 * size(bamfile, "GiB"))

    command <<<
        bamtools ~{cmd} -in ~{bamfile} -out ~{prefix}.bam ~{args}
    >>>

    output {
        File bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             2,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9.beta"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task DeduplicateBam {

    meta {
        description: "Utility to drop (occationally happening) duplicate records in input BAM"
    }

    parameter_meta {
        aligned_bam: "input BAM file"
        aligned_bai: "input BAM index file"
        same_name_as_input: "if true, output BAM will have the same name as input BAM, otherwise it will have the input basename with .dedup suffix"
        runtime_attr_override: "override default runtime attributes"
    }

    input {
        File aligned_bam
        File aligned_bai

        Boolean same_name_as_input = true

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3 * ceil(size(aligned_bam, "GB"))

    String base = basename(aligned_bam, ".bam")
    String prefix = if (same_name_as_input) then base else (base + ".dedup")

    command <<<
        echo "==========================================================="
        echo "collecting duplicate information"
        time \
            samtools view -@ 1 "~{aligned_bam}" | \
            awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' | \
            sort | uniq -d \
            > "~{aligned_bam}".duplicates.txt
        echo "==========================================================="
        echo "de-duplicating"
        time python3 /opt/remove_duplicate_ont_aln.py \
            --prefix "~{prefix}" \
            --annotations "~{aligned_bam}".duplicates.txt \
            "~{aligned_bam}"
        echo "==========================================================="
        echo "DONE"
        samtools index "~{prefix}.bam"
    >>>

    output {
        File corrected_bam = "~{prefix}.bam"
        File corrected_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.10"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Cat {

    meta {
        description: "Utility to concatenates a group of files into a single output file, with headers in the first line if has_header is true. If has_header is false, the script concatenates the files without headers."
    }

    parameter_meta {
        files:      "text files to combine"
        has_header: "files have a redundant header"
        out:        "[default-valued] output filename"
    }

    input {
        Array[File] files
        Boolean has_header = false
        String out = "out.txt"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(files, "GB"))

    command <<<
        set -euo pipefail

        HAS_HEADER=~{true='1' false='0' has_header}

        if [ HAS_HEADER == 1 ]; then
            ((head -1 ~{files[0]}) && (cat ~{sep=' ' files} | xargs -n 1 tail -n +2)) > ~{out}
        else
            cat ~{sep=' ' files} > ~{out}
        fi
    >>>

    output {
        File combined = out
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ComputeGenomeLength {

    meta {
        description: "Utility to compute the length of a genome from a FASTA file"
    }

    parameter_meta {
        fasta:  "FASTA file"
    }

    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(fasta, "GB"))

    command <<<
        set -euo pipefail

        samtools dict ~{fasta} | \
            grep '^@SQ' | \
            awk '{ print $3 }' | \
            sed 's/LN://' | \
            awk '{ sum += $1 } END { print sum }' > length.txt
    >>>

    output {
        Float length = read_float("length.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ListFilesOfType {

    meta {
        description: "Utility to list files of a given type in a directory"
    }

    parameter_meta {
        gcs_dir:  "input directory"
        suffixes: "suffix(es) for files"
        recurse:  "if true, recurse through subdirectories"
    }

    input {
        String gcs_dir
        Array[String] suffixes
        Boolean recurse = false

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1
    String in_dir = sub(gcs_dir, "/$", "")

    command <<<
        set -ex
        gsutil ls ~{true='-r' false='' recurse} ~{in_dir} > temp.txt
        grep -E '(~{sep="|" suffixes})$' temp.txt > files.txt || touch files.txt
        if [ ! -s files.txt ]; then echo "None found" && exit 1; fi
    >>>

    output {
        Array[String] files = read_lines("files.txt")
        File manifest = "files.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task StopWorkflow {

    meta {
        description: "Utility to stop a workflow"
    }

    parameter_meta {
        reason: "reason for stopping"
    }

    input {
        String reason
    }
    command <<<
        echo -e "Workflow explicitly stopped because \n  ~{reason}." && exit 1
    >>>
    runtime {docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"}
}

task InferSampleName {
    meta {
        description: "Infer sample name encoded on the @RG line of the header section. Fails if multiple values found, or if SM ~= unnamedsample."
    }

    parameter_meta {
        bam: {
            localization_optional: true,
            description: "BAM file"
        }
    }

    input {
        File bam
        File bai
    }



    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} > header.txt
        if ! grep -q '^@RG' header.txt; then echo "No read group line found!" && exit 1; fi

        grep '^@RG' header.txt | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g' | sort | uniq > sample.names.txt
        if [[ $(wc -l sample.names.txt) -gt 1 ]]; then echo "Multiple sample names found!" && exit 1; fi
        if grep -iq "unnamedsample" sample.names.txt; then echo "Sample name found to be unnamedsample!" && exit 1; fi
    >>>

    output {
        String sample_name = read_string("sample.names.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker:         "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task CheckOnSamplenames {

    meta {
        description: "Makes sure the provided sample names are all same, i.e. no mixture of sample names"
    }

    parameter_meta {
        sample_names: "sample names"
    }

    input {
        Array[String] sample_names
    }

    command <<<
        set -eux
        n_sm=$(sort ~{write_lines(sample_names)} | uniq | wc -l | awk '{print $1}')
        if [[ ${n_sm} -gt 1 ]]; then echo "Sample mixture!" && exit 1; fi
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker:         "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

# todo: hook this in to all tasks using LOCAL ssds
task ComputeAllowedLocalSSD {
    # This exists because of the following error message
    #   Task PBFlowcell.ShardLongReads:NA:1 failed. The job was stopped before the command finished. PAPI error code 3.
    #   Execution failed: allocating: creating instance: inserting instance: Number of local SSDs for an instance of type custom-8-15360
    #   should be one of [0, 1, 2, 3, 4, 5, 6, 7, 8, 16, 24], while [9] is requested.
    meta {
        description: "Compute the number of LOCAL ssd's allowed by Google"
    }

    parameter_meta {
        intended_gb: "intended number of GB"
    }

    input {
        Int intended_gb
    }
        Int raw = intended_gb / 375
    command <<<
        if [[ ~{raw} -lt 1 ]]; then  ## we are pushing the boundary here a bit, based on the assumption that input is a convervative estimate
            echo "1" > "result.txt"
        elif [[ ~{raw} -lt 9 ]]; then
            echo ~{raw} > "result.txt"
        elif [[ ~{raw} -lt 16  ]]; then
            echo "16" > "result.txt"
        elif [[ ~{raw} -lt 24  ]]; then
            echo "24" > "result.txt"
        else
            echo "Would request ~{raw} local SSDs, more than possible (24)." && exit 1
        fi
    >>>

    output {
        Int numb_of_local_ssd = read_int("result.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker:         "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task RandomZoneSpewer {

    meta {
        description: "Spews a random GCP zone"
    }

    parameter_meta {
        num_of_zones: "number of zones to spew"
    }

    input {
        Int num_of_zones
    }

    command <<<
        set -eux

        # by no means a perfect solution, but that's not desired anyway
        all_known_zones=("us-central1-a" "us-central1-b" "us-central1-c" "us-central1-f" "us-east1-b" "us-east1-c" "us-east1-d" "us-east4-a" "us-east4-b" "us-east4-c" "us-west1-a" "us-west1-b" "us-west1-c" "us-west2-a" "us-west2-b" "us-west2-c" "us-west3-a" "us-west3-b" "us-west3-c" "us-west4-a" "us-west4-b" "us-west4-c")
        for zone in "${all_known_zones[@]}"; do echo "${zone}" >> zones.txt; done

        shuf zones.txt | head -n ~{num_of_zones} | tr '\n' ' ' > "result.txt"
    >>>

    output {
        String zones = read_string("result.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker:         "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

# Get the current timestamp as a string.
# Levergaes the unix `date` command.
# You can enter your own date format string.
# The default date string is:
#     %Y%m%d_%H%M%S_%N
# which corresponds to a date of the following form:
# For August 10th, 2020 at 16:06:32.7091 EDT (20:06:32.7091 UTC):
#     20200810_200632_709100000
#
task GetCurrentTimestampString {

    meta {
        # The volatile keyword forces call caching to be turned off, which is
        # exactly what we want for this task.
        # For more info see: https://cromwell.readthedocs.io/en/stable/optimizations/VolatileTasks/
        volatile: true
        description: "Get the current timestamp as a string"
    }

    parameter_meta {
        date_format: "The date format string to use. See the unix `date` command for more info."
    }

    input {
        String date_format = "%Y%m%d_%H%M%S_%N"
    }

    String date_file = "the_date_file.txt"

    command {
        date +~{date_format} > ~{date_file}
        cat ~{date_file}
    }

    # ------------------------------------------------
    # Runtime settings:
     runtime {
         docker: "ubuntu:19.10"
         memory: "512 MB"
         disks: "local-disk 10 HDD"
         bootDiskSizeGb: "15"
         preemptible: 0
         cpu: 1
     }

    output {
        String timestamp_string   = read_string(date_file)
    }
}


task GetRawReadGroup {

    meta {
        description: "Get the raw read group from a bam file (assumed to have 1 read group only)"
    }

    parameter_meta {
        gcs_bam_path: "path to bam file in GCS"
        runtime_attr_override: "override the runtime attributes"
    }

    input {
        String gcs_bam_path

        RuntimeAttr? runtime_attr_override
    }

    String out_file = "raw_read_group.txt"

    command {
        set -x

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        # We need to escape the tabs and convert the spaces so that the read group will play nice with downstream processing:
        samtools view -H ~{gcs_bam_path} | grep -m1 '^@RG' | sed -e 's@\t@\\t@g' -e 's@ @_@g' > ~{out_file}

        echo "Raw Read Group:"
        cat ~{out_file}
    }

    output {
        String rg = read_string(out_file)
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.30"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task GetReadsInBedFileRegions {
    meta {
        desciption: "Get the reads from the given bam path which overlap the regions in the given bed file."
    }

    input {
        String gcs_bam_path
        File regions_bed

        String prefix = "reads"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        gcs_bam_path: "GCS URL to bam file from which to extract reads."
        regions_bed: "Bed file containing regions for which to extract reads."
        prefix:    "[default-valued] prefix for output BAM"
        runtime_attr_override: "Runtime attributes override struct."
    }

    Int disk_size = 2 * ceil(size([gcs_bam_path, regions_bed], "GB"))

    command <<<
        set -x
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        samtools view -@${np} -b -h -L ~{regions_bed} ~{gcs_bam_path} | samtools sort - > ~{prefix}.bam
        samtools index -@${np} ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.30"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MapToTsv {

    meta {
        description: "Convert a map to a tsv file"
    }

    parameter_meta {
        my_map: "The map to convert"
        name_of_file: "The name of the file to write to"
    }

    input {
        Map[String, Float] my_map
        String name_of_file
    }

    command <<<
        cp ~{write_map(my_map)} ~{name_of_file}
    >>>

    output {
        File result = "~{name_of_file}"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task CreateIGVSession {
    meta {
        description: "Create an IGV session given a list of IGV compatible file paths.  Adapted / borrowed from https://github.com/broadinstitute/palantir-workflows/blob/mg_benchmark_compare/BenchmarkVCFs ."
    }
    input {
        Array[String] input_bams
        Array[String] input_vcfs
        String reference_short_name
        String output_name

        RuntimeAttr? runtime_attr_override
    }

    Array[String] input_files = flatten([input_bams, input_vcfs])

    command {
        bash /usr/writeIGV.sh ~{reference_short_name} ~{sep=" " input_files} > "~{output_name}.xml"
    }

    output {
        File igv_session = "${output_name}.xml"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "quay.io/mduran/generate-igv-session_2:v1.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SplitContigToIntervals {
    meta {
        author: "Jonn Smith"
        notes: "Splits the given contig into intervals of the given size."
    }

    input {
        File ref_dict
        String contig
        Int size = 200000

        File ref_fasta
        File ref_fasta_fai

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2

    command <<<
        set -euo pipefail

        cat ~{ref_dict} | awk '{print $2,$3}' | grep '^SN' | sed -e 's@SN:@@' -e 's@LN:@@' | tr ' ' '\t' > genome.txt
        grep "~{contig}" genome.txt > genome.contig.txt

        bedtools makewindows -g genome.contig.txt -w ~{size} > ~{contig}.~{size}bp_intervals.bed

        # Make individual bed files from each line:
        while read line ; do
            start=$(echo "${line}" | cut -d $'\t' -f 2)
            end=$(echo "${line}" | cut -d $'\t' -f 3)
            echo "${line}" > ~{contig}.${start}-${end}.single_interval.bed
        done < ~{contig}.~{size}bp_intervals.bed
    >>>

    output {
        File full_bed_file = "~{contig}.~{size}bp_intervals.bed"
        Array[File] individual_bed_files = glob("*.single_interval.bed")
    }

    #########################
        RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ConvertReads {

    meta {
        description : "Convert reads from one format to another."
    }

    parameter_meta {
        reads: "Reads to be converted."
        output_format: "Output format."
    }

    input {
        File reads
        String output_format
    }

    Int disk_size = 3 * ceil(size(reads, "GB"))

    command <<<
        set -euo pipefail

        filename=~{reads}
        input_filetype=${filename*.}
        output_filetype=~{output_format}

        if [[ ($input_filetype == "fastq" || $input_filetype == "fq") && $output_filetype == "fasta" ]]; then
            echo "Converting $input_filetype to $output_filetype"
            seqkit fq2fa $filename -o tmp.out
        elif [ $input_filetype == $output_filetype ]; then
            echo "Input filetype is the output filetype"
            mv $filename tmp.out
        else
            echo "ConvertReads does not know how to convert $input_filetype to $output_filetype"
            exit 1
        fi

        mv tmp.out converted_reads.$output_filetype
    >>>

    output {
        File converted_reads = "converted_reads.~{output_format}"
    }

    runtime {
        cpu: 4
        memory: "8 GiB"
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 0
        docker: "quay.io/broad-long-read-pipelines/lr-pacasus:0.3.0"
    }
}
