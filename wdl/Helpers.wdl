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

    RuntimeAttr runtime_default = object {
                                      mem_gb: 4,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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

    python3 merge_vcfs.py -i ~{vcf_file} -a ~{annot_vcf_file} -o ~{combined_vcf_name} --cores ~{select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])} --mem ~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} \
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
    
    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

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
    
    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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
    
    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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
    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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
    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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
    
    python3 filter_intervals.py ~{mt_uri} ~{bed_file} ~{select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])} ~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} ~{bucket_id}
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

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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
    
    python3 filter_intervals.py ~{mt_uri} ~{bed_file} ~{select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])} ~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])}
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

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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

    python3 subset_samples.py ~{vcf_file} ~{samples_file} ~{select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])} ~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} ~{genome_build}
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

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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

    python3 merge_mts.py ~{sep=',' mt_uris} ~{cohort_prefix}_merged ~{select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])} ~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} ~{bucket_id}
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

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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

    python3 merge_hts.py ~{sep=',' ht_uris} ~{merged_filename} ~{select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])} ~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])}
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
    RuntimeAttr runtime_default = object {
                                      mem_gb: 4,
                                      disk_gb: ceil(base_disk_gb + (relatedness_size) * 5.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: docker_image
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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
        disk_gb: ceil(size(vcf, "GB")) + 20,
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