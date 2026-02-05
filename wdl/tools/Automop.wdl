version 1.0

import "Structs.wdl"

workflow Automop {
    input {
        String workspace_namespace
        String workspace_name
        String user
        Boolean dry_run
        String? root_directory

        RuntimeAttr? runtime_attr_mop
    }

    call MopTask {
        input:
            workspace_namespace = workspace_namespace,
            workspace_name = workspace_name,
            user = user,
            dry_run = dry_run,
            root_directory = root_directory,
            runtime_attr_override = runtime_attr_mop
    }

    output {
        File fissfc_log = MopTask.fissfc_log
    }
}

task MopTask {
    input {
        String workspace_namespace
        String workspace_name
        String user
        Boolean dry_run
        String? root_directory
        RuntimeAttr? runtime_attr_override
    }
    
    command <<<
        set -euo pipefail
        
        cat <<'EOF' > script.py
import subprocess
from datetime import datetime
from google.cloud import bigquery
import pytz

def main(workspace_namespace, workspace_name, user, root_directory):
    units = ['bytes', 'KiB', 'MiB', 'GiB', 'TiB', 'PiB']

    cmd = ['fissfc', '--yes', '--verbose', 'mop', '-w', workspace_name, '-p', workspace_namespace]
    
    if root_directory and root_directory != '' and root_directory != '/':
        cmd.extend(['-d', root_directory])
    
    ~{if dry_run then "cmd.append('--dry-run')" else ""}
    
    mop_process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    
    size_found = False
    run_successful = False
    with open('fissfc_log.log', 'w') as fissfc_log:
        for line in iter(mop_process.stdout.readline, b''):
            line = line.decode()
            fissfc_log.write(line)
            if line.startswith('Total Size: '):
                size_with_unit = line.rstrip()[12:]
                size, unit = size_with_unit.split(' ')
                size_in_bytes = int(float(size) * 1024 ** units.index(unit))
                size_found = True
            if line.startswith('No files to mop in'):
                size_in_bytes = 0
                print('No files to mop.')
                size_found = True
                run_successful = True
            if line.startswith('Operation completed over'):
                print('Mopping complete!')
                run_successful = True
    
    if not size_found:
        raise RuntimeError('No total deleted size found in fissfc output.')
    
    if not run_successful:
        raise RuntimeError('Did not receive "Operation completed" message from fissfc output.')
    
    mop_event = {
            'user': user,
            'datetime': datetime.now(pytz.timezone('US/Eastern')).strftime('%Y-%m-%d-%H-%M-%S'),
            'workspace_namespace': workspace_namespace,
            'workspace_name': workspace_name,
            'size_deleted': size_in_bytes
        }
    db = bigquery.Client(project='broad-dsde-methods-automop')
    db.insert_rows(db.get_table('broad-dsde-methods-automop.automop.mop_events'), [mop_event])


if __name__ == '__main__':
    main('~{workspace_namespace}', '~{workspace_name}', '~{user}', '~{default="" root_directory}')
EOF
        python script.py
    >>>

    output {
        File fissfc_log = "fissfc_log.log"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 16,
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 0,
        max_retries: 0,
        docker: "us.gcr.io/broad-dsde-methods/automop:0.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: select_first([runtime_attr.docker, default_attr.docker])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}