version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow DownloadAWSFile {
    input {
        String aws_path
        String gcs_folder
        String base_path

        String utils_docker
        RuntimeAttr? runtime_attr_override
    }

    # Strip the base S3 prefix and append to GCS folder to mirror the AWS path structure
    String relative_path = sub(aws_path, "^" + base_path, "")
    String output_gcs_path = sub(gcs_folder, "/+$", "") + "/" + relative_path

    call Helpers.TransferAWSToGCS {
        input:
            aws_path = aws_path,
            output_gcs_path = output_gcs_path,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_override
    }

    output {
        String gcs_path = output_gcs_path
    }
}
