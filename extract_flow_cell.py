from pathlib import Path
import pysam
from tqdm import tqdm

bams_file_list = Path("/sc/arion/projects/mscic1/techdev.incoming/DNA/all_625_samples_bams_v1.FREEZE/all_625_samples_bams_v1.FREEZE.fofn")

with bams_file_list.open() as list_file:
    for bam_path_str in tqdm(list_file, total=625):
        bam_path = Path(bam_path_str.strip().replace("mscic", "mscic1"))
        sample_name = bam_path.with_suffix("").name
        flowcells = set()
        with pysam.AlignmentFile(bam_path) as bam_file:
            for read in bam_file.head(1000):
                instrument, run, flowcell, lane, tile, x, y = read.query_name.split(":")
                flowcells.add(flowcell)
        print(sample_name + "," + "+".join(sorted(flowcells)))

