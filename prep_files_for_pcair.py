import socket, os
import hail as hl
from pathlib import Path
tmp_path = Path("/sc/arion/projects/mscic1/scratch/hail/tmp/").resolve()
hostname = socket.gethostname()
port = os.environ["SPARK_MASTER_PORT"]
hl.init(master=f"spark://{hostname}:{port}", tmp_dir=str(tmp_path), min_block_size=128)
sc = hl.spark_context()



hl.utils.java.info("hi I'm an end user and I'm writing to this log too")

mt_path = Path('/sc/arion/projects/mscic1/files/WGS/625_Samples.cohort.QC_filtered.sample_matched.mt')
hl.utils.java.info(f"loading matrix table from {mt_path!s}")
mt = hl.read_matrix_table(str(mt_path))
hl.utils.java.info(f"dropping duplicate sample PICR7147T1a and dropping timepoints from sample names")
mt = mt.filter_cols(mt.s != "PICR7147T1a")
mt = mt.annotate_cols(Subject_ID=mt.s.split("T")[0]).key_cols_by("Subject_ID").drop("s")
hl.utils.java.info("imputing sex for PLINK export")
imputed_sex = hl.impute_sex(mt.GT, male_threshold=0.6, female_threshold=0.5)
mt = mt.annotate_cols(is_female=imputed_sex[mt.Subject_ID].is_female)
hl.utils.java.info("filtering to MAF >= 0.01")
mt = mt.filter_rows(mt.gt_stats.AF.all(lambda af: af >= 0.01))
hl.utils.java.info("pruning LD")
pruned_variant_table = hl.ld_prune(mt.GT, r2=0.2, bp_window_size=500000)
mt = mt.filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))

plink_path = mt_path.with_suffix(".maf_and_ld")
hl.utils.java.info("exporting to PLINK format files {plink_path!s}.bed|bim|fam")
hl.export_plink(mt, str(plink_path), fam_id=mt.Subject_ID, ind_id=mt.Subject_ID, is_female=mt.is_female)
hl.utils.java.info("done!")
hl.stop()
