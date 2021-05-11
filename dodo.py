import sys, os, re, subprocess, socket, sysconfig, shutil, shlex, itertools, warnings, inspect
import attr, more_itertools
from typing import Union, ClassVar
from pathlib import Path

import rpy2.rinterface_lib.embedded
from doit.dependency import MD5Checker
from doit.exceptions import TaskFailed
from doit.task import clean_targets
from doit.action import CmdAction
import doit
from doit import create_after
import doit.globals
import hail as hl
from rpy2.robjects import r
import rpy2.robjects as ro
import rpy2.rinterface as ri
from decorator import decorate, decorator
import logging

scriptsdir = Path("../../scripts/WGS/").resolve()
sys.path.append(str(scriptsdir))
import build_design_matrix, race_prediction, hail_wgs

num_cpus = len(os.sched_getaffinity(0))


@ro.conversion.py2rpy.register(Path)
def convert_path(path):
    return ro.StrVector([str(path.resolve())])


@ro.conversion.rpy2py.register(ri.BoolSexpVector)
@ro.conversion.rpy2py.register(ri.FloatSexpVector)
@ro.conversion.rpy2py.register(ri.IntSexpVector)
def unpack_return_value(vector):
    if len(vector) == 1:
        value = vector[0]
        if bool(value) == value:
            return bool(value)
        try:
            if int(value) == value:
                return int(value)
        except ValueError:
            pass
        return value
    return vector


def wrap_r_function(funcname):
    pyfuncname = funcname.replace(".", "_")
    return ro.functions.wrap_r_function(r[funcname], pyfuncname)


r.source(scriptsdir / "seqarray_genesis.R")

vcf_path = Path('/sc/arion/projects/mscic1/techdev.incoming/DNA/all_625_samples_cohort_vcf/625_Samples.cohort.vcf.gz')

mt_path = Path(vcf_path.stem).with_suffix(".mt")
qc_path = mt_path.with_suffix(".QC_filtered.mt")
sample_matched_path = qc_path.with_suffix(".sample_matched.mt")
vep_path = sample_matched_path.with_suffix(".VEP.mt")
lof_filtered_path = vep_path.with_suffix(".LOF_filtered.mt")
gwas_filtered_path = sample_matched_path.with_suffix(".GWAS_filtered.mt")
ld_pruned_path = gwas_filtered_path.with_suffix(".LD_pruned.mt")
white_only_path = sample_matched_path.with_suffix(".WHITE_only.mt")
black_only_path = sample_matched_path.with_suffix(".BLACK_only.mt")
asian_only_path = sample_matched_path.with_suffix(".ASIAN_only.mt")
hispanic_only_path = sample_matched_path.with_suffix(".HISPANIC_only.mt")
white_gwas_path = white_only_path.with_suffix(".GWAS_filtered.mt")
black_gwas_path = black_only_path.with_suffix(".GWAS_filtered.mt")
asian_gwas_path = asian_only_path.with_suffix(".GWAS_filtered.mt")
hispanic_gwas_path = hispanic_only_path.with_suffix(".GWAS_filtered.mt")
white_ld_path = ld_pruned_path.with_suffix(".WHITE_only.mt")
black_ld_path = ld_pruned_path.with_suffix(".BLACK_only.mt")
asian_ld_path = ld_pruned_path.with_suffix(".ASIAN_only.mt")
hispanic_ld_path = ld_pruned_path.with_suffix(".HISPANIC_only.mt")

vcf_endpoints = {"full": sample_matched_path,
                 "lof": lof_filtered_path,
                 "gwas": gwas_filtered_path,
                 "ld": ld_pruned_path,
                 "white_full": white_only_path,
                 "black_full": black_only_path,
                 "hispanic_full": hispanic_only_path,
                 "asian_full": asian_only_path,
                 "white_gwas": white_gwas_path,
                 "black_gwas": black_gwas_path,
                 "hispanic_gwas": hispanic_gwas_path,
                 "asian_gwas": asian_gwas_path,
                 "white_ld": white_ld_path,
                 "black_ld": black_ld_path,
                 "hispanic_ld": hispanic_ld_path,
                 "asian_ld": asian_ld_path}

covariates_path = Path(
    "/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/Biobank_clinical_data_table_by_blood_sample_deidentified_UNCONSENTED.csv.gz")
design_matrix_path = Path(
    "/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/jordad05/625_Samples.cohort.QC_filtered.sample_matched.age_flowcell_PCAir_dmatrix.csv")

traits_of_interest = ["max_severity_moderate", "severity_ever_severe", "severity_ever_eod", "max_who",
        "severity_ever_increased", "severity_ever_decreased", "who_ever_increased", "who_ever_decreased", 
        "recovered", "highest_titer_irnt", "days_onset_to_encounter_log", "covid_encounter_days_log"]
bvl_traits = ["blood_viral_load_bitscore", "blood_viral_load_bitscore_log", "blood_viral_load_bitscore_percentile", "blood_viral_load_detected"]

class MD5DirectoryChecker(MD5Checker):
    """Just like the default MD5Checker, but works for directories too.
    For directories it recursively checks every file in the directory.
    If it's a file, defaults to MD5Checker behavior."""

    def check_modified(self, file_path, file_stat, state):
        if os.path.isfile(file_path):
            return super().check_modified(file_path, file_stat, state)
        for dirpath, dirnames, filenames in os.walk(file_path):
            for filename in filenames:
                child_path = os.path.join(dirpath, filename)
                child_file_stat = os.stat(child_path)
                try:
                    child_state = state[child_path]
                except KeyError:
                    return True
                if super().check_modified(child_path, child_file_stat, child_state):
                    return True
        return False

    def get_state(self, dep, current_state):
        if os.path.isfile(dep):
            return super().get_state(dep, current_state)
        new_state = {}
        for dirpath, dirnames, filenames in os.walk(dep):
            for filename in filenames:
                child_path = os.path.join(dirpath, filename)
                current_child_state = current_state if current_state is None else current_state.get(child_path)
                new_state[child_path] = super().get_state(child_path, current_child_state) or current_state[child_path]
        return new_state


DOIT_CONFIG = {'check_file_uptodate': MD5DirectoryChecker,
               'default_tasks': ['gwas_plots', 'race_prediction'] +
                                [f'build_vcf:{race}_full' for race in ("white", "black", "asian", "hispanic")]
               }


@attr.s(auto_attribs=True)
class bsub:
    time: str = "12:00"
    mem_gb: int = 8
    cpus: int = 1
    himem: bool = False
    queue: str = "premium"
    project: str = "acc_mscic"

    bsub_command_template: ClassVar[str] = "bsub -q {queue} " \
                                           "-P {project} " \
                                           "-W {time} " \
                                           "-n {cpus} " \
                                           "-R rusage[mem={mem_gb}G] " \
                                           "{'-R himem ' if himem else ''}" \
                                           "-J {job_name} " \
                                           "-oo {job_name}.%J.log " \
                                           "{'< ' if script else ''} "

    def bsubify_tasks(self, f, *args, **kwargs):
        if f.__name__.startswith("task_"):
            func_basename = f.__name__[5:]
        else:
            func_basename = None  # hopefully doit will deal appropriately with None where a basename is needed
        for task_dict in more_itertools.always_iterable(f(*args, **kwargs), base_type=dict):
            if task_dict["actions"] is None:
                yield task_dict
                continue
            basename = task_dict.get("basename", func_basename)
            if basename is None:
                raise ValueError("Task defined without task_ function name needs a basename")
            task_name = f"{basename}:{task_dict['name']}" if 'name' in task_dict else basename
            job_name = task_name.replace(":", "_")
            bsub_task = {
                "basename": f"bsub_{basename}",
                "actions": [BsubAction(self, task_dict["actions"][0])],
                "teardown": [f"bkill -J {job_name} 0"]
            }
            if "name" in task_dict:
                bsub_task["name"] = task_dict["name"]
            yield bsub_task

            bwait_task = task_dict.copy()
            bwait_task.update({
                'actions': [f"bwait -w 'done(%(job_id)d)' || sed -n '2!d;/Done$/!{{q1}}' {job_name}.%(job_id)d.log"]
                           + task_dict["actions"][1:],  # bsub only works on the first action, others are followup/cleanup
                'getargs': {'job_id': (f"bsub_{task_name}", 'job_id')},
                'setup': [f"bsub_{task_name}"]
            })
            yield bwait_task

    def get_bsub_command(self, cmd):
        return self.bsub_command_template.format_map(attr.asdict(self)) + + shlex.quote(cmd)

    def __call__(self, f):
        return decorate(f, self.bsubify_tasks)


class BsubAction(CmdAction):
    def __init__(self, bsub_obj: bsub, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.bsub_obj = bsub_obj

    def expand_action(self):
        expanded_action = super().expand_action()
        return self.bsub_obj.get_bsub_command(expanded_action)

    def execute(self, out=None, err=None):
        result = super().execute(out, err)
        if result is not None:
            return result
        match = re.match(r"^Job <(\d+)> is submitted to queue <\w+>.$", self.out)
        if match is None:
            action = self.expand_action()  # recreate action string for error message
            return TaskFailed(f"{action}")
        self.values["job_id"] = match.group(1)


class bsub_hail(bsub):
    hail_submit_script: ClassVar[str] = """ml spark/2.4.5
    ml -python
    export HAIL_HOME={sysconfig.get_path("purelib")}/hail
    export SPARK_LOG_DIR=/sc/arion/projects/mscic1/scratch/hail/logs/
    export SPARK_WORKER_DIR=/sc/arion/projects/mscic1/scratch/hail/worker/
    
    lsf-spark-submit.sh \
        --jars $HAIL_HOME/backend/hail-all-spark.jar \
        --conf spark.driver.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
        --conf spark.executor.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
        --conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
        --conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator \
        --executor-memory {mem_gb-4}G \
        --driver-memory {mem_gb-4}G """

    def get_bsub_command(self, cmd):
        return (self.bsub_command_template.format_map(attr.asdict(self)) +
                shlex.quote(self.hail_submit_script.format_map(attr.asdict(self)) + cmd))


def task_initialize_hail():
    os.environ["HAIL_HOME"] = str(Path(sysconfig.get_path("purelib")) / "hail")
    if "LSB_JOBID" in os.environ:
        hostname = socket.gethostname()
        os.environ["SPARK_LOG_DIR"] = str(Path("../../scratch/hail/logs/").resolve())
        os.environ["SPARK_WORKER_DIR"] = str(Path("../../scratch/hail/worker/").resolve())
        os.environ["SPARK16ASTER"] = socket.gethostname()
        os.environ["SPARK_MASTER_PORT"] = "6311"
        os.environ["SPARK_PID_DIR"] = f"/tmp/{os.environ['LSB_JOBID']}_{os.environ['LSB_JOBINDEX']}"
        os.environ["SPARK_WORKER_MEMORY"] = "20G"
        os.environ["SPARK_DAEMON_MEMORY"] = "20G"
        return {"actions": ["lsf-start-spark.sh",
                            (hail_wgs.start_hail, [f"spark://{hostname}:6311"])],
                "teardown": ["lsf-stop-spark.sh", hl.stop]
                }
    else:
        return {
            "actions": [(hail_wgs.start_hail, [f"local[{num_cpus}]"])],
            "teardown": [hl.stop]
        }


def clean_dir_targets(task):
    for target in task.targets:
        target = Path(target)
        if target.exists():
            if target.is_dir():
                shutil.rmtree(target)
            else:
                warnings.warn(f"clean_dir_targets found target {target!s} still existing, but it is not a directory")


@bsub_hail(cpus=256)
def task_vcf2mt():
    return {
        "actions": [f"{scriptsdir / 'hail_wgs.py'} convert-vcf-to-mt {vcf_path} {mt_path}"],
        "targets": [mt_path],
        "file_dep": [vcf_path],
        "clean": [clean_dir_targets]
    }


@bsub_hail(cpus=256)
def task_qc():
    return {
        "actions": [f"{scriptsdir / 'hail_wgs.py'} run-hail-qc {mt_path}"],
        "targets": [qc_path,
                    mt_path.with_suffix(".sample_qc.callrate_hist.png"),
                    mt_path.with_suffix(".sample_qc.gq_hist.png"),
                    mt_path.with_suffix(".sample_qc.dp_hist.png"),
                    mt_path.with_suffix(".sample_qc.dp_v_callrate_scatter.png"),
                    mt_path.with_suffix(".variant_qc.callrate_hist.png"),
                    mt_path.with_suffix(".variant_qc.gq_hist.png"),
                    mt_path.with_suffix(".variant_qc.dp_hist.png"),
                    mt_path.with_suffix(".variant_qc.dp_v_callrate_scatter.png")],
        "file_dep": [mt_path]
        "clean": [clean_targets, clean_dir_targets]
    }


@bsub_hail(cpus=256)
def task_match_samples():
    return {
        "actions": [f"{scriptsdir / 'hail_wgs.py'} match-samples {covariates_path} {qc_path}"],
        "targets": [sample_matched_path],
        "file_dep": [covariates_path, qc_path],
        "clean": [clean_dir_targets]
    }


@bsub_hail(cpus=256)
def task_mt2plink():
    for endpoint_name, endpoint_mtfile in vcf_endpoints.items():
        for name, mtfile in [(f"{endpoint_name}_chr{chrom}", endpoint_mtfile.with_suffix(f".chr{chrom}.mt")) for chrom in list(range(1,23)) + ["X"]] + [(endpoint_name, endpoint_mtfile)]:
            yield {
                "name": name,
                "actions": [f"{scriptsdir / 'hail_wgs.py'} convert-mt-to-plink {mtfile}",
                            f"sed -i s/^chr// {mtfile.with_suffix('.bim')}"],
                "file_dep": [mtfile],
                "targets": [mtfile.with_suffix(".bed"),
                            mtfile.with_suffix(".bim"),
                            mtfile.with_suffix(".fam")],
                "clean": True
            }


def task_king():
    return {
        "actions": ["ml king && king "
                    f'-b {sample_matched_path.with_suffix(".bed")} '
                    f'--kinship --cpus {num_cpus} '
                    f'--prefix {sample_matched_path.with_suffix("")}'],
        "file_dep": [sample_matched_path.with_suffix(".bed"),
                     sample_matched_path.with_suffix(".bim"),
                     sample_matched_path.with_suffix(".fam")],
        "targets": [sample_matched_path.with_suffix(".kin0"),
                    str(sample_matched_path.with_suffix("")) + "X.kin",
                    # pathlib can't handle a suffix that doesn't start with .
                    str(sample_matched_path.with_suffix("")) + "X.kin0"],
        "clean": True
    }


@bsub_hail(cpus=256)
def task_gwas_filter():
    for subset, input_path in [("all", sample_matched_path),
                               ("white", white_only_path),
                               ("black", black_only_path),
                               ("hispanic", hispanic_only_path),
                               ("asian", asian_only_path)]:
        yield {
            "name": subset,
            "actions": [f"{scriptsdir / 'hail_wgs.py'} gwas-filter {input_path}"],
            "file_dep": [input_path],
            "targets": [input_path.with_suffix(".GWAS_filtered.mt")],
            "clean": [clean_dir_targets]
        }


@bsub_hail(cpus=256)
def task_ld_prune():
    for subset, input_path in [("all", gwas_filtered_path),
                               ("white", white_gwas_path),
                               ("black", black_gwas_path),
                               ("hispanic", hispanic_gwas_path),
                               ("asian", asian_gwas_path)]:
        yield {
            "name": subset,
            "actions": [f"{scriptsdir / 'hail_wgs.py'} ld-prune {input_path}"],
            "file_dep": [input_path],
            "targets": [input_path.with_suffix(".LD_pruned.mt")],
            "clean": [clean_dir_targets]
        }


def task_initialize_r():
    return {
        "actions": [wrap_r_function("start_cluster")],
        "teardown": [wrap_r_function("stop_cluster")]
    }


def task_build_snp_gds():
    for name, mtfile in vcf_endpoints.items():
        output_gds = mtfile.with_suffix(".snp.gds")
        yield {
            "name": name,
            "actions": [(wrap_r_function("build_snp_gds"), [mtfile.with_suffix("")])],
            "file_dep": [mtfile.with_suffix(".bed"),
                         mtfile.with_suffix(".bim"),
                         mtfile.with_suffix(".fam")],
            "targets": [output_gds],
            "setup": ["initialize_r"],
            "clean": True
        }


def task_pcair():
    for race, inpath, outpath in [("all", ld_pruned_path, sample_matched_path),
                                  ("white", white_ld_path, white_only_path),
                                  ("black", black_ld_path, black_only_path),
                                  ("hispanic", hispanic_ld_path, hispanic_only_path),
                                  ("asian", asian_ld_path, asian_only_path)]:
        yield {
            "name": race,
            "actions": [(wrap_r_function("run_pcair"), [inpath.with_suffix(".snp.gds"), outpath.with_suffix("")])],
            "targets": [outpath.with_suffix(".PCAir.RDS"),
                        outpath.with_suffix(".PCAir.txt")] +
                       [outpath.with_suffix(f".PC{i}v{j}.pdf")
                        for i, j in itertools.combinations(range(1, 11), 2)],
            "file_dep": [inpath.with_suffix(".snp.gds"),
                         sample_matched_path.with_suffix(".kin0")],
            "setup": ["initialize_r"],
            "clean": True
        }


def task_race_prediction():
    return {"actions": [race_prediction.predict_races],
            "targets": [sample_matched_path.with_suffix(".race_and_PCA.csv"),
                        "WHITE.indiv_list.txt",
                        "BLACK.indiv_list.txt",
                        "ASIAN.indiv_list.txt",
                        "HISPANIC.indiv_list.txt"] +
                       [sample_matched_path.with_suffix(f".PC{i}v{j}.predicted_races.pdf")
                        for i, j in itertools.combinations(range(1, 11), 2)],
            "file_dep": [covariates_path, sample_matched_path.with_suffix(".PCAir.txt")],
            "clean": True
            }


@bsub_hail(cpus=256)
def task_split_races():
    for name, mtfile in [("full", sample_matched_path),
                         ("ld", ld_pruned_path)]:
        for race in ("WHITE", "BLACK", "HISPANIC", "ASIAN"):
            listfile = f"{race}.indiv_list.txt"
            outfile = mtfile.with_suffix(f".{race}_only.mt")
            yield {
                "name": f"{race.lower()}_{name}",
                "actions": [f"{scriptsdir / 'hail_wgs.py'} subset-mt-samples {mtfile} {listfile} {outfile}"],
                "file_dep": [mtfile, listfile],
                "targets": [outfile],
                "clean": [clean_dir_targets]
            }


def task_design_matrix():
    return {
        "actions": [(build_design_matrix.build_design_matrix, [covariates_path, design_matrix_path])],
        "targets": [design_matrix_path],
        "file_dep": [covariates_path, "flowcells.csv", sample_matched_path.with_suffix(".PCAir.txt"),
                     "MSCIC_blood_viral_load_predictions.csv", scriptsdir / "build_design_matrix.py"],
        "clean": ["rm {design_matrix_path!s} *.dist.png"]
    }


def get_phenotypes_list():
    return build_design_matrix.all_phenotypes


def task_pcrelate():
    return {
        "actions": [(wrap_r_function("run_pcrelate"), [sample_matched_path.with_suffix("")])],
        "targets": [sample_matched_path.with_suffix(".PCRelate.RDS")],
        "file_dep": [ld_pruned_path.with_suffix(".snp.gds"),
                     sample_matched_path.with_suffix(".PCAir.RDS")],
        "setup": ["initialize_r"],
        "clean": True
    }


@bsub_hail(cpus=256)
def task_mt2vcfshards():
    for name, mtfile in vcf_endpoints.items():
        output_vcf_dir = mtfile.with_suffix(".shards.vcf.bgz")
        yield {
            "name": name,
            "actions": [f"{scriptsdir / 'hail_wgs.py'} convert-mt-to-vcf-shards {mtfile} {vcf_path}"],
            "file_dep": [mtfile],
            "targets": [output_vcf_dir],
            "clean": [clean_dir_targets]
        }


def task_build_vcf():
    for name, mtfile in vcf_endpoints.items():
        vcf_shards_dir = mtfile.with_suffix(".shards.vcf.bgz")
        output_vcf = mtfile.with_suffix(".vcf.bgz")
        output_tbi = mtfile.with_suffix(".vcf.bgz.tbi")
        yield {
            "name": name,
            "actions": [f"ml bcftools && bcftools concat --naive -Oz -o {output_vcf!s} {vcf_shards_dir!s}/part-*.bgz",
                        f"ml htslib && tabix {output_vcf!s}"],
            "file_dep": [vcf_shards_dir],
            "targets": [output_vcf, output_tbi],
            "clean": True
        }


@bsub(mem_gb=16, cpus=128)
def task_null_model():
    for phenotype in get_phenotypes_list():
        yield {
            "basename": "null_model",
            "name": phenotype,
            "actions": [
                f"ml openmpi && mpirun --mca mpi_warn_on_fork 0 Rscript {scriptsdir / 'mpi_null_model_exhaustive.R'!s} {sample_matched_path.with_suffix('').resolve()!s} {phenotype}"],
            "file_dep": [design_matrix_path, sample_matched_path.with_suffix(".PCRelate.RDS")],
            "targets": [sample_matched_path.with_suffix(f".{phenotype}.null.RDS")],
            "setup": ["initialize_r"],
            "clean": True
        }
    yield {
         # for dependencies and cleaning
         "name": "all",
         "actions": None,
         "task_dep": [f"null_model:{phenotype}" for phenotype in get_phenotypes_list()],
         "clean": [f"rm {sample_matched_path.with_suffix('.*.null.RDS')!s}"]
     }
    yield {
         "name": "traits_of_interest",
         "actions": None,
         "task_dep": [f"null_model:{phenotype}" for phenotype in traits_of_interest]
     }
    yield {
         "name": "blood_viral_load",
         "actions": None,
         "task_dep": ["null_model:blood_viral_load_*"]
     }


def get_succeeded_phenotypes():
    for phenotype in get_phenotypes_list():
        try:
            result = r.readRDS(sample_matched_path.with_suffix(f".{phenotype}.null.RDS"))
        except rpy2.rinterface_lib.embedded.RRuntimeError:
            continue
        if result != ro.NULL and result.rx2['converged'] is True:
            yield phenotype


def gwas_to_run():
    return {
        'task_dep': [f"run_gwas:{phenotype}" for phenotype in get_succeeded_phenotypes()]
    }


def task_gwas_to_run():
    return {
        "actions": [gwas_to_run],
        "task_dep": ["null_model:all"]
    }


@bsub(mem_gb=16, cpus=256)
def task_vcf2gds_shards():
    for name, mtfile in vcf_endpoints.items():
        vcf_shards_dir = mtfile.with_suffix(".shards.vcf.bgz")
        gds_shards_dir = mtfile.with_suffix(".shards.seq.gds")
        yield {
            "name": name,
            "actions": [
                f"ml openmpi && mpirun --mca mpi_warn_on_fork 0 Rscript {scriptsdir / 'mpi_vcf2gds.R'!s} {vcf_shards_dir!s}"],
            "file_dep": [vcf_shards_dir],
            "targets": [gds_shards_dir],
            "clean": [clean_dir_targets]
        }


@bsub(cpus=128, mem_gb=16)
def task_run_gwas():
    for phenotype in get_phenotypes_list():
        yield {
            "name": phenotype,
            "actions": [
                f"ml openmpi && mpirun --mca mpi_warn_on_fork 0 Rscript {scriptsdir / 'mpi_genesis_gwas.R'!s} {sample_matched_path.with_suffix('').resolve()!s} {phenotype}"],
            "targets": [Path(f"{phenotype}.GENESIS.assoc.txt").resolve()],
            "file_dep": [gwas_filtered_path.with_suffix(".shards.seq.gds"),
                         sample_matched_path.with_suffix(f".{phenotype}.null.RDS")],
            "setup": ["initialize_r"],
            "clean": True
        }
    yield {
        "name": "all",
        "actions": None,
        "calc_dep": ["gwas_to_run"],
        "clean": ["rm *.GENESIS.assoc.txt"]
    }
    yield {
         "name": "traits_of_interest",
         "actions": None,
         "task_dep": [f"run_gwas:{phenotype}" for phenotype in traits_of_interest]
     }
    yield {
         "name": "blood_viral_load",
         "actions": None,
         "task_dep": [f"run_gwas:{phenotype}" for phenotype in bvl_traits]
     }


def task_gwas_plots():
    phenotypes_list = get_phenotypes_list()
    yield {
        "name": "all",
        "actions":  [(wrap_r_function("make_gwas_plots"), [phenotypes_list])],
        "task_dep": ["run_gwas:all"],
        "calc_dep": ["gwas_to_run"],
        "setup":    ["initialize_r"],
        "clean":    ["rm *.GENESIS.qq.png *.GENESIS.manhattan.png"]
    }
    yield {
        "name": "blood_viral_load",
        "actions": [(wrap_r_function("make_gwas_plots"), [bvl_traits])],
        "task_dep": ["run_gwas:blood_viral_load"],
        "setup": ["initialize_r"],
        "clean": ["rm *.GENESIS.qq.png *.GENESIS.manhattan.png"]
    }


@bsub_hail(cpus=256)
def task_vep():
    return {
        "actions": [f"{scriptsdir / 'hail_wgs.py'} run-vep {sample_matched_path}"],
        "file_dep": [sample_matched_path, "/sc/arion/projects/mscic1/files/WGS/vep/vep_config_script.json"],
        "targets": [vep_path]
        "clean": [clean_dir_targets]
    }


@bsub_hail(cpus=256)
def task_lof_filter():
    return {
        "actions": [f"{scriptsdir / 'hail_wgs.py'} filter-lof-hc {vep_path}"],
        "file_dep": [vep_path],
        "targets": [lof_filtered_path],
        "clean": [clean_dir_targets]
    }


def task_build_seq_gds():
    for name, mtfile in vcf_endpoints.items():
        gds_shards_dir = mtfile.with_suffix(".shards.seq.gds")
        output_gds = mtfile.with_suffix(".seq.gds")
        yield {
            "name":     name,
            "actions":  [wrap_r_function("build_seq_gds")],
            "file_dep": [gds_shards_dir],
            "targets":  [output_gds],
            "setup":    ["initialize_r"]
            }


def task_run_smmat():
    for phenotype in get_phenotypes_list():
        yield {
            "name": phenotype,
            "actions": [(wrap_r_function("run_smmat"), [lof_filtered_path.with_suffix(".seq.gds"),
                                                        sample_matched_path.with_suffix(f".{phenotype}.null.RDS"),
                                                        phenotype])],
            "targets": [f"{phenotype}.GENESIS.SMMAT.assoc.txt",
                        f"{phenotype}.GENESIS.SMMAT.manhattan.png",
                        f"{phenotype}.GENESIS.SMMAT.qq.png"],
            "file_dep": [lof_filtered_path.with_suffix(".seq.gds"),
                         sample_matched_path.with_suffix(f".{phenotype}.null.RDS")],
            "setup": ["initialize_r"],
            "clean": True
        }
    yield {
         "name": "traits_of_interest",
         "actions": None,
         "task_dep": [f"run_smmat:{phenotype}" for phenotype in traits_of_interest]
     }

    
def task_split_chromosomes():
    for name, mtfile in vcf_endpoints.items():
        yield {
            "name": name,
            "actions":  [(hail_wgs.split_chromosomes, [mtfile])],
            "file_dep": [mtfile],
            "targets":  [mtfile.with_suffix(f".chr{chr}.mt") for chr in list(range(1,23)) + ['X']],
            "clean":    [clean_dir_targets],
            "setup":    ["initialize_hail"]
        }


@bsub
def task_ld_scores():
    ldsc_path = gwas_filtered_path.with_suffix(".ld_chr")
    ldsc_path.mkdir(exist_ok=True)
    for chrom in list(range(1,23)) + ["X"]:
        bfile_prefix = gwas_filtered_path.with_suffix(f".chr{chrom}")
        yield {
            "name":    f"chr{chrom}",
            "actions": ["ml ldsc && " \
                        f"ldsc.py --bfile {bfile_prefix} " \
                        f"--l2 --ld-wind-kb 1000 --out {ldsc_path!s}/{chrom}"],
            "targets": [ldsc_path / f"{chrom}.l2.M",
                        ldsc_path / f"{chrom}.l2.M_5_50",
                        ldsc_path / f"{chrom}.l2.ldscore.gz"],
            "file_dep": [gwas_filtered_path.with_suffix(f".chr{chrom}.bed"),
                         gwas_filtered_path.with_suffix(f".chr{chrom}.bim"),
                         gwas_filtered_path.with_suffix(f".chr{chrom}.fam")],
            "clean": True
            }
    yield {
            "name":     "all",
            "actions":  None,
            "task_dep": [f"ld_scores:chr{chrom}" for chrom in list(range(1,23)) + ["X"]]
            }


@bsub
def task_munge_sumstats():
    for phenotype in get_phenotypes_list():
        assoc_file = Path(f"{phenotype}.GENESIS.assoc.txt").resolve()
        yield {
                "name": phenotype,
                "file_dep": [assoc_file],
                "targets": [assoc_file.with_suffix(".sumstats.gz")],
                "actions": ["ml ldsc && " \
                            f"munge_sumstats.py --sumstats {assoc_file!s} " \
                            f"--out {assoc_file.with_suffix('')!s} " \
                            "--snp variant.id --N-col n.obs --frq freq " \
                            "--a1 effect.allele --a2 other.allele " \
                            "--p Score.pval --signed-sumstats Est,0"],
                "clean": True
                }
    yield {
         "name": "traits_of_interest",
         "actions": None,
         "task_dep": [f"munge_sumstats:{phenotype}" for phenotype in traits_of_interest]
     }


@bsub(mem_gb=16)
def task_ld_score_regression():
    ldscore_path = gwas_filtered_path.with_suffix(".ld_chr")
    for trait1, trait2 in itertools.permutations(get_phenotypes_list(), 2):
        sumstats1 = Path(f"{trait1}.GENESIS.assoc.sumstats.gz").resolve()
        sumstats2 = Path(f"{trait2}.GENESIS.assoc.sumstats.gz").resolve()
        yield {
            "name": f"{trait1}.{trait2}",
            "actions": ["ml ldsc && " \
                        f"ldsc.py --rg {sumstats1!s},{sumstats2!s} "\
                        f"--ref-ld-chr {ldscore_path!s}/ "\
                        f"--w-ld-chr {ldscore_path!s}/ "\
                        f"--out {trait1}.{trait2}.assoc.rg"],
            "file_dep": [sumstats1, sumstats2, ldscore_path] + [ldscore_path / f"{chr}.l2.ldscore.gz" for chr in list(range(1,23)) + ["X"]],
            "targets":  [f"{trait1}.{trait2}.assoc.rg.log"],
            "clean": True
            }
    yield {
            "name": "traits_of_interest",
            "actions": None,
            "task_dep": [f"ld_score_regression:{trait1}.{trait2}" for trait1, trait2 in itertools.combinations(traits_of_interest, 2)]
            }

@bsub(mem_gb=20)
def task_metaxcan_harmonize():
    gwas_parsing_script = scriptsdir / "summary-gwas-imputation" / "src" / "gwas_parsing.py"
    metadata_file = Path("../../resources/metaxcan_data/reference_panel_1000G/variant_metadata.txt.gz").resolve()
    for phenotype in get_phenotypes_list():
        assoc_file = Path(f"{phenotype}.GENESIS.assoc.txt").resolve()
        harmonized_file = assoc_file.with_suffix(".metaxcan_harmonized.txt")
        yield {
            "name": phenotype,
            "actions": [f"python {gwas_parsing_script} "
                        "-separator ' ' "
                        f"-gwas_file {assoc_file!s} "
                        f"-snp_reference_metadata {metadata_file!s} METADATA "
                        "-output_column_map variant.id variant_id "
                        "-output_column_map other.allele non_effect_allele "
                        "-output_column_map effect.allele effect_allele "
                        "-output_column_map Est effect_size "
                        "-output_column_map Est.SE standard_error "
                        "-output_column_map Score.pval pvalue "
                        "-output_column_map chr chromosome --chromosome_format "
                        "-output_column_map pos position "
                        "-output_column_map n.obs sample_size "
                        "-output_column_map freq frequency "
                        "-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size "
                        f"-output {harmonized_file}"
                        ],
            "targets": [harmonized_file],
            "file_dep": [assoc_file],
            "clean": True
        }


@bsub
def task_s_predixcan():
    s_predixcan_script = scriptsdir / "MetaXcan" / "software" / "SPrediXcan.py"
    gtex_models_path = Path("../../resources/metaxcan_data/models/eqtl/mashr").resolve()
    output_dir = Path("./spredixcan_results").resolve()
    output_dir.mkdir(exist_ok=True)
    model_names = []
    for model in gtex_models_path.glob("*.db"):
        model_name = model.with_suffix("").name
        model_names.append(model_name)
        for phenotype in get_phenotypes_list():
            assoc_file = Path(f"{phenotype}.GENESIS.assoc.metaxcan_harmonized.txt").resolve()
            output_file = output_dir / f"{phenotype}.{model_name}.csv"
            yield {
                "name": f"{phenotype}_{model_name}",
                "actions": [f"python {s_predixcan_script!s} "
                            f"--gwas_file {assoc_file!s} "
                            "--snp_column panel_variant_id "
                            "--chromosome_column chromosome "
                            "--position_column position "
                            "--effect_allele_column effect_allele "
                            "--non_effect_allele_column non_effect_allele "
                            "--beta_column effect_size "
                            "--se_column standard_error "
                            "--pvalue_column pvalue "
                            "--model_db_snp_key varID "
                            "--keep_non_rsid "
                            "--additional_output "
                            "--overwrite "
                            "--throw "
                            f"--model_db_path {model!s} "
                            f"--covariance {model.with_suffix('.txt.gz')!s} "
                            f"--output_file {output_file!s}"],
                "file_dep": [assoc_file],
                "targets": [output_file],
                "clean": True
            }
    yield {
        "name": "traits_of_interest",
        "actions": None,
        "task_dep": [f"s_predixcan:{phenotype}_{model_name}" for phenotype, model_name in itertools.product(traits_of_interest, model_names)]
    }

@bsub(mem_gb=8)
def task_s_multixcan():
    s_multixcan_script = scriptsdir / "MetaXcan" / "software" / "SMulTiXcan.py"
    metaxcan_data_dir = Path("../../resources/metaxcan_data").resolve()
    gtex_models_path = metaxcan_data_dir / "models" / "eqtl" / "mashr"
    snp_covariance_file = metaxcan_data_dir / "models" / "gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz"
    predixcan_output_dir = Path("./spredixcan_results").resolve()
    model_names = [model.with_suffix("").name for model in gtex_models_path.glob("*.db")]
    for phenotype in get_phenotypes_list():
        assoc_file = Path(f"{phenotype}.GENESIS.assoc.metaxcan_harmonized.txt").resolve()
        output_file = predixcan_output_dir / f"{phenotype}.smultixcan.txt"
        yield {
            "name": phenotype,
            "actions": [f"python {s_multixcan_script!s} "
                        f"--models_folder {gtex_models_path!s} "
                        '--models_name_pattern "mashr_(.*).db" '
                        f"--snp_covariance {snp_covariance_file!s} "
                        f"--metaxcan_folder {predixcan_output_dir!s} "
                        f'--metaxcan_filter "{phenotype}.mashr_(.*).csv" '
                        '--metaxcan_file_name_parse_pattern "(.*).mashr_(.*).csv" '
                        f"--gwas_file {assoc_file!s} "
                        "--snp_column panel_variant_id "
                        "--chromosome_column chromosome "
                        "--position_column position "
                        "--effect_allele_column effect_allele "
                        "--non_effect_allele_column non_effect_allele "
                        "--beta_column effect_size "
                        "--se_column standard_error "
                        "--pvalue_column pvalue "
                        "--model_db_snp_key varID "
                        "--keep_non_rsid "
                        "--cutoff_condition_number 30 "
                        "--throw "
                        f"--output {output_file!s}"],
            "file_dep": [assoc_file] + [predixcan_output_dir / f"{phenotype}.{model_name}.csv" for model_name in model_names],
            "targets": [output_file],
            "clean": True
        }
    yield {
        "name": "traits_of_interest",
        "actions": None,
        "task_dep": [f"s_multixcan:{phenotype}" for phenotype in traits_of_interest]
    }

if __name__ == "__main__":
    import doit

    doit.run(globals())
