import inspect
import itertools
import os
import re
import shlex
import shutil
import socket
import subprocess
import sys
import sysconfig
import warnings
from pathlib import Path
from typing import ClassVar, List, Sequence, Optional, Dict

import attr
import hail as hl
import more_itertools
import rpy2.rinterface as ri
import rpy2.rinterface_lib.embedded
import rpy2.robjects as ro
import wrapt
from doit.action import CmdAction
from doit.dependency import MD5Checker
from doit.exceptions import TaskFailed
from doit.task import clean_targets
from rpy2.robjects import r

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
exome_bed_path = Path("padded_Twist_ComprehensiveExome_targets_hg38.bed").resolve()

mt_path = Path(vcf_path.stem).with_suffix(".mt")
qc_path = mt_path.with_suffix(".QC_filtered.mt")
sample_matched_path = qc_path.with_suffix(".sample_matched.mt")
vep_path = sample_matched_path.with_suffix(".VEP.mt")
lof_filtered_path = vep_path.with_suffix(".LOF_filtered.mt")
gwas_filtered_path = sample_matched_path.with_suffix(".GWAS_filtered.mt")
rare_filtered_path = sample_matched_path.with_suffix(".rare_filtered.mt")
exome_filtered_path = sample_matched_path.with_suffix(".exome_filtered.mt")
ld_pruned_path = gwas_filtered_path.with_suffix(".LD_pruned.mt")

base_paths = {
    "full": sample_matched_path,
    "lof": lof_filtered_path,
    "gwas": gwas_filtered_path,
    "rare": rare_filtered_path,
    "exome": exome_filtered_path,
    "ld": ld_pruned_path
}
all_subset_paths = base_paths.copy()
for race in ("white", "black", "hispanic", "asian"):
    for subset in {"full", "ld"}:
        all_subset_paths[f"{race}_{subset}"] = all_subset_paths[subset].with_suffix(f".{race.upper()}_only.mt")
    for subset, suffix in [ ("gwas", "GWAS_filtered"),
                            ("lof", "LOF_filtered"),
                            ("rare", "rare_filtered"),
                            ("exome", "exome_filtered")]:
        all_subset_paths[f"{race}_{subset}"] = all_subset_paths[f"{race}_full"].with_suffix(f".{suffix}.mt")

class TaskDecorator:
    def __new__(cls, *args, **kwargs):
        obj = super().__new__(cls)
        if len(args) == 1 and len(kwargs) == 0 and callable(args[0]):
            obj.__init__()
            return obj(args[0])
        else:
            obj.__init__(*args, **kwargs)
            return obj

    @property
    def consumed_args(self):
        return []

    def generate_tasks(self, wrapped, instance, args, kwargs):
        yield from self.iter_transformed_tasks(wrapped, *args, **kwargs)

    def get_new_signature(self,  wrapped):
        original_signature = inspect.signature(wrapped)
        new_params = []
        for param_name, param in original_signature.parameters.items():
            if param_name not in self.consumed_args:
                new_params.append(param)
        new_signature = original_signature.replace(parameters=new_params)
        return str(new_signature)

    def extract_basename(self, wrapped, task_dict):
        if "basename" in task_dict:
            return task_dict["basename"]
        elif wrapped.__name__.startswith("task_"):
            return wrapped.__name__[5:]
        else:
            raise ValueError("Can't figure out basename")

    def iter_transformed_tasks(self, wrapped, label=None, *args, **kwargs):
        for task_dict in more_itertools.always_iterable(wrapped(*args, **kwargs), dict):
            task_dict["basename"] = self.extract_basename(wrapped, task_dict)
            if label is not None:
                task_dict["name"] = f"{task_dict['name']}_{label}" if "name" in task_dict else label
            yield task_dict

    def __call__(self, wrapped):
        return wrapt.decorator(self.generate_tasks, adapter=wrapt.adapter_factory(self.get_new_signature))(wrapped)


@attr.s(auto_attribs=True)
class each_subset(TaskDecorator):
    subsets: Sequence[str] = tuple(base_paths.keys())
    use_races: bool = True
    use_chroms: bool = True
    include_all_races: bool = True
    include_all_chroms: bool = True

    consumed_args: ClassVar[List[str]] = ['mtfile', 'race', 'subset', 'chrom']

    def generate_tasks(self, wrapped, instance, args, kwargs):
        races: List[Optional[str]] = ['white', 'black', 'hispanic', 'asian']
        chroms: List[Optional] = [str(chrom) for chrom in range(1,23)] + ["X"]
        func_signature = inspect.signature(wrapped)
        for subset in self.subsets:
            for race in races + [None]:
                if race is not None and not self.use_races:
                    continue
                label_without_chrom = f"{race}_{subset}" if race is not None else subset
                path_without_chrom = all_subset_paths[label_without_chrom]
                for chrom in chroms + [None]:
                    if chrom is not None and not self.use_chroms:
                        continue
                    path = path_without_chrom.with_suffix(f".chr{chrom}.mt") if chrom is not None else path_without_chrom
                    if len(self.subsets) == 1:
                        # drop subset label
                        if self.use_races:
                            label = race if race is not None else "all"
                            if chrom is not None:
                                label += f"_chr{chrom}"
                        else:
                            label = chrom if chrom is not None else "all"
                    else:
                        label = f"{label_without_chrom}_chr{chrom}" if chrom is not None else label_without_chrom
                    task_args = {}
                    for arg, value in [('subset', subset), ('race', race), ('chrom', chrom), ('mtfile', path)]:
                        if arg in func_signature.parameters:
                            task_args[arg] = value
                    for task_dict in self.iter_transformed_tasks(wrapped, label, *args, **task_args, **kwargs):
                        basename = task_dict["basename"]
                        name = task_dict["name"]
                        if (chrom is not None or self.include_all_chroms) and (race is not None or self.include_all_races):
                            yield task_dict

                        if chrom is None and self.use_chroms:
                            yield {
                                "basename": basename,
                                "name": f"{name}_chrom_split",
                                "actions": None,
                                "task_dep": [f"{basename}:{name}_chr{chrom}" for chrom in chroms]
                            }
                        if race is None and self.use_races:
                            yield {
                                "basename": basename,
                                "name": f"{name}_race_split",
                                "actions": None,
                                "task_dep": [f"{basename}:{race}_{name}" for race in races]
                            }
                        if race is None and chrom is None and self.use_races and self.use_chroms:
                            yield {
                                "basename": basename,
                                "name": f"{name}_race_and_chrom_split",
                                "actions": None,
                                "task_dep": [f"{basename}:{race}_{name}_chr{chrom}" for race, chrom in itertools.product(races, chroms)]
                            }


class each_race(TaskDecorator):
    def __init__(self, **path_args):
        self.path_args = path_args

    @property
    def consumed_args(self):
        return self.path_args.keys()

    def generate_tasks(self, wrapped, instance, args, kwargs):
        races: List[Optional[str]] = ['white', 'black', 'hispanic', 'asian']
        for race in races + [None]:
            if race is None:
                label = "all"
                resolved_path_args = {key: all_subset_paths[subset] for key, subset in self.path_args.items()}
            else:
                label = race
                resolved_path_args = {key: all_subset_paths[f"{race}_{subset}"] for key, subset in self.path_args.items()}
            for task_dict in self.iter_transformed_tasks(wrapped, label, *args, **resolved_path_args, **kwargs):
                basename = task_dict["basename"]
                name = task_dict["name"]
                yield task_dict

                if race is None:
                    name = name[:-3]  # removes "all"
                    yield  {
                        "basename": basename,
                        "name": f"{name}race_split",
                        "actions": None,
                        "task_dep": [f"{basename}:{name}{race}" for race in races]
                    }


def get_phenotypes_list():
    return build_design_matrix.all_phenotypes


traits_of_interest = ["max_severity_moderate", "severity_ever_severe", "severity_ever_eod", "max_who",
        "severity_ever_increased", "severity_ever_decreased", "who_ever_increased", "who_ever_decreased",
        "recovered", "highest_titer_irnt", "days_onset_to_encounter_log", "covid_encounter_days_log"]
bvl_traits = ["blood_viral_load_bitscore", "blood_viral_load_bitscore_log", "blood_viral_load_bitscore_percentile", "blood_viral_load_detected"]


@attr.s(auto_attribs=True)
class each_phenotype(TaskDecorator):
    use_succeeded_for_all: bool = True

    consumed_args: ClassVar[List[str]] = ["phenotype"]
    phenotype_lists: ClassVar[Dict[str,List[str]]] = {"traits_of_interest": traits_of_interest,
                                                      "blood_viral_load": bvl_traits}

    @classmethod
    def phenotypes_to_run(cls, basename):
        return {
            'task_dep': [f"{basename}:{phenotype}" for phenotype in get_succeeded_phenotypes()]
        }

    def generate_tasks(self, wrapped, instance, args, kwargs):
        for phenotype in get_phenotypes_list():
            yield from self.iter_transformed_tasks(wrapped, phenotype, phenotype=phenotype, *args, **kwargs)
        generated_phenotypes_to_run = {}
        for task_dict in self.iter_transformed_tasks(wrapped, None, phenotype="dummy", *args, **kwargs):  # to get task names
            basename = task_dict["basename"]
            name = task_dict.get("name")
            all_traits_task = {
                "basename": basename,
                'name': f"{name}_all" if name is not None else 'all',
                'actions': None
            }
            if self.use_succeeded_for_all:
                if generated_phenotypes_to_run.get(basename, False):
                    yield {
                        'basename': basename,
                        'name': f"_phenotypes_to_run",
                        'actions': [(self.phenotypes_to_run, [basename])],
                        'task_dep': ['null_model:all']
                    }
                    generated_phenotypes_to_run[basename] = True
                all_traits_task['calc_dep'] = [f"{basename}:_phenotypes_to_run"]
            else:
                all_traits_task['task_dep'] = [f"{basename}:{phenotype}" for phenotype in get_phenotypes_list()]
            yield all_traits_task
            for list_name, phenotypes_list in self.phenotype_lists.items():
                yield {
                    'basename': basename,
                    'name': f"{name}_{list_name}" if name is not None else list_name,
                    "actions": None,
                    "task_dep": [f"{basename}:{phenotype}" for phenotype in phenotypes_list]
                }


class phenotype_pairs(each_phenotype):
    consumed_args: ClassVar[List[str]] = ["trait1", "trait2"]

    @classmethod
    def phenotypes_to_run(cls, basename):
        return {
            'task_dep': [f"{basename}:{trait1}.{trait2}" for trait1, trait2 in itertools.combinations(get_succeeded_phenotypes(), 2)]
        }

    def generate_tasks(self, wrapped, instance, args, kwargs):
        for trait1, trait2 in itertools.permutations(get_phenotypes_list(), 2):
            yield from self.iter_transformed_tasks(wrapped, f"{trait1}.{trait2}", trait1=trait1, trait2=trait2, *args, **kwargs)
        generated_phenotypes_to_run = {}
        for task_dict in self.iter_transformed_tasks(wrapped, None, trait1="dummy", trait2="dummy", *args, **kwargs):  # to get task names
            basename = task_dict["basename"]
            name = task_dict.get("name")
            all_traits_task = {
                "basename": basename,
                'name': f"{name}_all" if name is not None else 'all',
                'actions': None
            }
            if self.use_succeeded_for_all:
                if generated_phenotypes_to_run.get(basename, False):
                    yield {
                        'basename': basename,
                        'name': f"_phenotypes_to_run",
                        'actions': [(self.phenotypes_to_run, [basename])],
                        'task_dep': ['null_model:all']
                    }
                    generated_phenotypes_to_run[basename] = True
                all_traits_task['calc_dep'] = [f"{basename}:_phenotypes_to_run"]
            else:
                all_traits_task['task_dep'] = [f"{basename}:{phenotype}" for phenotype in get_phenotypes_list()]
            yield all_traits_task
            for list_name, phenotypes_list in self.phenotype_lists.items():
                yield {
                    'basename': basename,
                    'name': f"{name}_{list_name}" if name is not None else list_name,
                    "actions": None,
                    "task_dep": [f"{basename}:{trait1}.{trait2}" for trait1, trait2 in itertools.combinations(phenotypes_list, 2)]
                }


class phenotypes_by_tissues(each_phenotype):
    @staticmethod
    def get_tissues_list():
        gtex_models_path = Path("../../resources/metaxcan_data/models/eqtl/mashr").resolve()
        tissues = []
        for model in gtex_models_path.glob("*.db"):
            model_name = model.with_suffix("").name
            tissues.append(model_name[6:])  # "mashr_"
        return tissues

    consumed_args: ClassVar[List[str]] = ["phenotype", "tissue"]

    @classmethod
    def phenotypes_to_run(cls, basename):
        return {
            'task_dep': [f"{basename}:{phenotype}_{tissue}" for phenotype, tissue in itertools.product(get_succeeded_phenotypes(), cls.get_tissues_list())]
        }

    def generate_tasks(self, wrapped, instance, args, kwargs):
        for phenotype, tissue in itertools.product(get_succeeded_phenotypes(), self.get_tissues_list()):
            yield from self.iter_transformed_tasks(wrapped, f"{phenotype}_{tissue}", phenotype=phenotype, tissue=tissue, *args, **kwargs)
        generated_phenotypes_to_run = {}
        for task_dict in self.iter_transformed_tasks(wrapped, None, phenotype="dummy", tissue="dummy", *args, **kwargs):  # to get task names
            basename = task_dict["basename"]
            name = task_dict.get("name")
            all_traits_task = {
                "basename": basename,
                'name': f"{name}_all" if name is not None else 'all',
                'actions': None
            }
            if self.use_succeeded_for_all:
                if generated_phenotypes_to_run.get(basename, False):
                    yield {
                        'basename': basename,
                        'name': f"_phenotypes_to_run",
                        'actions': [(self.phenotypes_to_run, [basename])],
                        'task_dep': ['null_model:all']
                    }
                    generated_phenotypes_to_run[basename] = True
                all_traits_task['calc_dep'] = [f"{basename}:_phenotypes_to_run"]
            else:
                all_traits_task['task_dep'] = [f"{basename}:{phenotype}_{tissue}" for phenotype, tissue in itertools.product(get_phenotypes_list(), self.get_tissues_list())]
            yield all_traits_task
            for list_name, phenotypes_list in self.phenotype_lists.items():
                yield {
                    'basename': basename,
                    'name': f"{name}_{list_name}" if name is not None else list_name,
                    "actions": None,
                    "task_dep": [f"{basename}:{phenotype}_{tissue}" for phenotype, tissue in itertools.product(phenotypes_list, self.get_tissues_list())]
                }


covariates_path = Path(
    "/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/Biobank_clinical_data_table_by_blood_sample_deidentified_UNCONSENTED.csv.gz")
design_matrix_path = Path(
    "/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/jordad05/625_Samples.cohort.QC_filtered.sample_matched.age_flowcell_PCAir_dmatrix.csv")

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
               'default_tasks': ['gwas_plots:traits_of_interest', 'race_prediction'] +
                                [f'build_vcf:{race}_full' for race in ("white", "black", "asian", "hispanic")]
               }


@attr.s(auto_attribs=True)
class bsub(TaskDecorator):
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

    def generate_tasks(self, wrapped, instance, args, kwargs):
        for task_dict in self.iter_transformed_tasks(wrapped, *args, **kwargs):
            if task_dict["actions"] is None or task_dict.get("name", "").startswith("_"):
                yield task_dict
                continue
            basename = task_dict["basename"]
            task_name = f"{basename}:{task_dict['name']}" if 'name' in task_dict else basename
            job_name = task_name.replace(":", "_")
            bsub_action = BsubAction(job_name, self, task_dict["actions"][0])
            bsub_task = {
                "basename": f"bsub_{basename}",
                "actions": [bsub_action],
                "teardown": [bsub_action.kill_me]
            }
            if "name" in task_dict:
                bsub_task["name"] = task_dict["name"]
            yield bsub_task

            bwait_task = task_dict.copy()
            bwait_task.update({
                'actions': [f"bwait -w 'done(%(job_id)s)' || sed -n '2!d;/Done$/!{{q1}}' {job_name}.%(job_id)s.log"]
                           + task_dict["actions"][1:],  # bsub only works on the first action, others are followup/cleanup
                'getargs': {'job_id': (f"bsub_{task_name}", 'job_id')},
                'setup': [f"bsub_{task_name}"]
            })
            yield bwait_task

    def get_bsub_invocation(self, job_name):
        return f"bsub -q {self.queue} " \
                    f"-P {self.project} " \
                    f"-W {self.time} " \
                    f"-n {self.cpus} " \
                    f"-R rusage[mem={self.mem_gb}G] " \
                    f"{'-R himem ' if self.himem else ''}" \
                    f"-J {job_name} " \
                    f"-oo {job_name}.%J.log "

    def format_bsub_command(self, cmd, job_name):
        return self.get_bsub_invocation(job_name) + shlex.quote(cmd)


class BsubAction(CmdAction):
    def __init__(self, job_name: str, bsub_obj: bsub, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.job_name = job_name
        self.bsub_obj = bsub_obj

    def expand_action(self):
        expanded_action = super().expand_action()
        return self.bsub_obj.format_bsub_command(expanded_action, self.job_name)

    def execute(self, out=None, err=None):
        result = super().execute(out, err)
        if result is not None:
            return result
        match = re.match(r"^Job <(\d+)> is submitted to queue <\w+>.$", self.out)
        if match is None:
            action = self.expand_action()  # recreate action string for error message
            return TaskFailed(f"{action}")
        self.values["job_id"] = match.group(1)

    def kill_me(self):
        try:
            job_id = self.values["job_id"]
        except KeyError:
            return
        subprocess.run(f"bkill {job_id}", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)


class bsub_hail(bsub):
    def format_bsub_command(self, cmd, job_name):
        hail_submit_script = f"""ml spark/2.4.5 java
        ml -python
        export HAIL_HOME={sysconfig.get_path("purelib")}/hail
        export SPARK_LOCAL_DIRS=/local/tmp/
        export SPARK_LOG_DIR=/sc/arion/projects/mscic1/scratch/hail/logs/
        export SPARK_WORKER_DIR=/sc/arion/projects/mscic1/scratch/hail/worker/
        export LD_PRELOAD=/usr/lib64/libgslcblas.so
        
        lsf-spark-submit.sh \
        --jars $HAIL_HOME/backend/hail-all-spark.jar \
        --conf spark.driver.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
        --conf spark.executor.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
        --conf spark.executor.extraJavaOptions='-Djava.io.tmpdir=/local/tmp' \
        --conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
        --conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator \
        --executor-memory {self.mem_gb-4}G \
        --driver-memory {self.mem_gb-4}G {cmd}"""

        return self.get_bsub_invocation(job_name) + shlex.quote(hail_submit_script)


def clean_dir_targets(task):
    for target in task.targets:
        target = Path(target)
        if target.exists():
            if target.is_dir():
                shutil.rmtree(target)
            else:
                warnings.warn(f"clean_dir_targets found target {target!s} still existing, but it is not a directory")


@bsub_hail(cpus=128)  # takes about 45 minutes on 128 cores
def task_vcf2mt():
    return {
        "actions": [f"{scriptsdir / 'hail_wgs.py'} convert-vcf-to-mt {vcf_path} {mt_path}"],
        "targets": [mt_path],
        "file_dep": [vcf_path],
        "clean": [clean_dir_targets]
    }


@bsub_hail(cpus=128)  # takes about 40 minutes on 128 cores
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
        "file_dep": [mt_path],
        "clean": [clean_targets, clean_dir_targets]
    }


@bsub_hail(cpus=128)  # takes about 15 minutes on 128 cores
def task_match_samples():
    return {
        "actions": [f"{scriptsdir / 'hail_wgs.py'} match-samples {covariates_path} {qc_path}"],
        "targets": [sample_matched_path],
        "file_dep": [covariates_path, qc_path],
        "clean": [clean_dir_targets]
    }


@bsub_hail(cpus=128)  # takes about 10 minutes on 128 cores (for full)
@each_subset
def task_mt2plink(mtfile):
    return {
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


@bsub_hail(cpus=128) # takes about 10 minutes on 128 cores (for all)
@each_race(input_path='full')
def task_gwas_filter(input_path):
    return {
        "actions": [f"{scriptsdir / 'hail_wgs.py'} gwas-filter {input_path}"],
        "file_dep": [input_path],
        "targets": [input_path.with_suffix(".GWAS_filtered.mt")],
        "clean": [clean_dir_targets]
    }


@bsub_hail(cpus=128) # takes about 10 minutes on 128 cores (for all)
@each_race(input_path="full")
def task_rare_filter(input_path):
    return {
        "actions": [f"{scriptsdir / 'hail_wgs.py'} rare-filter {input_path}"],
        "file_dep": [input_path],
        "targets": [input_path.with_suffix(".rare_filtered.mt")],
        "clean": [clean_dir_targets]
    }


@bsub_hail(cpus=128, mem_gb=12) # takes about 10 minutes on 128 cores (for all)
@each_race(input_path="full")
def task_exome_filter(input_path):
    output_path = input_path.with_suffix(".exome_filtered.mt")
    return {
        "actions": [f"{scriptsdir / 'hail_wgs.py'} restrict-to-bed {input_path} {exome_bed_path} {output_path}"],
        "file_dep": [input_path],
        "targets": [output_path],
        "clean": [clean_dir_targets]
    }


@bsub_hail(cpus=128)
@each_race(input_path="full")
def task_ld_prune(input_path):  # takes about 20 minutes on 128 cores (for all)
    return {
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


@each_subset
def task_plink2snpgds(mtfile):
    output_gds = mtfile.with_suffix(".snp.gds")
    return {
        "actions": [(wrap_r_function("build_snp_gds"), [mtfile.with_suffix("")])],
        "file_dep": [mtfile.with_suffix(".bed"),
                     mtfile.with_suffix(".bim"),
                     mtfile.with_suffix(".fam")],
        "targets": [output_gds],
        "setup": ["initialize_r"],
        "clean": True
    }


@each_race(inpath="ld", outpath="full")
def task_pcair(inpath, outpath):
    return {
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


@bsub_hail(cpus=128)
@each_subset(subsets=["full", "ld"], use_chroms=False, include_all_races=False)
def task_split_races(subset, race):
    if race is None:
        race = "all"
    mtfile = all_subset_paths[subset]
    listfile = f"{race.upper()}.indiv_list.txt"
    outfile = mtfile.with_suffix(f".{race.upper()}_only.mt")
    return {
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


def task_pcrelate():
    return {
        "actions": [(wrap_r_function("run_pcrelate"), [sample_matched_path.with_suffix("")])],
        "targets": [sample_matched_path.with_suffix(".PCRelate.RDS")],
        "file_dep": [ld_pruned_path.with_suffix(".snp.gds"),
                     sample_matched_path.with_suffix(".PCAir.RDS")],
        "setup": ["initialize_r"],
        "clean": True
    }


@bsub_hail(cpus=128)  # takes about 20 minutes on 128 cores (for full)
@each_subset
def task_mt2vcfshards(mtfile):
    output_vcf_dir = mtfile.with_suffix(".shards.vcf.bgz")
    return {
        "actions": [f"{scriptsdir / 'hail_wgs.py'} convert-mt-to-vcf-shards {mtfile} {vcf_path}"],
        "file_dep": [mtfile],
        "targets": [output_vcf_dir],
        "clean": [clean_dir_targets]
    }


@each_subset
def task_build_vcf(mtfile):
    vcf_shards_dir = mtfile.with_suffix(".shards.vcf.bgz")
    output_vcf = mtfile.with_suffix(".vcf.bgz")
    output_tbi = mtfile.with_suffix(".vcf.bgz.tbi")
    return {
        "actions": [f"ml bcftools && bcftools concat --naive -Oz -o {output_vcf!s} {vcf_shards_dir!s}/part-*.bgz",
                    f"ml htslib && tabix {output_vcf!s}"],
        "file_dep": [vcf_shards_dir],
        "targets": [output_vcf, output_tbi],
        "clean": True
    }


@bsub(mem_gb=16, cpus=128)
@each_phenotype(use_succeeded_for_all=False)
def task_null_model(phenotype):
    return {
        "actions": [
            f"ml openmpi && mpirun --mca mpi_warn_on_fork 0 Rscript {scriptsdir / 'mpi_null_model_exhaustive.R'!s} {sample_matched_path.with_suffix('').resolve()!s} {phenotype}"],
        "file_dep": [design_matrix_path, sample_matched_path.with_suffix(".PCRelate.RDS")],
        "targets": [sample_matched_path.with_suffix(f".{phenotype}.null.RDS")],
        "setup": ["initialize_r"],
        "clean": True
    }


def get_succeeded_phenotypes():
    for phenotype in get_phenotypes_list():
        try:
            result = r.readRDS(sample_matched_path.with_suffix(f".{phenotype}.null.RDS"))
        except rpy2.rinterface_lib.embedded.RRuntimeError:
            continue
        if result != ro.NULL and result.rx2['converged'] is True:
            yield phenotype


@bsub(mem_gb=16, cpus=128)
@each_subset
def task_vcf2seqgds_shards(mtfile):
    vcf_shards_dir = mtfile.with_suffix(".shards.vcf.bgz")
    gds_shards_dir = mtfile.with_suffix(".shards.seq.gds")
    return {
        "actions": [
            f"ml openmpi && mpirun --mca mpi_warn_on_fork 0 Rscript {scriptsdir / 'mpi_vcf2gds.R'!s} {vcf_shards_dir!s}"],
        "file_dep": [vcf_shards_dir],
        "targets": [gds_shards_dir],
        "clean": [clean_dir_targets]
    }


@bsub(cpus=128, mem_gb=16)
@each_phenotype
def task_run_gwas(phenotype):
    return {
        "actions": [
            f"ml openmpi && mpirun --mca mpi_warn_on_fork 0 Rscript {scriptsdir / 'mpi_genesis_gwas.R'!s} {sample_matched_path.with_suffix('').resolve()!s} {phenotype}"],
        "targets": [Path(f"{phenotype}.GENESIS.assoc.txt").resolve()],
        "file_dep": [gwas_filtered_path.with_suffix(".shards.seq.gds"),
                     sample_matched_path.with_suffix(f".{phenotype}.null.RDS")],
        "setup": ["initialize_r"],
        "clean": True
    }


@bsub
@each_phenotype
def task_gwas_plots(phenotype):
    assoc_file = Path(f"{phenotype}.GENESIS.assoc.txt").resolve()
    return {
        "actions": [f"Rscript {scriptsdir / 'make_gwas_plot.R'} {phenotype}"],
        "file_dep": [assoc_file],
        "targets": [assoc_file.with_suffix("").with_suffix(".qq.png"),
                    assoc_file.with_suffix("").with_suffix(".manhattan.png")],
        "clean": True
    }


@bsub_hail(cpus=128)
def task_vep():
    return {
        "actions": [f"{scriptsdir / 'hail_wgs.py'} run-vep {sample_matched_path}"],
        "file_dep": [sample_matched_path, "/sc/arion/projects/mscic1/files/WGS/vep/vep_config_script.json"],
        "targets": [vep_path],
        "clean": [clean_dir_targets]
    }


@bsub_hail(cpus=128)
def task_lof_filter():
    return {
        "actions": [f"{scriptsdir / 'hail_wgs.py'} filter-lof-hc {vep_path}"],
        "file_dep": [vep_path],
        "targets": [lof_filtered_path],
        "clean": [clean_dir_targets]
    }


@bsub(cpus=64)
@each_subset
def task_vcf2seqgds_single(mtfile):
    vcf_shards_dir = mtfile.with_suffix(".shards.vcf.bgz")
    output_gds = mtfile.with_suffix(".seq.gds")
    return {
        "actions":  [f"Rscript {scriptsdir / 'seqvcf2gds.R'} {vcf_shards_dir} {output_gds}"],
        "file_dep": [vcf_shards_dir],
        "targets":  [output_gds]
        }


@bsub(cpus=128, mem_gb=16)
@each_subset(subsets=["lof", "rare"], use_chroms=False, use_races=False)
@each_phenotype
def task_run_smmat(phenotype, subset, mtfile):
    return {
        "actions": [f"ml openmpi && mpirun --mca mpi_warn_on_fork 0 Rscript {scriptsdir / 'mpi_genesis_smmat.R'} {mtfile.with_suffix('.seq.gds')} "
                                                f"{sample_matched_path.with_suffix('')}.{phenotype}.null.RDS "
                                                f"{phenotype}.{subset}"],
        "targets": [f"{phenotype}.{subset}.GENESIS.SMMAT.assoc.txt",
                    f"{phenotype}.{subset}.GENESIS.SMMAT.manhattan.png",
                    f"{phenotype}.{subset}.GENESIS.SMMAT.qq.png"],
        "file_dep": [mtfile.with_suffix(".seq.gds"),
                     sample_matched_path.with_suffix(f".{phenotype}.null.RDS")],
        "clean": True
    }


@bsub_hail(cpus=128)
@each_subset(use_chroms=False)
def task_split_chromosomes(mtfile):
    return {
        "actions":  [f"{scriptsdir / 'hail_wgs.py'} split-chromosomes {mtfile}"],
        "file_dep": [mtfile],
        "targets":  [mtfile.with_suffix(f".chr{chr}.mt") for chr in list(range(1,23)) + ['X']],
        "clean":    [clean_dir_targets]
    }


@bsub
@each_subset(subsets=["gwas"], use_chroms=True, include_all_chroms=False)
def task_ld_scores(chrom):
    ldsc_path = gwas_filtered_path.with_suffix(".ld_chr")
    ldsc_path.mkdir(exist_ok=True)
    bfile_prefix = gwas_filtered_path.with_suffix(f".chr{chrom}")
    return {
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


@bsub
@each_phenotype
def task_munge_sumstats(phenotype):
    assoc_file = Path(f"{phenotype}.GENESIS.assoc.txt").resolve()
    return {
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


@bsub(mem_gb=16)
@phenotype_pairs
def task_ld_score_regression(trait1, trait2):
    ldscore_path = gwas_filtered_path.with_suffix(".ld_chr")
    sumstats1 = Path(f"{trait1}.GENESIS.assoc.sumstats.gz").resolve()
    sumstats2 = Path(f"{trait2}.GENESIS.assoc.sumstats.gz").resolve()
    return {
        "actions": ["ml ldsc && " \
                    f"ldsc.py --rg {sumstats1!s},{sumstats2!s} "\
                    f"--ref-ld-chr {ldscore_path!s}/ "\
                    f"--w-ld-chr {ldscore_path!s}/ "\
                    f"--out {trait1}.{trait2}.assoc.rg"],
        "file_dep": [sumstats1, sumstats2, ldscore_path] + [ldscore_path / f"{chr}.l2.ldscore.gz" for chr in list(range(1,23)) + ["X"]],
        "targets":  [f"{trait1}.{trait2}.assoc.rg.log"],
        "clean": True
        }

@bsub(mem_gb=20)
@each_phenotype
def task_metaxcan_harmonize(phenotype):
    gwas_parsing_script = scriptsdir / "summary-gwas-imputation" / "src" / "gwas_parsing.py"
    metadata_file = Path("../../resources/metaxcan_data/reference_panel_1000G/variant_metadata.txt.gz").resolve()
    assoc_file = Path(f"{phenotype}.GENESIS.assoc.txt").resolve()
    harmonized_file = assoc_file.with_suffix(".metaxcan_harmonized.txt")
    return {
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
@phenotypes_by_tissues
def task_s_predixcan(phenotype, tissue):
    s_predixcan_script = scriptsdir / "MetaXcan" / "software" / "SPrediXcan.py"
    model = Path("../../resources/metaxcan_data/models/eqtl/mashr").resolve() / f"mashr_{tissue}.db"
    output_dir = Path("./spredixcan_results").resolve()
    output_dir.mkdir(exist_ok=True)
    assoc_file = Path(f"{phenotype}.GENESIS.assoc.metaxcan_harmonized.txt").resolve()
    output_file = output_dir / f"{phenotype}.mashr_{tissue}.csv"
    return {
        "name": f"{phenotype}_{tissue}",
        "actions": [f"python {s_predixcan_script} "
                    f"--gwas_file {assoc_file} "
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
                    f"--model_db_path {model} "
                    f"--covariance {model.with_suffix('.txt.gz')} "
                    f"--output_file {output_file}"],
        "file_dep": [assoc_file],
        "targets": [output_file],
        "clean": True
    }

@bsub(mem_gb=8)
@each_phenotype
def task_s_multixcan(phenotype):
    s_multixcan_script = scriptsdir / "MetaXcan" / "software" / "SMulTiXcan.py"
    metaxcan_data_dir = Path("../../resources/metaxcan_data").resolve()
    gtex_models_path = metaxcan_data_dir / "models" / "eqtl" / "mashr"
    snp_covariance_file = metaxcan_data_dir / "models" / "gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz"
    predixcan_output_dir = Path("./spredixcan_results").resolve()
    model_names = [model.with_suffix("").name for model in gtex_models_path.glob("*.db")]
    assoc_file = Path(f"{phenotype}.GENESIS.assoc.metaxcan_harmonized.txt").resolve()
    output_file = predixcan_output_dir / f"{phenotype}.smultixcan.txt"
    return {
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


if __name__ == "__main__":
    import doit

    doit.run(globals())
