import pandas as pd
from typing import Optional, Union, TextIO
from os import PathLike

CLINICAL_COVARIATES_PATH = "/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/Biobank_clinical_data_table_by_blood_sample_deidentified_UNCONSENTED.csv.gz"

def load_covariates_file_time_index(covariate_file: Union[str, PathLike, TextIO] = CLINICAL_COVARIATES_PATH)\
            -> pd.DataFrame:
   """Loads the clinical covariates file indexed by Subject_ID and Event_Date"""
   covariate_table = pd.read_csv(covariate_file,
                                 index_col=["Subject_ID", "Event_Date"],
                                 parse_dates=["Event_Date"])
   return covariate_table


def load_covariates_file_simple_index(covariate_file: Union[str, PathLike, TextIO] = CLINICAL_COVARIATES_PATH, 
                                      id_column: str = "Corrected_Blood_Sample",
                                      index: Optional[pd.Index] = None) -> pd.DataFrame:
    """Loads the clinical covariates file and indexes it according to a specified
    column (by default, "Corrected_Blood_Sample"), without aggregating values (i.e. just takes the
    first value for every index). Optionally also reindexes to match an existing index."""
    covariate_table = pd.read_csv(covariate_file)
    covariate_table = covariate_table.dropna(subset=[id_column])\
                                     .drop_duplicates(id_column)\
                                     .set_index(id_column)
    if index is not None:
        covariate_table = covariate_table.reindex(index)
    return covariate_table


def write_plink_table(table: pd.DataFrame,
                      table_file: Union[str, PathLike, TextIO], 
                      double_id: bool = True,
                      const_fid: Optional[str] = None,
                      id_delim: Optional[str] = None,
                      missing="-9"):
    """Writes a PLINK text format data file with FID+IID columns. Options for processing the ID columns
    match PLINK's options for converting VCFs:
    * double_id prints index twice
    * const_fid uses {const_fid} as FID, index as IID
    * id_delim tries to split index by id_delim, otherwise errors if one of the other two options isn't set.
    """
    if not double_id and const_fid is None and id_delim is None:
        raise ValueError("At least one of double_id, const_fid, or id_delim must be set")
    if double_id and const_fid is not None:
        raise ValueError("Can't use both double_id and const_fid")
    
    if id_delim:
        split_index = table.index.str.split(id_delim)
        fid = split_index.str.get(0)
        iid = split_index.str.get(1)
        if double_id or const_fid is not None:
            fid = fid.where(split_index.str.len() == 2, const_fid if const_fid is not None else table.index)
            iid = iid.where(split_index.str.len() == 2, table.index)
        elif split_index.str.len().ne(2).any():
            raise ValueError(f"Delimiter {id_delim} not found in all IDs, but no backup given")
        new_index = pd.MultiIndex.from_arrays([fid, iid])
    elif const_fid:
        new_index = pd.MultiIndex.from_product([[const_fid], table.index])
    elif double_id:
        new_index = pd.MultiIndex.from_arrays([table.index, table.index])
    new_index.names = ["FID", "IID"]
    reindexed_table = table.set_index(new_index)
    reindexed_table.to_csv(table_file, sep="\t", na_rep=missing)


def load_plink_table(table_file: Union[str, PathLike, TextIO], 
                     double_id: bool = True,
                     const_fid: Optional[str] = None,
                     id_delim: Optional[str] = None, *args, **kwargs):
    """Loads a PLINK text format data file with FID+IID columns. Options for processing the ID columns
    match PLINK's options for converting VCFs:
    * double_id discards FID, but throws AssertionError if FID != IID (unless id_delim is set)
    * const_fid discards FID, but throws AssertionError if FID != const_fid (unless id_delim is set)
    * id_delim uses FID{id_delim}IID
    If none of these options are set, it uses a two-level multiindex with FID and IID.
    """
    if double_id and const_fid is not None:
        raise ValueError("Can't use double_id and const_fid at the same time")
    table = pd.read_csv(table_file, delim_whitespace=True, *args, **kwargs)
    if double_id:
        assert id_delim is not None or table.FID.equals(table.IID)
        table.set_index(table.FID.where(table.FID == table.IID, table.FID.str.cat(table.IID, id_delim)), inplace=True)
    elif const_fid is not None:
        assert id_delim is not None or table.FID.eq(const_fid).all()
        table.set_index(table.IID.where(table.FID == const_fid, table.FID.str.cat(table.IID, id_delim)), inplace=True)
    elif id_delim is not None:
        table.set_index(table.FID.str.cat(table.IID, id_delim), inplace=True)
    else:
        table.set_index(["FID", "IID"], inplace=True)

    if table.index.name is not None:
        # so index was set with one of the interpretation options
        del table.index.name
        table.drop(["FID", "IID"], axis=1, inplace=True)

    return table

def grouped_unconflicting_values(column: pd.Series, 
                                 raise_on_conflict: bool = False, 
                                 groupby_level: Optional[Union[int, str]] = 0, 
                                 *args, **kwargs) -> pd.Series:
    notnull_grouped_column = column.dropna().groupby(level=groupby_level, *args, **kwargs)
    if raise_on_conflict and notnull_grouped_column.nunique().gt(1).any():
        raise ValueError(f"Conflicting values in {column.name}")
    grouped_column = column.groupby(level=groupby_level, *args, **kwargs)
    values = grouped_column.first()
    values = values.mask(notnull_grouped_column.nunique().gt(1))
    return values

def format_boolean_for_plink(in_series: pd.Series) -> pd.Series:
    out_series = in_series.fillna(-1).astype(int) + 1
    return out_series

def grouped_ever_true(column: pd.Series, 
                         groupby_level: Optional[Union[int, str]] = 0, 
                         *args, **kwargs) -> pd.Series:
    grouped_column = column.groupby(level=groupby_level, *args, **kwargs)
    isnull_grouped_column = column.isnull().groupby(level=groupby_level, *args, **kwargs)
    values = grouped_column.any()
    values = values.mask(isnull_grouped_column.all())
    return format_boolean_for_plink(values)
