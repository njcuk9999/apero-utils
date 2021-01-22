"""
Log class to hold log information for APERO products

@author: vandal
"""
from typing import List, Optional, Union, Tuple

import numpy as np
import pandas as pd

class Log:
    """Log information for DRS products using pandas dataframe"""

    def __init__(self,
                 full_df: pd.DataFrame,
                 recipe: Optional[Union[str, List]] = None):
        """__init__.

        The input is a dataframe to avoid parsing file for every object
        created.

        :param full_df: Dataframe with all log entries to consider
        :type full_df: pd.DataFrame
        :param recipe: Recipe(s) to filter out log entries
        :type recipe: Optional[Union[str, List]]
        """
        # Make list of recipe
        if isinstance(recipe, str):
            recipe = [recipe]
        elif isinstance(recipe, list):
            pass
        elif recipe is not None:
            raise TypeError('recipe must be a string or a list')

        if recipe is not None:
            self._df = full_df.loc[full_df.RECIPE.isin(recipe)].copy()
        else:
            self._df = full_df.copy()

    @property
    def df(self) -> pd.DataFrame:
        """df.

        :rtype: pd.DataFrame
        """
        return self._df

    @property
    def tot_num(self) -> int:
        """tot_num.

        :rtype: int
        """
        return len(self.df)

    @property
    def num_night(self) -> pd.Series:
        """num_night.

        Count entries per night (directory)

        :rtype: pd.Series
        """
        return self.df.groupby('DIRECTORY').size()

    @property
    def master_mask(self) -> pd.Series:
        """master_mask.

        :rtype: pd.Series
        """
        return self.df['RUNSTRING'].str.contains('--master')

    @property
    def master_nights(self) -> np.ndarray:
        """master_nights.

        :rtype: np.ndarray
        """
        return self.master_mask[self.master_mask].index.values

    @property
    def master_num(self) -> int:
        """master_recipe_num_logfits.

        :rtype: int
        """
        return self.master_nights.size

    @property
    def num(self) -> int:
        """recipe_num_logfits.

        :rtype: int
        """
        return self.log_tot_num - self.master_recipe_num_logfits

    @property
    def qc_failed(self) -> pd.DataFrame:
        """qc_failed.

        :rtype: pd.DataFrame
        """
        return self.df[~self.df.PASSED_ALL_QC]

    @property
    def ended_false(self) -> pd.DataFrame:
        """ended_false.

        :rtype: pd.DataFrame
        """
        return self.df[~self.df.ENDED]

    def get_align_count(self, other) -> Tuple[pd.Series]:
        """get_align_count.

        Return count per night, aligned with other series that may have
        other index dimensions (usually a serise of output files).

        :param other_count: other series to align count index with
        :rtype: Tupe[pd.DataFrame]
        """
        log_num_align, other_align = self.num_night.align(other, fill_value=0)

        # Filling might switch to float: go back to int (it's a count)
        log_num_align = log_num_align.astype(int)
        other_align = other_align.astype(int)

        return (log_num_align, other_align)
