"""
Check if bad pixel calib worked fine.

Tests preformed
    check1: how many recipes were run (cal_badpix_{instrument} in log.fits)?
            how many in the master directory?
    check2: how many of each output do we have?
            #output1: {ODOMETER_CODE}_pp_badpixel.fits
            #output2: {ODOMETER_CODE}_pp_bmap.fits
    check3: how many of each unique output do we have?
            #output1: {ODOMETER_CODE}_pp_badpixel.fits
            #output2: {ODOMETER_CODE}_pp_bmap.fits
    stop1: check3 == check1?
    check4: using the log.fits how many entry failed one or more QC?
            Which nights? Which QC?
    check5: using the log.fits how many entry failed to finish? Which nights?
            Why (using the ERRORS and LOGFILE columns)?
    check6: how many entry BADPIX and BKGRDMAP in
            master_calib_{INSTRUMENT}.txt?
    check7: for each calib entry how many are in the calibDB?
    stop2: check7 == check6?

@author: charles
"""
import os
from typing import List, Optional, Union

import glob
import numpy as np
import pandas as pd

from .apero_tests import Test, CalibTest
import apero_tests_func as atf


class BadPixTest(CalibTest):
    """BadPixTest."""

    def __init__(self,
                 inst: str = 'SPIROU',
                 setup: Optional[str] = None,
                 logdir: Union[str, list] = 'night'):
        """__init__.

        :param inst:
        :type inst: str
        :param setup:
        :type setup: Optional[str]
        :param logdir:
        :type logdir: Union[str, list]
        """

        super().__init__(inst=inst, setup=setup, logdir=logdir)

        self._name = 'Bad Pixel Correction Recipe Test #1'

        self._output_list = ['*_pp_badpixel.fits', '*_pp_bmap.fits']

        # calibDB entries list
        self._calibdb_list = ['BADPIX', 'BKGRDMAP']
        # private pd.Series of output files
        self._output_files = self._gen_output_files()

        self._recipe = 'cal_badpix_{}'.format(self.instrument.lower())

        # Handle logs (this creates log_df and missing_logfits properties)
        self._log_df, self._missing_logs = self._gen_log_df()

    @property
    def name(self):
        """name."""
        return self._name

    @property
    def output_list(self) -> List[str]:
        """List of output string patterns

        :return: output_list
        :rtype: list[str]
        """
        return self._output_list

    @property
    def recipe(self) -> List[str]:
        """Recipe name

        :return: output_list
        :rtype: list[str]
        """
        return self._recipe

    @property
    def calibdb_list(self) -> List[str]:
        """List of calibDB entries

        :return: calibdb_list
        :rtype: list[str]
        """
        return self._calibdb_list

    @property
    def output_files(self) -> pd.Series:
        """output_files.

        :return: output_files
        :rtype: pd.Series
        """
        return self._output_files

    def runtest(self):
        """runtest."""

        print(self.name)
        # output list

        # Total number of log entries for this recipe, and per night
        night_recipe_num_logfits = self.log_df.groupby('DIRECTORY').size()
        tot_recipe_num_logfits = self.log_df.shape[0]  # Check 1
        log_recipe_num = self.log_df.groupby('DIRECTORY').size()

        # Number of outputs for each night for each output pattern
        ind_name = ['PATTERN', 'DIRECTORY']
        output_files_num = self.output_files.groupby(ind_name).size()

        # Align count from log and output (duplicate log on patterns and add
        # missing nights with 0)
        log_num_align, output_num_align = log_recipe_num.align(output_files_num,
                                                               fill_value=0)
        # Filling might switch to float
        log_num_align = log_num_align.astype(int)
        output_num_align = output_num_align.astype(int)

        out_missing_mask = (log_num_align > 0) & (output_num_align == 0)
        output_missing = output_num_align.index.to_frame()[out_missing_mask]
        output_missing = output_missing.reset_index(drop=True)

        # Checking duplicates
        output_dup_mask = log_num_align > output_num_align
        if output_dup_mask.any():

            output_dup_ind = output_dup_mask[output_dup_mask]
            output_dup_ind = output_dup_ind.index .get_level_values(
                                                                'DIRECTORY'
                                                                ).unique()
            # Check if '--master' in runstring
            master_mask = self.log_df.loc[output_dup_ind,
                                          'RUNSTRING'
                                          ].str.contains('--master')
            master_nights = master_mask[master_mask].index.values
            master_recipe_num_logfits = master_nights.size
            if master_recipe_num_logfits > 0:
                comments_check1 = ('An additional {0} recipe with '
                        '--master=True was called in the master '
                        'directory {1}.'.format(master_recipe_num_logfits,
                                                master_nights)
                        )

            # Check for non-master duplicates for each pattern
            # NOTE: This replaces [night_]output{N}_dup and
            master_mask_group = master_mask.astype(int).groupby('DIRECTORY')
            master_num_night = master_mask_group.sum()
            dup_num_night = master_mask_group.size()
            dup_not_master = dup_num_night - master_num_night
            log_dup_align, output_dup_align = dup_not_master.align(
                                            output_num_align[output_dup_mask]
                                            )
            true_dup = log_dup_align - output_dup_align  # number of dups
            true_dup = true_dup[true_dup > 0]  # Keep only non-zero
            true_dup.name = 'COUNT'
            trued_dup = true_dup.reset_index()

        # Check passed QC and ENDED, can then access easily later
        # NOTE: both are empty if all when well
        log_qc_failed = self.log_df[~self.log_df.PASSED_ALL_QC]
        log_ended_false = self.log_df[~self.log_df.ENDED]
        num_logfits_qc_failed = log_qc_failed.size
        num_logfits_ended_false = log_ended_false.size


        # Get number of files in:
        # - logfits (minus master number), check 1
        # - Number of output (all and unique) for each pattern, check 2 and 3
        recipe_num_logfits = tot_recipe_num_logfits - master_recipe_num_logfits
        output_num_all = self.output_files.groupby('PATTERN').size()
        output_base = self.output_files.apply(os.path.basename)
        output_num_unique = output_base.groupby('PATTERN').nunique()

        # check4
        # Report how many and which QC failed
        if num_logfits_qc_failed == 0:
            comments_check4 = ''
            inspect_check4 = ''
        else:
            comments_check4 = 'One or more recipe have failed QC.'
            data_dict_check4 = {
                    'Night': log_qc_failed.reset_index().DIRECTORY.tolist,
                    'QC_STRING': log_qc_failed.QC_STRING,
                                }
            inspect_check4 = atf.inspect_table(
                    'badpixel_test1',
                    'check4',
                    data_dict_check4,
                    'Nights that Failed Quality Control'
                    )

        # check5
        # Report how many and which did not end
        if num_logfits_ended_false == 0:
            comments_check5 = ''
            inspect_check5 = ''
        else:
            comments_check5 = 'One or more recipe have failed to finish.'
            data_dict_check5 = {
                    'Night': log_ended_false.reset_index().DIRECTORY.tolist,
                    'ERRORS': log_ended_false.ERRORS,
                    'LOGFILE': log_ended_false.LOGFILE,
                                }
            inspect_check5 = atf.inspect_table(
                    'badpixel_test1',
                    'check5',
                    data_dict_check5,
                    'Nights that Failed to Finish'
                    )


        # Check if number of unique outputs equals number in logs
        if (output_num_unique == recipe_num_logfits).all():
            color_stop1 = 'Lime'
            result_stop1 = 'Yes'
            comment_stop1 = ''
            inspect_stop1 = ''
        elif (output_num_unique < recipe_num_logfits).any():

            color_stop1 = 'Yellow'
            result_stop1 = 'No'

            # if missing output
            # NOTE: could there be both missing output and duplicate?
            if output_missing.shape[0] > 0:

                comment_stop1 = 'One or more output were not produced.'
                # NOTE: might need arrays instead of list
                data_dict_stop1 = {
                        'Night': output.DIRECTORY.tolist(),
                        'File name': output.PATTERN.tolist(),
                                   }
                inspect_stop1 = atf.inspect_table(
                        'badpixel_test1',
                        'stop1',
                        data_dict_stop1,
                        'Missing Outputs in {0}'.format(reduced_path)
                        )
            # if duplicates
            else:
                comment_stop1 = ('Recipe called multiple times producing the '
                                 'same outputs.')
                data_dict_stop1 = {
                        'Night': true_dup.DIRECTORY.tolist(),
                        'File name': true_dup.PATTERN.tolist(),
                        'Occurrence': true_dup.COUNT.tolist(),
                                   }

                inspect_msg = ('Same Bad Pixel Correction Recipe Called Twice '
                               'or More Producing the Same Outputs in {0}'
                               ).format(reduced_path)
                inspect_stop1 = atf.inspect_table(
                        'badpixel_test1',
                        'stop1',
                        data_dict_stop1,
                        inspect_msg
                        )

        else:
            color_stop1 = 'Red'
            result_stop1 = 'No'
            comment_stop1 = ('The number of unique output files should always '
                             'be smaller or equal to the number of recipe '
                             'called.')
            inspect_stop1 = ''

        # =====================================================================
        # Inspect calibDB
        # =====================================================================
        # Read master_calib_{instument}.txt
        master_calib_path = os.path.join(  # Full path: dir and file
                                self.calibdb_path,
                                'master_calib {}.txt'.format(self.instrument)
                                )
        f = open(master_calib_path, "r")
        master_calib_txt = f.read()

        # Find where DRS processed starts
        nprocessed = master_calib_txt.index('# DRS processed')
        index_start = master_calib_txt[:nprocessed].count('\n')

        # Load dataframe and keep only badpix entries
        master_calib_df = pd.read_csv(
                master_calib_path,
                header=index_start-1,  # Keep first entry with -1
                delimiter=' ',
                usecols=[0, 1, 2, 3],
                names=['KEY', 'MASTER', 'NIGHT', 'FILE'],
                index_col=0,  # Use KEY as index
                )
        master_calib_df.MASTER = master_calib_df.MASTER.astype(bool)
        master_calib_df = master_calib_df.loc[self.calibdb_list]

        # Get total number of entries and of files per output in calibd
        tot_num_entry = master_calib_df.groupby('KEY').size()  # check 6
        master_calib_mask = master_calib_df.MASTER
        master_calib_group = master_calib_df[master_calib_mask].groupby('KEY')
        master_num_entry = master_calib_group.size()
        output_num_entry = tot_num_entry - master_num_entry
        out_patterns = pd.Series(self.output_list)
        sep = os.path.sep
        calib_ind = pd.Index(self.calibdb_list, name='KEY')
        output_calib_paths = pd.Series(self.calibdb_path + sep + out_patterns,
                                       index=calib_ind
                                       )
        output_calibdb = output_calib_paths.apply(glob.glob).explode()
        output_num_calibdb = output_calibdb.groupby('KEY').size()  # check 7

        # Generate comment on extra master entries
        comments_check6 = ('An additional {0} {1} and {2} {3} with master = 1 '
                           'are in the master_calibDB_{4}.'.format(
                               master_num_entry[self.calibdb_list[0]],
                               self._calibdb_list[0],
                               master_num_entry[self.calibdb_list[1]],
                               self._calibdb_list[1],
                               self.instrument)
                           )


        # Checking for missing output for any pattern
        # i.e. files in db entries but not output
        output_calibdb_nopath = output_calibdb.apply(os.path.basename)
        calib_missing_mask = ~master_calib_df.FILE.isin(output_calibdb_nopath)
        if calib_missing_mask.any():

            # NOTE: replaces {night,day}_output{1,2}_missing
            missing_calibdb_output = master_calib_df[calib_missing_mask]

        # Check if some files have duplicaltes in db
        if (output_num_calibdb < output_num_entry).any():

            # Get all duplicates
            calib_dup_mask = master_calib_df.duplicated(keep=False)
            master_mask = master_calib_df.MASTER
            calib_dup = master_calib_df[calib_dup_mask & ~master_mask]
            calib_dup = calib_dup.set_index('FILE', append=True)
            ind_names = calib_dup.index.names  # get file and key to count
            calib_dup['COUNT'] = calib_dup.groupby(ind_names).size()
            calib_dup = calib_dup.reset_index('FILE')
            calib_dup = calib_dup.drop_duplicates()

        # stop2
        if (output_num_calibdb == output_num_entry).all():
            color_stop2 = 'Lime'
            result_stop2 = 'Yes'
            comment_stop2 = ''
            inspect_stop2 = ''

        elif (output_num_calibdb < output_num_entry).any():

            color_stop2 = 'Yellow'
            result_stop2 = 'No'

            # if missing output
            # NOTE: Could we have both missing and dup?
            if calib_missing_mask.any():

                comment_stop2 = 'One or more outputs are not in the calibDB.'
                data_dict_stop2 = {
                        'Night': missing_calibdb_output.NIGHT.tolist(),
                        'File name': missing_calibdb_output.FILE.tolist(),
                                   }
                inspect_stop2 = atf.inspect_table(
                                'badpixel_test1',
                                'stop2',
                                data_dict_stop2,
                                'Missing Output in {0}'.format(calibDB_path)
                                )
            # if duplicates
            else:

                comment_stop2 = ('Some entries in master_calib_{0}.txt are '
                                 'identical.').format(self.instrument)
                data_dict_stop2 = {
                        'Night': calib_dup.NIGHT.tolist(),
                        'File name': calib_dup.FILE.tolist(),
                        'Occurrence': calib_dup.COUNT.tolist(),
                                   }
                inspect_stop2 = atf.inspect_table(
                        'badpixel_test1',
                        'stop2',
                        data_dict_stop2,
                        ('Duplicate Entries in '
                         'master_calib_{0}.txt').format(self.instrument)
                        )

        else:
            color_stop2 = 'Red'
            result_stop2 = 'No'
            comment_stop2 = ('The calibDB should not have more output files '
                             'than what was produced.')
            inspect_stop2 = ''


        # Build badpixel_test1.html
        html_text = f"""
        <html>


        <head>
        <title>APERO Tests</title>
        <style>
        table {{
          width:75%;
        }}
        table, th, td {{
          border: 1px solid black;
          border-collapse: collapse;
        }}
        th, td {{
          padding: 15px;
          text-align: left;
        }}
        #t01 tr:nth-child(even) {{
          background-color: #eee;
        }}
        #t01 tr:nth-child(odd) {{
         background-color: #fff;
        }}
        #t01 th {{
          background-color: white;
          color: black;
        }}
        </style>
        </head>

        <body>

        <h3>Bad Pixel Correction Recipe Test #1</h3>
        <p><b>Setup: {setup}</b><br>
        <p><b>Instrument: {instrument}</b><br>
        <p><b>Date: {date}</b><br>
        <br>
        <p>Script: cal_badpix_{instrument.lower()}.py<br>
        <p>Output files: {', '.join(output_list)}<br>
        <p>Calibration database entry: {', '.join(calibDB_entry_list)}<br>
        <p><a href='https://github.com/njcuk9999/apero-drs#83-bad-pixel-correction-recipe'>Link</a> to Bad Pixel Correction Recipe description</p>
        <br></br>

        <table id="t01">

          <colgroup>
             <col span="1" style="width: 5%;">
             <col span="1" style="width: 55%;">
             <col span="1" style="width: 15%;">
             <col span="1" style="width: 20%;">
             <col span="1" style="width: 5%;">
          </colgroup>
          <tr>
            <th>Check</th>
            <th>Description</th>
            <th>Result</th>
            <th>Comments</th>
            <th>Details</th>
          </tr>
          <tr>
            <td>1</td>
            <td># of time cal_badpix_{instrument.lower()}.py was called</td>
            <td>{recipe_num_logfits}</td>
            <td>{comments_check1}</td>
            <td></td>
          </tr>
          <tr>
            <td>2</td>
            <td># of outputs in {reduced_path}</td>
            <td>{output_list[0]}: {output1_num}<br>{output_list[1]}: {output2_num}</td>
            <td></td>
            <td></td>
          </tr>
          <tr>
            <td>3</td>
            <td># of unique outputs in {reduced_path}</td>
            <td>{output_list[0]}: {output1_num_unique}<br>{output_list[1]}: {output2_num_unique}</td>
            <td></td>
            <td></td>
          </tr>
          <tr>
            <td></td>
            <td>Check 3 == Check 1?</td>
            <td bgcolor={color_stop1}>{result_stop1}</td>
            <td>{comment_stop1}</td>
            <td>{inspect_stop1}</td>
          </tr>
          <tr>
            <td>4</td>
            <td># of entry in {reduced_path}/log.fits that failed one or more QC</td>
            <td>{num_logfits_qc_failed}</td>
            <td>{comments_check4}</td>
            <td></td>
          </tr>
          <tr>
            <td>5</td>
            <td># of entry in {reduced_path}/log.fits that failed to finish</td>
            <td>{num_logfits_ended_false}</td>
            <td>{comments_check5}</td>
            <td></td>
          </tr>
          <tr>
            <td>6</td>
            <td># of {', '.join(calibDB_entry_list)} entry in {calibDB_path}/master_calib_{instrument}.txt</td>
            <td>{calibDB_entry_list[0]}: {output1_num_entry}<br>{calibDB_entry_list[1]}: {output2_num_entry}</td>
            <td>{comments_check6}</td>
            <td></td>
          </tr>
          <tr>
            <td>7</td>
            <td># of outputs in {calibDB_path}</td>
            <td>{output_list[0]}: {output1_num_calibDB}<br>{output_list[1]}: {output2_num_calibDB}</td>
            <td></td>
            <td></td>
          </tr>
          <tr>
            <td></td>
            <td>Check 7 == Check 6?</td>
            <td bgcolor={color_stop2}>{result_stop2}</td>
            <td>{comment_stop2}</td>
            <td>{inspect_stop2}</td>
          </tr>
        </table>


        </body>
        </html>
        """

        with open('badpixel_test1/badpixel_test1.html', 'w') as f:
            f.write(html_text)


if __name__ == '__main__':
    test = BadPixTest()
    test.runtest()
