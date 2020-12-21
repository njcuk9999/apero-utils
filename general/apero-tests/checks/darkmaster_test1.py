"""
Check if dark master calib worked fine.

Tests performed:
    check1: how many recipes were run? (cal_dark_master_{instrument} in
            other/log.fits)
    check2: how many of each output do we have?
            output1: {ODOMETER_CODE}_pp_dark_master.fits
    check3: how many of each unique output do we have?
            output1: {ODOMETER_CODE}_pp_dark_master.fits
    stop1: check3 == check1?
    check4: using the log.fits how many entry failed one or more QC?
            Which nights? Which QC?
    check5: using the log.fits how many entry failed to finish? Which nights?
            Why (using the ERRORS and LOGFILE columns)?
    check6: how many entry DARKM in master_calib_{INSTRUMENT}.txt?
    check7: for each calib entry how many are in the calibDB?
    stop2: check7 == check6?

@author: charles
"""
import os
from typing import Optional
import numpy as np

from checks.apero_tests import Test
import apero_tests_func as atf


class DarkMTest(Test):
    """DarkMTest."""

    def __init__(self, inst: str = 'SPIROU', setup: Optional[str] = None):
        """__init__.

        :param inst:
        :type inst: str
        :param setup:
        :type setup: Optional[str]
        """

        super().__init__(inst=inst, setup=setup)

        self._name = 'Dark Master Recipe Test #1'

        # output list
        self._output_list = ['*_pp_dark_master.fits']

        # dataframe of all outputs
        self._output_df = self._get_output_df()

        # calibDB entries list
        self._calibdb_path = ['DARKM']

    @property
    def name(self) -> str:
        """name.

        :return: name
        :rtype: str
        """
        return self._name

    @property
    def output_list(self) -> list[str]:
        """List of output string patterns

        :return: output_list
        :rtype: list[str]
        """
        return self._output_list

    @property
    def calibdb_list(self) -> list[str]:
        """List of calibDB entries

        :return: calibdb_list
        :rtype: list[str]
        """
        return self._calibdb_path

    def runtest(self):
        """runtest."""

        # inspect all reduced_nights and other/log.fits
        output1 = []  # Check 2
        nights_output1 = []

        # =====================================================================
        # Gather all outputs in reduced night dirs
        # =====================================================================
        for i in range(len(reduced_nights)):

            output1_list_files = atf.list_files('{0}/{1}'.format(reduced_path,
                                                reduced_nights[i]),
                                                files=output_list[0][1:])

            output1.extend(output1_list_files)
            nights_output1.extend(reduced_nights[i] * len(output1_list_files))

        # TODO: get from output df
        output1_num = len(output1)  # check2
        output1_num_unique = len(np.unique(output1))  # check3

        # =====================================================================
        # Load logfits info (single file in other)
        # =====================================================================
        if os.path.isfile('{0}/other/log.fits'.format(reduced_path)):
            logfits = atf.log_fits('{0}/other/log.fits'.format(reduced_path))

            index_recipe = logfits.recipe == 'cal_dark_master_{0}'.format(
                                                                    instrument.lower()
                                                                    )
            recipe_num_logfits = sum(index_recipe)  # check1

            # checking for missing output
            if (recipe_num_logfits > output1_num_unique and len(np.unique(logfits.nights[index_recipe])) > len(np.unique(nights_output1))):
                night_output1_missing = logfits.nights[index_recipe][np.in1d(logfits.nights[index_recipe], ~np.unique(nights_output1))]
                output1_missing = output_list * len(night_output1_missing)

            # checking for duplicates
            if recipe_num_logfits > output1_num_unique and len(np.unique(logfits.nights[index_recipe])) == len(np.unique(nights_output1)):

                nights_logfits_unique, return_counts_output1 = np.unique(logfits.nights[index_recipe], return_counts=True)

                night_output1_dup = nights_logfits_unique[return_counts_output1 > 1]
                output1_dup = output_list * len(night_output1_dup)

            indexQCfalse = np.array(logfits.indexQCfalse) * np.array(index_recipe)  # QC false + correct recipe
            indexENDEDfalse = np.array(logfits.indexENDEDfalse) * np.array(index_recipe)  # ENDED false + correct recipe

            num_logfits_QCfalse = sum(indexQCfalse)  # check4
            num_logfits_ENDEDfalse = sum(indexENDEDfalse)  # check5

            nights_logfits_QCfalse = logfits.nights[indexQCfalse]  # check4
            QCstr_logfits_QCfalse = logfits.QCstr[indexQCfalse]    # check4

            nights_logfits_ENDEDfalse = logfits.nights[indexENDEDfalse]    # check5
            ERRORS_logfits_ENDEDfalse = logfits.ERRORS[indexENDEDfalse]    # check5
            LOGFILE_logfits_ENDEDfalse = logfits.LOGFILE[indexENDEDfalse]  # check5

        # missing log.fits
        else:
            missing_logfits = '{0}/other/log.fits'.format(reduced_path)


        # check4
        if num_logfits_QCfalse == 0:
            comments_check4 = ''
            inspect_check4 = ''
        else:
            comments_check4 = 'One or more recipe have failed QC.'
            data_dict_check4 = {'Night': nights_logfits_QCfalse,
                                'QC_STRING': QCstr_logfits_QCfalse,
                                }
            inspect_check4 = atf.inspect_table('darkmaster_test1',
                                               'check4',
                                               data_dict_check4,
                                               'Nights that Failed Quality Control'
                                               )

        # check5
        if num_logfits_ENDEDfalse == 0:
            comments_check5 = ''
            inspect_check5 = ''
        else:
            comments_check5 = 'One or more recipe have failed to finish.'
            data_dict_check5 = {'Night': nights_logfits_ENDEDfalse,
                                'ERRORS': ERRORS_logfits_ENDEDfalse,
                                'LOGFILE': LOGFILE_logfits_ENDEDfalse,
                                }
            inspect_check5 = atf.inspect_table('darkmaster_test1',
                                               'check5',
                                               data_dict_check5,
                                               'Nights that Failed to Finish')


        # stop1
        if output1_num_unique == recipe_num_logfits:
            color_stop1 = 'Lime'
            result_stop1 = 'Yes'
            comment_stop1 = ''
            inspect_stop1 = ''
        elif output1_num_unique < recipe_num_logfits:

            color_stop1 = 'Yellow'
            result_stop1 = 'No'

            # if missing output
            if len(output1_missing) > 0:
                comment_stop1 = 'One or more output were not produced.'
                data_dict_stop1 = {'Night': night_output1_missing,
                                   'File name': output1_missing,
                                   }
                inspect_stop1 = atf.inspect_table('darkmaster_test1',
                                                  'stop1',
                                                  data_dict_stop1,
                                                  'Missing Outputs in {0}'.format(
                                                      reduced_path
                                                      ),
                                                  )
            # if duplicates
            else:
                comment_stop1 = ('Recipe called multiple times producing the same '
                                 'outputs.')
                data_dict_stop1 = {'Night': night_output1_dup,
                                   'File name': output1_dup,
                                   }
                inspect_stop1 = atf.inspect_table(
                        'darkmaster_test1',
                        'stop1',
                        data_dict_stop1,
                        ('Same Dark Master Recipe Called Twice or More '
                         'Producing the Same Outputs in {0}'
                         ).format(reduced_path)
                         )
        else:
            color_stop1 = 'Red'
            result_stop1 = 'No'
            comment_stop1 = ('The number of unique output files should always be '
                             'smaller or equal to the number of recipe called.'
                             )
            inspect_stop1 = ''

        # Inspect calibDB
        f = open("{0}/master_calib_{1}.txt".format(calibDB_path, instrument), "r")
        master_calib_txt = f.read()
        nprocessed = master_calib_txt.index('# DRS processed')
        index_start = master_calib_txt[:nprocessed].count('\n')
        key_col, master_col, night_col, file_col = np.genfromtxt(
                "{0}/master_calib_{1}.txt".format(calibDB_path, instrument),
                delimiter=' ',
                unpack=True,
                usecols=(0, 1, 2, 3),
                skip_header=index_start,
                dtype=str)
        master_col = master_col.astype(dtype = bool) # str to bool

        index_key_output1 = key_col == calibDB_entry_list[0]

        output1_num_entry = len(key_col[index_key_output1]) # check6

        output1_calibDB = atf.list_files("{0}".format(calibDB_path),
                                         files=output_list[0][1:])
        output1_num_calibDB = len(output1_calibDB)  # check7


        # checking for missing output

        output1_missing_calibDB = []

        if sum(np.in1d(file_col[index_key_output1], output1_calibDB)) < len(file_col[index_key_output1]): 

            index_output1_missing = ~np.in1d(
                    file_col[index_key_output1],
                    output1_calibDB)

            night_output1_missing_calibDB = night_col[index_key_output1][index_output1_missing]
            output1_missing_calibDB = file_col[index_key_output1][index_output1_missing]


        # checking for duplicates

        if len(output1_calibDB) < output1_num_entry:

            (file_col_output1_unique,
             index_output1_dup,
             return_counts_output1) = np.unique(
                         file_col[index_key_output1],
                         return_index=True,
                         return_counts=True
                         )

            count_mask1 = return_counts_output1 > 1

            night_output1_dup_calibDB = night_col[index_key_output1][index_output1_dup][count_mask1]
            output1_dup_calibDB = file_col_output1_unique[count_mask1]


        # stop2
        if output1_num_calibDB == output1_num_entry:
            color_stop2 = 'Lime'
            result_stop2 = 'Yes'
            comment_stop2 = ''
            inspect_stop2 = ''

        elif output1_num_calibDB < output1_num_entry:

            color_stop2 = 'Yellow'
            result_stop2 = 'No'

            # if missing output
            if len(output1_missing_calibDB) > 0:

                comment_stop2 = 'One or more outputs are not in the calibDB.'
                data_dict_stop2 = {'Night': night_output1_missing_calibDB,
                                   'File name': output1_missing_calibDB,
                                   }
                inspect_stop2 = atf.inspect_table('darkmaster_test1',
                                                  'stop2',
                                                  data_dict_stop2,
                                                  'Missing Output in {0}'.format(
                                                      calibDB_path)
                                                  )
            # if duplicates
            else:

                comment_stop2 = ('Some entries in master_calib_{0}.txt are '
                                 'identical.').format(instrument)
                data_dict_stop2 = {'Night': night_output1_dup_calibDB,
                                   'File name': output1_dup_calibDB,
                                   'Occurrence': return_counts_output1[count_mask1]
                                   }
                inspect_stop2 = atf.inspect_table(
                        'darkmaster',
                        'stop2',
                        data_dict_stop2,
                        ('Duplicate Entries in '
                         'master_calib_{0}.txt').format(instrument)
                        )

        else:
            color_stop2 = 'Red'
            result_stop2 = 'No'
            comment_stop2 = ('The calibDB should not have more output files than what '
                             'was produced.')
            inspect_stop2 = ''


        # Build darkmaster_test1.html
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

        <h3>Dark Master Recipe Test #1</h3>
        <p><b>Setup: {setup}</b><br>
        <p><b>Instrument: {instrument}</b><br>
        <p><b>Date: {date}</b><br>
        <br>
        <p>Script: cal_dark_master_{instrument.lower()}.py<br>
        <p>Output files: {output_list}<br>
        <p>Calibration database entry: {calibDB_entry_list}<br>
        <p><a href='https://github.com/njcuk9999/apero-drs#82-dark-master-recipe'>Link</a> to Dark Master Recipe description</p>
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
            <td># of time cal_dark_master_{instrument.lower()}.py was called</td>
            <td>{recipe_num_logfits}</td>
            <td></td>
            <td></td>
          </tr>
          <tr>
            <td>2</td>
            <td># of outputs in {reduced_path}</td>
            <td>{output_list[0]}: {output1_num}</td>
            <td></td>
            <td></td>
          </tr>
          <tr>
            <td>3</td>
            <td># of unique outputs in {reduced_path}</td>
            <td>{output_list[0]}: {output1_num_unique}</td>
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
            <td># of entry in {reduced_path}/other/log.fits that failed one or more QC</td>
            <td>{num_logfits_QCfalse}</td>
            <td>{comments_check4}</td>
            <td></td>
          </tr>
          <tr>
            <td>5</td>
            <td># of entry in {reduced_path}/other/log.fits that failed to finish</td>
            <td>{num_logfits_ENDEDfalse}</td>
            <td>{comments_check5}</td>
            <td></td>
          </tr>
          <tr>
            <td>6</td>
            <td># of {calibDB_entry_list[0]} entry in {calibDB_path}/master_calib_{instrument}.txt</td>
            <td>{output1_num_entry}</td>
            <td></td>
            <td></td>
          </tr>
          <tr>
            <td>7</td>
            <td># of outputs in {calibDB_path}</td>
            <td>{output_list[0]}: {output1_num_calibDB}</td>
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

        with open('darkmaster_test1/darkmaster_test1.html', 'w') as f:
            f.write(html_text)


if __name__ == '__main__':
    test = DarkMTest()
    test.runtest()
