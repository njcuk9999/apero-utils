#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-07-29 at 10:23

@author: cook
"""
import os

import pandas as pd
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
# path to the xls file
path_to_file = '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero/lang/databases/language.xls'

out_path = '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero/core/lang/tables'

# Define preable
PREAMBLE = """
#!/usr/bin/env python
# -*- coding: utf-8 -*-
\"\"\"
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-07-29 at 09:10

@author: cook
\"\"\"
from apero.base import base
from apero.core.lang import drs_lang_list


# =============================================================================
# Define variables
# =============================================================================
__NAME__ = 'apero.lang.tables.{FILENAME}'
__PACKAGE__ = base.__PACKAGE__
__INSTRUMENT__ = 'None'
__version__ = base.__version__
__author__ = base.__author__
__date__ = base.__date__
__release__ = base.__release__

# =============================================================================
# Define functions
# =============================================================================
# Get the language list
langlist = drs_lang_list.LanguageList()

"""

# Define end of code text
END = """

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # print 'Hello World!'
    print("Hello World!")

# =============================================================================
# End of code
# =============================================================================
"""

# Define how to add a new item
PYCODE = """
# =============================================================================
# {KEYWORD} 
# =============================================================================
item = langlist.create('{KEYWORD}', kind='{KIND}')
{ITEMS}
item.arguments = 'None'
item.comment = '{COMMENT}'
langlist.add(item)"""
# define how to add a language instance
PYVALUE = "item.value['{LANG}'] = '{VALUE}'"

LANGUAGES = ['ENG', 'FR']

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def clean_text(text):
    text = text.replace("'", "\\'")
    text = text.replace('"', '\\"')
    text = text.replace('\n', '\\n')
    text = text.replace('\t', '\\t')
    return text


def clean_kind(text):
    if text in ['HELP', 'TEXT']:
        return text
    else:
        return f'{text.replace(" ", "_").lower()}-code'

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":

    # storage of items strings
    item_strings = dict()
    item_string_dict = dict()
    # ----------------------------------------------------------------------
    # Step 1: Load the excel sheet
    # ----------------------------------------------------------------------
    xls = pd.read_excel(path_to_file, sheet_name=None)

    # ----------------------------------------------------------------------
    # Step 2: Get the sheets from the excel sheet
    # ----------------------------------------------------------------------
    sheets = xls.keys()

    # ----------------------------------------------------------------------
    # Step 3: loop around sheets
    # ----------------------------------------------------------------------
    for sheet in sheets:
        # don't do blank sheets
        if 'BLANK' in sheet:
            continue
        # print progress
        print('Processing sheet {0}'.format(sheet))
        # storage for items
        item_strings[sheet] = []
        item_string_dict[sheet] = dict()
        # convert all columns to strings
        for col in xls[sheet].columns:
            xls[sheet][col] = xls[sheet][col].astype(str)

        # loop around rows
        for row in tqdm(xls[sheet].index):

            # -----------------------------------------------------------------
            # Step 4: Get the keyname, kind, keydesc, arguments, eng, fr
            #     from sheet
            # -----------------------------------------------------------------
            keyname = clean_text(xls[sheet]['KEYNAME'][row])

            # deal with empty
            if keyname.upper() in ['NULL', 'NONE', 'NAN', '']:
                continue

            kind = clean_text(xls[sheet]['KIND'][row])
            kind = clean_kind(kind)

            keydesc = clean_text(xls[sheet]['KEYDESC'][row])
            if keydesc.upper() in ['NULL', 'NONE', 'NAN', '']:
                keydesc = ''

            arguments = clean_text(xls[sheet]['ARGUMENTS'][row])
            if arguments.upper() in ['NULL', 'NONE', 'NAN', '']:
                arguments = 'None'

            # deal with languages
            item_dict = dict()
            # loop around languages
            for lang in LANGUAGES:
                # get value
                value = clean_text(xls[sheet][lang].iloc[row])
                # deal with empty
                if value.upper() in ['NULL', 'NONE', 'NAN', '']:
                    continue
                # push into storage
                item_dict[lang] = value


            # -----------------------------------------------------------------
            # Step 5: Construct the items language string
            # -----------------------------------------------------------------
            langs = []
            for lang in LANGUAGES:
                if lang not in item_dict:
                    continue
                langs.append(PYVALUE.format(LANG=lang, VALUE=item_dict[lang]))

            if len(langs) == 0:
                continue

            # -----------------------------------------------------------------
            # Step 6: Construct the items string
            # -----------------------------------------------------------------
            keydict = dict(KEYWORD=keyname, KIND=kind, ITEMS='\n'.join(langs),
                           COMMENT=keydesc)

            item = PYCODE.format(**keydict)




            # deal with code already in another sheet
            used = False

            if sheet not in ['HELP', 'TEXT']:
                for master in ['HELP', 'TEXT']:
                    if keyname in item_string_dict[master]:
                        if item == item_string_dict[master][keyname]:
                            used = True

            if not used:
                # append to list
                item_strings[sheet].append(item)
                item_string_dict[sheet][keyname] = item

        # construct filename
        if sheet == 'HELP':
            filename = 'default_help.py'
        elif sheet == 'TEXT':
            filename = 'default_text.py'
        elif 'TEXT_' in sheet:
            instrument = sheet.lower().replace('text_', '')
            filename = f'{instrument}_text.py'
        else:
            instrument = sheet.lower().replace('help_', '')
            filename = f'{instrument}_help.py'

        # get all lines
        lines = PREAMBLE.format(FILENAME=filename)
        lines += '\n'.join(item_strings[sheet])
        lines += END
        # construct absolute path
        path = os.path.join(out_path, filename)
        # save this sheets items to a file
        with open(path, 'w') as file:
            file.write(lines)



# =============================================================================
# End of code
# =============================================================================
