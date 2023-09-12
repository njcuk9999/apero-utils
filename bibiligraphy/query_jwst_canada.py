import ads
import numpy as np
from astropy.table import Table
import pandas as pd
import openpyxl

# can be either a single year or a range of years with a dash
year_search = '2023'
# If it's a JWST paper, then we'll have one of these keywords in the abstract
abstract_search = 'NIRISS','NIRSpec','MIRI','NIRCAM','JWST'

# maximum number of rows to return, increase if you think you're missing papers
# You'll get a warning if the code reaches this limit anyway
Nmax_rows = 5000

# rejection keys if present in affiliation
rejection_keys = ['Villanueva','Canada-France-Hawaii']

# token for etienne.artigau@umontreal.ca ---> do not get too greedy with my token!
ads.config.token = 'MOMaLYC2B6SA4nhNvgMnHeXhhNTqNkzX21ujXniI'

outname='JWST_CANADA_{}_{}.xlsx'.format(year_search,'-'.join(abstract_search))

# You probably should not edit below this line
#-------------------------------------------------------------------------------

# properly format query
query = "abstract:({}) year:({})".format(' OR '.join(abstract_search),year_search)
print('Query to be passed : {} '.format(query))


papers = ads.SearchQuery(q=query, rows=Nmax_rows, sort="date", fl = ['title','author',
                                                                     'abstract','aff','alternate_bibcode',
                                                                     'pub','pubdate','citation','database','bibcode'])
# convert to a list
papers = list(papers)

# track the number of Canadian papers
N_Canada = 0

# Placeholder for the list of Canadian authors
Canadian_authors = []
Canadian_institutions = []


all_authors = []
all_publishers = []
all_pubdates = []
all_abstracts = []
all_titles = []
all_affiliations = []
all_canadian_authors = []
all_canadian_affiliations = []
all_urls = []

empty_url = '=HYPERLINK("https://ui.adsabs.harvard.edu/abs/{0}/abstract","{0}")'

# You (most likely) have too many papers
if len(papers) == Nmax_rows:
    print('WARNING: maximum number of rows ({}) reached. Some papers may be missing.'.format(Nmax_rows))
    print('Increase Nmax_rows in query_jwst_canada.py')
else:
    # loop over papers, we do not know yet if this is a canadian paper
    for paper in papers:
        # sanity check for occurences where no affiliation is given
        if type(paper.aff) == type(None):
            continue

        for i in range(len(paper.aff)):
            for key in rejection_keys:
                if key in paper.aff[i]:
                    paper.aff[i] = ''


        # Does it include a Canadian author?
        if 'CANADA' not in ''.join(paper.aff).upper():
            continue

        # Does it include a Canadian author?
        N_Canada += 1
        # print the paper title and author info
        print('\n')
        print('-'*80)
        print('Title: {}'.format(paper.title[0]))
        print('\n')
        print('Abstract: {}'.format(paper.abstract))
        print('\n')
        print('All authors: {}'.format(', '.join(paper.author)))
        print('Publisher : {}'.format(paper.pub))
        print('\n')
        print('Authors with Canadian affiliations : ')

        can_author_pub = []
        can_aff_pub = []
        for i in range(len(paper.author)):
            affs = paper.aff[i]
            # sometimes there are multiple affiliations for a single author
            affs = affs.split(';')
            for aff in affs:
                if ('CANADA' in aff.upper()):
                    print('\t{}\n\t {}'.format(paper.author[i],aff))

                    can_author_pub.append(paper.author[i])
                    can_aff_pub.append(aff)

                    # add to list of Canadian authors
                    if paper.author[i] not in Canadian_authors:
                        Canadian_authors.append(paper.author[i])
                        # add to list of Canadian institutions
                        Canadian_institutions.append(aff)

            can_author_pub = list(np.unique(can_author_pub))

        # append relevant info to lists
        # use ' • ' as a separator for multiple authors and affiliations
        all_authors.append(' • '.join(paper.author))
        all_publishers.append(paper.pub)
        all_pubdates.append(paper.pubdate)
        all_abstracts.append(paper.abstract)
        all_titles.append(paper.title[0])
        all_affiliations.append(' • '.join(paper.aff))
        all_canadian_authors.append(' • '.join(can_author_pub))
        all_canadian_affiliations.append(' • '.join(np.unique(can_aff_pub)))
        all_urls.append(empty_url.format(paper.bibcode))

# create a pandas dataframe with all the info
df = pd.DataFrame({'authors': all_authors, 'all_publishers': all_publishers, 'all_pubdates': all_pubdates,'url':all_urls,
                   'all_abstracts':all_abstracts, 'all_titles': all_titles, 'all_affiliations': all_affiliations,
                   'all_canadian_authors': all_canadian_authors, 'all_canadian_affiliations': all_canadian_affiliations})
# save to excel file, first create a sheet with the dataframe
sheetname = 'mySheet'
with pd.ExcelWriter(outname) as writer:
    if not df.index.name:
        df.index.name = 'Index'
    df.to_excel(writer, sheet_name=sheetname)
# then add a table to the sheet
# open the workbook
wb = openpyxl.load_workbook(filename = outname)
tab = openpyxl.worksheet.table.Table(displayName="df", ref=f'A1:{chr(len(df.columns)+64)}{len(df)+1}')
wb[sheetname].add_table(tab)
wb.save(outname)

# order the name of authors alphabetically
Canadian_authors.sort()
# convert to numpy array
Canadian_authors = np.array(Canadian_authors)
# order the name of institutions alphabetically
Canadian_institutions = np.array(Canadian_institutions)
#
print('Canadian authors [n={}] : '.format(len(Canadian_authors)))
print('\n'.join(Canadian_authors))
Canadian_institutions = np.unique(Canadian_institutions)
print('Canadian institutions [n={}] : '.format(len(Canadian_institutions)))
print('\n'.join(Canadian_institutions))

print('\n')
print('Total number of papers : {}'.format(len(papers)))
print('Number of Canadian papers : {}'.format(N_Canada))
