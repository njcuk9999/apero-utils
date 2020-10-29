import os


test_list_short = ['preprocessing_test1', 'darkmaster_test1', 'badpixel_test1',
                   'localisation_test1', 'shapemaster_test1', 'shape_test1',
                   'flat_test1', 'thermal_test1', 'masterleak_test1',
                   'leak_test1', 'masterwavelength_test1', 'wavelength_test1',
                   'extraction_test1', 'extraction_test2', 'extraction_test3',
                   'maketellu_test1', 'fittellu_test1', 'maketemplate_test1',
                   'ccf_test1']

test_list_long = ['Preprocessing Recipe Test #1', 'Dark Master Recipe Test #1',
                  'Bad Pixel Correction Recipe Test #1',
                  'Localisation Recipe Test #1', 'Shape Master Recipe Test #1',
                  'Shape (per night) Recipe Test #1',
                  'Flat/Blaze Correction Recipe Test #1',
                  'Thermal Correction Test #1',
                  'Master Leak Correction Recipe Test #1',
                  'Leak (per night) Correction Test #1',
                  'Master Wavelength Solution Recipe Test #1',
                  'Wavelength Solution (per night) Test #1',
                  'Extraction Recipe Test #1', 'Extraction Recipe Test #2',
                  'Extraction Recipe Test #3', 'Make Telluric Recipe Test #1',
                  'Fit Telluric Recipe Test #1',
                  'Make Template Recipe Test #1', 'CCF Recipe Test #1']

n = len(test_list_short)  # number of tests
n = 1
for i in range(n):

    if not os.path.isdir(test_list_short[i]):
        os.system('mkdir {0}'.format(test_list_short[i]))

    print('test {0}/{1}'.format(i+1,n))   
    print('running {0}\n'.format(test_list_short[i]))    
    os.system('python {0}.py'.format(test_list_short[i]))


print('all tests done')   


# build table element

html_table = []
for i in range(n):
    html_table.append("""
  <tr>
    <td><a href='{0}/{0}.html'>{1}</a></td>
  </tr>
""".format(test_list_short[i], test_list_long[i]))
html_table = "".join(html_table)


# build main .html doc

html_text = f"""
<html>

<head>
<title>APERO Tests</title>
<style>
table {{
  width:25%;
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

<img src='images/apero_logo.png' alt='APERO'>

<p>A PipelinE to Reduce Observations</p>

<table id="t01">
  <tr>
    <th>Recipe Test</th> 
  </tr>
{html_table}
</table>

</body>
</html>
"""


with open('apero_tests.html', 'w') as f:
    f.write(html_text)
