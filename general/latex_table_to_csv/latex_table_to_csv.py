import re
import csv
import pandas as pd


def latex_longtable_to_csv(latex_file, csv_file):
    with open(latex_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    in_table = False
    rows = []
    in_content = False

    for line in lines:
        line = line.strip()
        if r"\begin{longtable}" in line:
            in_table = True
            continue
        if r"\end{longtable}" in line:
            break

        if line.startswith('\endhead'):
            in_content = True
            continue
        if not in_content:
            continue

        if line.startswith('\\') or line.startswith('"'):
            continue

        if len(line) == 0:
            continue

        if in_table:
            # Remove LaTeX commands and extract table content
            line = re.sub(r"\\(hline|cline\{.*?\})", "", line)
            line = re.sub(r"\\multicolumn\{\d+\}\{.*?\}\{(.*?)\}", r"\1", line)
            line = re.sub(r"\\href\{.*?\}\{(.*?)\}", r"\1", line)  # Keep only Y in \href{X}{Y}
            line = line.replace(r"\\" , "")  # Remove end-of-row markers
            columns = [col.strip() for col in line.split("&")]
            if columns:
                rows.append(columns)

    # Write to CSV
    with open(csv_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerows(rows)

def list_of_tellurics(csv_file, python_file):
    # Don't use first row as column names
    df = pd.read_csv("output.csv", header=None)

    # tstar columns
    cols = [1,2,3,4]
    # storage for telluric stars
    tstars_all = []

    for col in cols:

        tstars_all += list(df[col])

    # make unique
    tstars = list(set(tstars_all))


# Example Usage
latex_longtable_to_csv("table.tex", "output.csv")


list_of_tellurics("output.csv", "output.py")