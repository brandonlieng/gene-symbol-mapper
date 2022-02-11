# gene-symbol-mapper

This script annotates an Excel file containing a column of gene symbols (e.g. IDH1 or ALDH2). Genes are cross-referenced to external databases using flat-file-formatted data from HGNC and Reactome. These flat files are required for the functioning of the script and are included in this repository. With a CSV file of pathways from Reactome (also included here), one can specify which pathways to report in the annotations.

To use this script, download the contents of this repository and from a terminal enter:

`python3 gene_symbol_mapper.py -i path_to_input_excel -o path_to_output_excel -n path_to_pathway_csv_file`
