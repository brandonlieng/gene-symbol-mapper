# gene-symbol-mapper

This script annotates an Excel file containing a column of gene symbols (e.g. IDH1 or ALDH2). Genes are cross-referenced to external databases using flat-file-formatted data from HGNC and Reactome. These flat files are required for the functioning of the script and are included in this repository. With a CSV file of pathways from Reactome (also included here), one can specify which pathways should be reported in the annotations.

To use this script, download the contents of this repository and from a terminal enter:

`python3 gene_symbol_mapper.py -i path_to_input_excel -o path_to_output_excel -n path_to_pathway_csv_file`

## References

Thank you to the HGNC and Reactome for providing access to their data files.

Marc Gillespie, Bijay Jassal, Ralf Stephan, Marija Milacic, Karen Rothfels, Andrea Senff-Ribeiro, Johannes Griss, Cristoffer Sevilla, Lisa Matthews, Chuqiao Gong, Chuan Deng, Thawfeek Varusai, Eliot Ragueneau, Yusra Haider, Bruce May, Veronica Shamovsky, Joel Weiser, Timothy Brunson, Nasim Sanati, Liam Beckman, Xiang Shao, Antonio Fabregat, Konstantinos Sidiropoulos, Julieth Murillo, Guilherme Viteri, Justin Cook, Solomon Shorser, Gary Bader, Emek Demir, Chris Sander, Robin Haw, Guanming Wu, Lincoln Stein, Henning Hermjakob, Peter D’Eustachio, The reactome pathway knowledgebase 2022, *Nucleic Acids Research*, 2021;, gkab1028, https://doi.org/10.1093/nar/gkab1028

Tweedie S, Braschi B, Gray KA, Jones TEM, Seal RL, Yates B, Bruford EA. **Genenames.org: the HGNC and VGNC resources in 2021**. Nucleic Acids Res. PMID: 33152070 PMCID: PMC7779007 DOI: 10.1093/nar/gkaa980
