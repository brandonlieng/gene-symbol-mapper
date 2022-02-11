import argparse
import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm


def get_gene_info(gene_symbols, hgnc_data):
    gene_info = [["name", "ensembl_gene_id",
                  "uniprot_ids","enzyme_id"]]

    # Would be faster as a vectorized operation but we need to search through
    # both the current and previous symbol columns for potential matches
    for symbol in tqdm(gene_symbols):
        result_row = hgnc_data[hgnc_data["symbol"] == symbol]
        if result_row.empty:
            search_aliases = [symbol in aliases.split("|") if not \
                              pd.isna(aliases) else False \
                              for aliases in hgnc_data["alias_symbol"]]

            # Aliases can be assigned to more than one gene
            if any(search_aliases) and sum(search_aliases) == 1:
                result_row = hgnc_data[search_aliases]

        if not result_row.empty:
            assert result_row.shape[0] == 1, "More than one result for symbol"
            gene_info.append([result_row["name"].values[0],
                               result_row["ensembl_gene_id"].values[0],
                               result_row["uniprot_ids"].values[0],
                               result_row["enzyme_id"].values[0]])
        else:
            gene_info.append([np.nan] * 4)

    return gene_info


def get_pathway_mappings(ensembl_gene_ids, reactome_data, pathways_to_report):
    pathway_mappings = [["ensembl_gene_id", "reactome_ids", "event_names"]]

    if pathways_to_report is not None:
        included_pathways = pathways_to_report[pathways_to_report["Include"]]
        collapsed_pathway_ids = dict(zip(included_pathways["PathwayID"],
                                         included_pathways["Collapse"]))
        collapsed_pathway_ids = {k: k if pd.isnull(collapsed_pathway_ids[k])\
                                 else collapsed_pathway_ids[k]\
                                 for k in collapsed_pathway_ids}
        pathway_id_name_pairs = dict(zip(included_pathways["PathwayID"],
                                         included_pathways["PathwayName"]))

    for ensembl_gene_id in tqdm(ensembl_gene_ids):
        if pd.isna(ensembl_gene_id):
            pathway_mappings.append([np.nan] * 3)
            continue

        gene_pathway_mappings = reactome_data[reactome_data["ensembl_gene_id"] \
                                              == ensembl_gene_id]

        # Remove pathways that are not part of the inclusion list
        gene_pathway_mappings = gene_pathway_mappings[ \
            [pathway in included_pathways["PathwayID"].values for pathway in \
             gene_pathway_mappings["reactome_id"]]
        ]

        if not gene_pathway_mappings.empty:
            # Collapse pathways to less-granular levels as per the incl. list
            gene_pathway_mappings["reactome_id"] = \
                gene_pathway_mappings["reactome_id"].map(collapsed_pathway_ids)

            # Fix names
            gene_pathway_mappings["event_name"] = \
                gene_pathway_mappings["reactome_id"].map(pathway_id_name_pairs)

            gene_pathway_mappings.drop_duplicates(subset=["ensembl_gene_id",
                                                          "reactome_id"],
                                                  inplace=True)
            pathway_mappings.append([
                ensembl_gene_id,
                gene_pathway_mappings["reactome_id"].str.cat(sep="|"),
                gene_pathway_mappings["event_name"].str.cat(sep="|")])
        else:
            pathway_mappings.append([np.nan] * 3)

    return pathway_mappings


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, help="file path to the Excel workbook")
    parser.add_argument("-o", type=str, help="desired output file path")
    parser.add_argument("-n", type=str,
                        help="name of the column holding gene symbols")
    parser.add_argument("-p", type=str,
                        help=("file path to the CSV file of Reactome pathway" +
                              " identifiers to include"))
    parser.add_argument("-g", type=str,
                        help="file path to the .txt-formatted HGNC database",
                        default="./hgnc_complete_set.txt")
    parser.add_argument("-m", type=str,
                        help=("file path to the .txt-formatted " +
                              "Ensembl2Reactome mapping file"),
                        default="./Ensembl2Reactome_Homo_sapiens.txt")

    flags = parser.parse_args()

    if None in [flags.i, flags.o, flags.n]:
        parser.print_help()
        exit()

    # Load Excel workbook holding gene symbols to convert/map
    print("Loading Excel workbook...")
    excel_data = pd.read_excel(flags.i)

    # Load HGNC data
    print("Loading HGNC data table...")
    hgnc_data = pd.read_csv(flags.g, sep="\t",
                            usecols=["hgnc_id", "symbol", "name", "locus_type",
                                     "status", "alias_symbol", "prev_symbol",
                                     "ensembl_gene_id", "uniprot_ids",
                                     "enzyme_id"]
                            )

    # Load Reactome pathway mapping file
    print("Loading Reactome data table...")
    reactome_data = pd.read_csv(flags.m, sep="\t",
                                names=["ensembl_gene_id", "reactome_id",
                                       "url", "event_name", "evidence_code",
                                       "species"]
                                )

    # Get gene info for each gene symbol in the Excel workbook
    print("Getting gene information for gene symbols...")
    gene_info = get_gene_info(excel_data[flags.n], hgnc_data)
    gene_info = pd.DataFrame(gene_info[1:], columns=gene_info[0])

    # If an in/exclusion file of pathways is provided, collapse pathway mappings
    # to the level reported, then get pathway mappings for each of the ENSEMBL
    # gene IDs
    print("Getting pathway mappings...")
    pathways_to_report = None
    if flags.p is not None:
        with open("./ReactomePathwaysBranches_14Nov2021.pkl", "rb") as b_file:
            pathway_branches = pickle.load(b_file)
        with open(flags.p, "r") as inc_fp:
            pathways_to_report = pd.read_csv(inc_fp)

    pathway_mappings = get_pathway_mappings(gene_info["ensembl_gene_id"],
                                            reactome_data,
                                            pathways_to_report)
    pathway_mappings = pd.DataFrame(pathway_mappings[1:],
                                    columns=pathway_mappings[0])

    # Consolidate the cross-referenced info into one dataframe and insert
    # the columns after the gene symbol column in the workbook data
    all_crossref_info = pd.concat([gene_info, pathway_mappings], axis=1)
    gene_symbols_col_idx = excel_data.columns.values.tolist().index(flags.n)
    output_data = pd.concat([excel_data.iloc[:, 0:gene_symbols_col_idx + 1],
                             all_crossref_info,
                             excel_data.iloc[:, gene_symbols_col_idx + 1:]],
                             axis=1
                            )

    # Save the output workbook
    print("Writing output Excel workbook...")
    output_data.to_excel(flags.o, index=False)
