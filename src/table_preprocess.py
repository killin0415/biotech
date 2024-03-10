import pandas as pd
import fastaparser
import re


def preprocess(data_path: str, saved_path: str, ac: str, length: str) -> None:
    data = {"NCBI GeneID":[], ac: [], length: []}
    with open(data_path) as fasta_file:
        parser = fastaparser.Reader(fasta_file, parse_method='quick')
        for seq in parser:
            par = re.search(r"(\w\w_\d+.\d*).+(\[\w+=\d+\])", seq.header)
            acr = par.group(1)
            idr = par.group(2)
            gene_id = re.search(r"(\d+)", idr).group(1)
            data["NCBI GeneID"].append(gene_id)
            data[ac].append(acr)
            data[length].append(len(seq.sequence))
    df = pd.DataFrame(data)
    df.to_csv(saved_path, index=False)

def matchProteinNRna(data_path: str, saved_path: str) -> None:
    data = {"Protein Accession": [], "Transcript Accession": []}
    with open(data_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if re.match(r".*CDS\t\d+\t\d+.*", line):
                pair = re.search(r"(ID=.+?;)(Parent=.+?;)", line)
                protein = pair.group(1)[7:-1]
                rna = pair.group(2)[11:-1]
                data["Protein Accession"].append(protein)
                data["Transcript Accession"].append(rna)
    df = pd.DataFrame(data)
    df.to_csv(saved_path, index=False)    


if __name__ == "__main__":
    preprocess("input\\Dsim\\protein.faa", "results\\protein_length_table.csv", "Protein Accession", "Protein_length")
    preprocess("input\\Dsim\\rna.fna", "results\\transcript_length_table.csv", "Transcript Accession", "Transcript_length")
    matchProteinNRna("input\\Dsim\\genomic.gff", "results\\protein_transcript_pair.csv")



