"""
Automated NCBI Sequence Retrieval via Entrez E-utilities
Author: Surena Ghadam
Project: Computational Biology Basics (02_hbb_mutation_and_entrez)

Pipeline Workflow:
1. Configure Entrez credentials with user email (NCBI compliance).
2. Automate accession retrieval using E-Search queries.
3. Fetch curated mRNA sequences directly from NCBI Nucleotide database via E-Fetch.
4. Save structured FASTA files into the project data directory.
"""

from pathlib import Path
from Bio import Entrez, SeqIO


ENTREZ_EMAIL = "surenaghadam2007@gmail.com"

def search_ncbi_accession(term: str, db: str = "nucleotide") -> str:
    """
    Search NCBI database using keywords and return the top matching accession.
    """
    Entrez.email = ENTREZ_EMAIL
    with Entrez.esearch(db=db, term=term, retmax=1, sort="relevance") as handle:
        search_results = Entrez.read(handle)

    id_list = search_results.get("IdList", [])
    if not id_list:
        raise ValueError(f"No records found for query: {term}")

    return id_list[0]


def fetch_sequence(uid_or_accession: str, db: str = "nucleotide") -> SeqIO.SeqRecord:
    """
    Fetch a full sequence record in FASTA format using accession or UID.
    """
    Entrez.email = ENTREZ_EMAIL
    with Entrez.efetch(db=db, id=uid_or_accession, rettype="fasta", retmode="text") as handle:
        record = SeqIO.read(handle, "fasta")
    return record


def main() -> None:
    base_dir = Path(__file__).resolve().parent.parent
    data_dir = base_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    # Standard curated targets for hemoglobin beta comparative analysis:
    # 1. Human (Homo sapiens) HBB
    # 2. Mouse (Mus musculus) Hbb-b1 (Gold standard comparative model)
    targets = {
        "human_hbb_mrna.fasta": "NM_000518.5",
        "mouse_hbb_mrna.fasta": "NM_008220.4",
    }

    print("=" * 65)
    print("NCBI ENTREZ CURATED GENOMIC INGESTION")
    print("=" * 65)

    for filename, accession in targets.items():
        try:
            print(f"[*] Fetching target: {accession} -> {filename}...")
            record = fetch_sequence(uid_or_accession=accession)
            output_path = data_dir / filename

            SeqIO.write(record, output_path, "fasta")

            print(f"[+] Download complete: {record.id}")
            print(f"    Description: {record.description}")
            print(f"    Length: {len(record.seq)} nt")
            print(f"    Saved: {output_path.name}\n")

        except Exception as err:
            print(f"[-] Error retrieving {accession}: {err}\n")

    print("=" * 65)
    print("Ingestion pipeline successfully finished.")


if __name__ == "__main__":
    main()
