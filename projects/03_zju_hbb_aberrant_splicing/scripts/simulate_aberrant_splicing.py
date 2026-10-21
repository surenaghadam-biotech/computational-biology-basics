"""
Aberrant Splicing Simulation for Human HBB IVS-II-654 (C>T) Mutation.

This module reproduces the molecular pathology described in:
Yang et al. (2023), "Therapeutic genome editing of an aberrant splice site in
beta-thalassemia by CRISPR/Cas9 with multiple sgRNAs", Genes & Diseases (PMC10425742).

Biological Mechanism:
    - Pathogenic variant: IVS-II-654 (C>T) in intron 2 of human HBB gene.
    - Activates a cryptic splice donor site and unleashes an upstream cryptic
      acceptor site at nucleotide 579 of intron 2.
    - Causes retention of a 73-nt pseudo-exon within mature beta-globin mRNA:
      mRNA_mutant = Exon_1 + Exon_2 + Pseudo_Exon (73 nt) + Exon_3.
    - The 73-nt insertion (+1 frameshift) induces a Premature Termination
      Codon (PTC), causing severe beta-thalassemia major phenotype.

Author: Surena Ghadam
Project: Computational Biology Research - Target: Zhejiang University (ZJU)
Target Accession: NCBI GenBank U01317.1 (Human beta-globin region)
"""

from Bio import Entrez,SeqIO
from Bio.Seq import Seq
# --- 1. Fetching Genomic Data from NCBI ---
def Entrez_gb(id) :
  Entrez.email="surenaghadam2007@gmail.com"
  with Entrez.efetch(db="nucleotide",id=id,rettype="gb",retmode="text") as handel :
    record=SeqIO.read(handel,"genbank")
  return record
# --- 2. Extracting Exons and Introns ---
def exons_extraxt(record):
  """Extract HBB gene exons with genomic coordinates."""
  exons=[]
  for feature in record.features:
    # Filter strictly for HBB gene exons to avoid beta-globin cluster interference
    if feature.type=="exon" and feature.qualifiers.get("gene")==["HBB"] :
      start=int(feature.location.start)
      end=int(feature.location.end)
      exon_seq=feature.extract(record.seq)
      exons.append({"start": start , "end":end ,"exon":exon_seq })
  return exons

def intron_finding(record,exons):
  """Derive introns dynamically based on exon boundaries."""
  intron=[]
  for i in range(len(exons)-1) :
      start=exons[i]["end"]
      end=exons[i+1]["start"]
      intron_seq=record.seq[start:end]
      intron.append({"intron_number": i + 1,"start": start , "end":end ,"intron":intron_seq})
  return intron
# --- 3. Translation and Mutation Simulation ---
def translation_seq(seq):
  """Translate nucleotide sequence to protein starting from ATG up to the first stop codon."""
  start = seq.find("ATG")
  if start == -1:
    return "No Start Codon (ATG) found!"
        
  cds = seq[start:]
  cds_trimmed = cds[: len(cds) - (len(cds) % 3)]
  pro = Seq(str(cds_trimmed)).translate(to_stop=True)
  return pro


def mutate_sequence(sequence, position, new_nucleotide):
  """Introduce a single nucleotide point mutation at a specific 1-based position."""
  seq_list = list(sequence)
  seq_list[position - 1] = new_nucleotide
  return "".join(seq_list)
# --- 4. Main Execution Pipeline ---
def main():
  # Load HBB genomic record (U01317.1)
  
  record=Entrez_gb("U01317.1") 
  exons=exons_extraxt(record)
  intron=intron_finding(record,exons)

  # Simulate IVS-II-654 (C>T) pathogenic mutation in Intron 2
  mutate_seq=mutate_sequence(intron[1]["intron"],654,"T")
  
  # Extract the 73 bp aberrant cryptic exon as described in ZJU study
  E73 = mutate_seq[579:652]
  print(f"Length of cryptic exon: {len(E73)} bp")

  
  # Assemble Wild-Type and Mutant mRNA transcripts
  E1=exons[0]["exon"]
  E2=exons[1]["exon"]
  E3=exons[2]["exon"]

  #splicing
  spliced_seq_mu=E1+E2+E73+E3
  spliced_seq_wd=E1+E2+E3


  # Translate transcripts into protein sequences
  pro_wd=translation_seq(spliced_seq_wd)
  proIVS=translation_seq(spliced_seq_mu)

  
  # --- 5. Formatted Output Results ---
  print("\n" + "=" * 70)
  print("       BIOINFORMATICS ANALYSIS REPORT: HBB IVS-II-654 SPLICING")
  print("=" * 70)
  print(f" [+] Wild-Type Protein Length : {len(pro_wd)} amino acids")
  print(f" [+] Mutant Protein Length    : {len(proIVS)} amino acids (Truncated via PTC)")
  print("-" * 70)
  print(f"\n[Wild-Type Protein Sequence]\n{pro_wd}\n")
  print(f"[Mutant Protein Sequence (Pathogenic)]\n{proIVS}")
  print("=" * 70)


if __name__ == "__main__":
    main()

