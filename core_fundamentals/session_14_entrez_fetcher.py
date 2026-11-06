"""
Session 14: Advanced Entrez Workflow.
Parsing and displaying CDS features with locations.

"""
from typing import Optional
from urllib.error import URLError
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from utils import get_ncbi_email

def fetch_genbank_record(accession_id)->Optional[SeqRecord] :

    Entrez.email=get_ncbi_email()
    handle=None
    try:
        handle=Entrez.efetch(db="nuccore",id=accession_id,rettype="gb",retmode="text")
        record=SeqIO.read(handle,"genbank")
        return record
    except URLError as err:
        print(f"Network issue: {err}")
        return None
    except ValueError as parse_err:
        print(f"[Parsing Error] Could not parse GenBank record: {parse_err}")
        return None
    finally:
        if handle:
            handle.close()

if __name__ == "__main__":
        
    record=fetch_genbank_record("NG_000007")
    if record:
        cds_count=0
        cds_found =False
        for i in record.features:
            if i.type=="CDS":
                cds_count+=1
                cds_found =True
                seq_i=i.extract(record.seq)
                print(f"\n CDS #{cds_count} | Location: {i.location} | Length: {len(seq_i)} bp")
        if not cds_found :
            print(" CDS not find")
        else:
            print(f"\n[+] Total CDS found: {cds_count} \n")
