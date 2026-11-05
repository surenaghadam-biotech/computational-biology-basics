"""
project
"""
import numpy as np

def zeros(seq1:str,seq2:str) -> np.ndarray:
        
    matrix=np.zeros((len(seq1),len(seq2)),dtype=int)
    for idx_1,char1 in enumerate(seq1):
        for idx_2,char2 in enumerate(seq2):
            if char1==char2:
                matrix[idx_1,idx_2]=1
    return matrix
                

def dotplot(matrix:np.ndarray,seq1:str,seq2:str)-> None:
    header="   "+ " ".join(list(seq2))
    print(header)
    print("  +"+"- "*(len(seq2)) )
    for idx_1,char1 in enumerate(seq1):
        rang=[]
        for idx_2 in range(len(seq2)):
            if matrix[idx_1,idx_2]==1:
                rang.append("●")
            else:
                rang.append(".")
        print(F"{char1} |"+ " ".join(rang))

def calculate(matrix: np.ndarray, seq1: str, seq2: str) :
    min_len=min(len(seq1),len(seq2))
    cont=np.trace(matrix[:min_len,:min_len])
    Mismatch=min_len-cont
    identity=(cont/min_len)*100
    print(f"cont:{cont}")
    print(f"Mismatch : {Mismatch}")
    print(f"identity : {identity}")

seq1="ATGGTT"#H
seq2="ATGGCT"#M
matrix=zeros(seq1,seq2)
dotplot(matrix,seq1,seq2)
calculate(matrix,seq1,seq2)
