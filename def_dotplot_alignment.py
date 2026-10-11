"""
Session 08: Pairwise Sequence Alignment using Dot Plot Matrix
Author: Surena Ghadam
Description: Visualizes sequence similarities, mutations, and indels using ASCII Dot Plot.
"""
def generate_dot_plot(seq1, seq2):
    print(" ",end=" ")
    for s2 in seq2 :
        print(s2,end=" ")
    print()

    for s1 in seq1:
        print(s1,end=" ")
        for s2 in seq2:
            if s2==s1 :
                print("*",end=" ")
            else: 
                print(" ",end=" ")
        print()


