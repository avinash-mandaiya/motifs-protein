## Motifs for secondary structure

The PCA analysis of motifs capturing the interaction between two peptide bond planes (C<sub>alpha</sub> - NH - CO - C<sub>alpha</sub>) shows the distinct secondary structure.  

![BB](https://github.com/user-attachments/assets/afd11d32-a73e-473f-8824-a23ef91a5e84)

These motifs can help form secondary structures in a constraint satisfaction formulation for protein folding.

## Rotamer Motifs

rotamer.py generates rotamer motifs for all the amino acids from a list of ~1000 high-resolution proteins. Rotamer motifs for each amino acid consist of all the non-hydrogen atoms except "O" of the carboxyl group. These rotamer motifs are then visualized using PCA (principal component analysis) to understand the degrees of freedom and compare with the theoretical understanding based on the hybridization geometry of each atom in the rotamer. 

Points that form a circular path represent a single rotational degree of freedom. Due to the steric hindrance, the density of points is larger in some areas of the circular region.

![rot_PCA](https://github.com/user-attachments/assets/800d6bee-93e5-40d7-b344-cf1d4a0b9565)


