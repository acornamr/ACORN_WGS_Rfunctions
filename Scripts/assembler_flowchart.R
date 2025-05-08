# This script is used to create flowchart diagram for the poster about the automated assembling pipeline 
# of the container.
# Written and maintained by Tung Trinh
# March 13th 2025
# Version 1.0
# For more informaiton, please contact to tungts@oucru.org 
####################################################################################################
##### required packages and built-in functions #####
library(DiagrammeR)

### create flow chart of data process
jpeg(filename = "ESCMID2025/assembler_flowchart.jpeg", res = 300, unit = "in", height = 4, width = 4.5)
grViz("
digraph {

  # graph attributes
  graph [overlap = true]

  # node attributes
  node [shape = box,
        fontname = Arial,
        color = '#00ADAD',
        fill = 'white'
        size = 3
        #style = filled
        ]

  # edge attributes
  edge [color = gray]

  # node statements
  A [label = 'Input raw Reads\n(.fastq.gz format)']
  B [label = 'Quality control using FastQC']
  C1 [label = 'Trimming reads using Trimmomatic']
  C2 [label = 'Estimating genome size and read length\nusing kmc']
  D1 [label = 'Contamination Check with confindr']
  D2 [label = 'Filtering the read by length with Nanoq\n(Min length 1000)']
  E1 [label = 'Count the number of reads and estimate\ngenome size using Mash']
  E2 [label = 'Reduce FASTQ files to a sensible depth']
  F1 [label = 'Downsample reads']
  F2 [label = 'Assemble with Flye']
  G1 [label = 'Merge reads using Flash where the insert size is small']
  G2 [label = 'Polish assembly with MEdaka']
  H1 [label = 'Assemble reads using SPAdes']
  H2 [label = 'Remove contigs that are too short, too low coverage,\n or pure homopolymers']
  I1 [label = 'Assess species identification using bactinspector']
  I2 [label = 'Reorient contigs from final FASTA using dnaapler']
  J1 [label = 'Assess assembly quality using Quast']
  J2 [label = 'Assess assembly quality using Quast']
  K2 [label = 'Summarise QC metrics with QuailFyr']
  K1 [label = 'Summarise QC metrics with QuailFyr']
  
  
  # edge statements
  A->B 
  B->C1 [label = 'Short read input' ]
  B->C2 [label = 'Long read input']
  C1->D1  
  C2->D2
  D1->E1
  D2->E2
  E1->F1
  E2->F2
  F1->G1
  F2->G2
  G1->H1
  G2->H2
  H1->I1 
  H2->I2 
  I1->J1
  I2->J2
  J1->K1
  J2->K2
}")
dev.off()
