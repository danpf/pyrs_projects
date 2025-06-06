This will be a pytorch based machine learning algorithm to build 'mutation' strings for a protein sequence, given a list of possible reference sequences.

The inputs to this pipeline will be:

- one or more DNA sequences that we will define mutation strings for (we can also call them query sequences)
- a single reference dna sequence that we will define mutation strings against
- **Note** DNA sequences sometimes have extra data at the 5' and/or 3' ends, so the protein sequence is not necessarily a direct translation of the DNA sequence.

The output of this pipeline will be:

- a DNA mutstring

DNA Mutation strings can be defined with the format:

- mutation: {reference NCP}-{reference NCP position}-{query NCP}
  - Example: Reference seq: ACGTACGTACGT Query seq: ACGGACGTACGT should give: T-4-G
- insertion: i{reference NCP position to append AFTER}-{inserted NCP sequence}
  - Example: Reference seq: ACTT Query seq: ACTACCCCT should give: i3-ACCCC
- deletion: d{reference NCP position to start deleting}-{reference NCP position to stop deleting (inclusive)}
  - Example: Reference seq: AACCGGTT Query seq: AATT should give: d3-6

