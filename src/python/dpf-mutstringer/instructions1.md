

## Todo later:
# Citations

Ratnasingham, Sujeevan, and Paul D N Hebert. “bold: The Barcode of Life Data System (http://www.barcodinglife.org).” Molecular ecology notes vol. 7,3 (2007): 355-364. doi:10.1111/j.1471-8286.2007.01678.x
IGNORE FROM HERE
- protein:

- a single reference protein sequence that we will define mutation strings against
You MUST remember that during protein translation STOP signals will result in a deletion of the entirety of the protein sequence downstream of the signal.  You should use the reading frame determinined by the protein sequence to determine if something is a STOP signal or not.

Protein Mutation strings can be defined with the format:

- mutation: {reference AA}-{reference AA position}-{query AA}
  - Example: Reference seq: YYYQY Query seq: YYYTY should give: Q-4-T
- insertion: i{reference AA position to append AFTER}-{inserted AA sequence}
  - Example: Reference seq: YYYQY Query seq: YYYQHITYY should give: i4-HITY
- deletion: d{reference AA position to start deleting}-{reference AA position to stop deleting (inclusive)}
  - Example: Reference seq: YYYQY Query seq: QY should give: d1-3
UNTIL HERE


To accomplish this, lets new build the code that will actually be simulating the data.
You should create code that will on the fly, and optionally write to a tsv database, generate a reference DNA sequence,
then use that to generate a query sequence along with an associated mutation string.

I would like you to complete the code for this file (named src/dpf_mutstringer/simulate_data.py)
