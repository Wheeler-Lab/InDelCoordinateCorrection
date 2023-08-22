from bisect import bisect_right
from collections.abc import Sequence
import pandas as pd
from typing import List, Union, Dict

import geffa

# Give the snakemake object a type...
snakemake: object

# Read the pilon changes file
df = pd.read_table(snakemake.input.changes, sep=' ', header=None, names=['old_coordinate', 'new_coordinate', 'old', 'new'])
# Split the coordinates into contig name and position
split = df.old_coordinate.str.split(':', expand=True)
df.loc[:, 'contig'] = split.iloc[:, 0]
df.loc[:, 'old_coordinate'] = split.iloc[:, 1]
split = df.new_coordinate.str.split(':', expand=True)
df.loc[:, 'new_coordinate'] = split.iloc[:, 1]

# Function to categorise the each row (deletion, insertion or snp)
# Also calculates the change in length (negative for deletion, positive for insertion, zero for snps)
def categorize_change(row):
    lold = len(row.old)
    lnew = len(row.new)
    diff = lnew - lold
    if diff == 0:
        if lold == 1:
            category = 'snp'
        elif lold == 0:
            raise ValueError('No change!?!')
        else:
            category = 'indel'
    elif diff < 0:
        category = 'deletion'
    else:
        category = 'insertion'
    return pd.Series({'difference': diff, 'category': category})

# Incorporate the categorisation into the table
df = pd.concat([df, df.replace('.', '').apply(categorize_change, axis=1)], axis=1)
df.loc[:, 'old_coordinate_int'] = df.old_coordinate.map(lambda a: int(a.split('-')[0]))

# Class that takes a list of change coordinates and associated shifts and allows us to 
# then transform old feature coordinates into new ones.
class PieceWiseConstantShift:
    def __init__(self, coordinates: Sequence[int], shifts: Sequence[int]) -> None:
        self.coordinates: list[int] = [0]
        self.shifts: list[int] = [0]
        self.cumshifts: list[int] = [0]
        
        # Calculate cumulative shifts
        for c, s in zip(coordinates, shifts):
            if s != 0: # Filter out snps, they don't shift anything
                self.coordinates.append(c)
                self.shifts.append(s)
                self.cumshifts.append(self.cumshifts[-1] + s)

    # This takes a coordinate and shifts it according to the changes given above
    def __call__(self, x: Union[int,List[int]]) -> Union[int,List[int]]:
        if isinstance(x, Sequence):
            return [self(y) for y in x]
        # Find the index of the shift to the left of the given coordinate
        index = bisect_right(self.coordinates, x)-1
        shift = self.shifts[index]
        cumshift = self.cumshifts[index]
        coordinate = self.coordinates[index]
        # If the old coordinate is within a deletion or insertion region, we need to shift it beyond
        if x < coordinate - shift:
            new = coordinate
        else:
            new = x
        new += cumshift
        try:
            # If we shifted beyond the next change, we need to correct for that
            if new > self.coordinates[index+1]:
                new = self.coordinates[index+1]
        except IndexError:
            # We were at the last entry, so nothing to see here.
            pass
        return new

# Build a translator for each contig
translations:Dict[str, PieceWiseConstantShift] = {
    contig: PieceWiseConstantShift(g.old_coordinate_int, g.difference)
    for contig, g in df.groupby('contig')
}

# TODO: Validate - mark incomplete CDS as pseudogenes!!!

# Load the parental GFF file and the new sequence fasta file
gff = geffa.GffFile(snakemake.input.gff, fasta_file=snakemake.input.fasta, ignore_unknown_feature_types=True)
for seqreg in gff.sequence_regions.values():
    # Iterate over contigs
    try:
        # Load the correct translator
        trans = translations[seqreg.name]
    except KeyError:
        # If we have no reads in a contig, it won't show up in the translation table
        continue

    marked_for_deletion = []
    # Iterate over all features in the contig
    for feature in seqreg.node_registry.values():
        # Translate the start and end coordinates
        start, end = feature.start, feature.end
        feature.start, feature.end = trans((feature.start, feature.end))
        try:
            # Validate the feature with the new coordinates
            feature.validate()
        except Exception as e:
            # The new feature isn't correct!
            # Print some debug info
            print(feature)
            print(feature.sequence_region.sequence[feature.start-5:feature.end+5])
            print(f"{start}->{feature.start}, {end}->{feature.end}")
            msg = str(e)
            # Delete the feature if it's a feature that has been corrupted by an indel
            if ("SLAS feature needs to be of length 2" in msg) or ("__len__() should return >= 0" in msg):
                print("Deleting this feature because it's been destroyed by an indel.")
                marked_for_deletion.append(feature)
            elif "Exon parent needs to be an RNA" in msg:
                print("Invalid exon found - skipping because it's not that severe.")
            else:
                # Fail otherwise, don't know what to do.
                print(list(zip(trans.coordinates, trans.shifts, trans.cumshifts)))
                raise
    # Delete the features that were marked for deletion
    # (couldn't do that earlier, can't change a container while iterating over it)
    for feature in marked_for_deletion:
        feature.delete()
# Save the new annotation GFF
gff.save(snakemake.output.gff, include_sequences=True)
