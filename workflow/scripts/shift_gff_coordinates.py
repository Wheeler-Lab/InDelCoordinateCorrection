from bisect import bisect_right
from collections.abc import Sequence
import pandas as pd
from typing import List, Union, Dict

import geffa

snakemake: object

df = pd.read_table(snakemake.input.changes, sep=' ', header=None, names=['old_coordinate', 'new_coordinate', 'old', 'new'])
split = df.old_coordinate.str.split(':', expand=True)
df.loc[:, 'contig'] = split.iloc[:, 0]
df.loc[:, 'old_coordinate'] = split.iloc[:, 1]
split = df.new_coordinate.str.split(':', expand=True)
df.loc[:, 'new_coordinate'] = split.iloc[:, 1]
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
df = pd.concat([df, df.replace('.', '').apply(categorize_change, axis=1)], axis=1)
df.loc[:, 'old_coordinate_int'] = df.old_coordinate.map(lambda a: int(a.split('-')[0]))

class PieceWiseConstantShift:
    def __init__(self, coordinates: Sequence[int], shifts: Sequence[int]) -> None:
        self.coordinates: list[int] = [0]
        self.shifts: list[int] = [0]
        self.cumshifts: list[int] = [0]
        for c, s in zip(coordinates, shifts):
            if s != -100:
                self.coordinates.append(c)
                self.shifts.append(s)
                self.cumshifts.append(self.cumshifts[-1] + s)

    def __call__(self, x: Union[int,List[int]]) -> Union[int,List[int]]:
        if isinstance(x, Sequence):
            return [self(y) for y in x]
        index = bisect_right(self.coordinates, x)-1
        shift = self.shifts[index]
        cumshift = self.cumshifts[index]
        coordinate = self.coordinates[index]
        if x < coordinate - shift:
            new = coordinate
        else:
            new = x
        new += cumshift
        try:
            if new > self.coordinates[index+1]:
                print("WARNING: Shifting across shiftpoints!")
                new = self.coordinates[index+1]
        except IndexError:
            pass
        return new

translations:Dict[str, PieceWiseConstantShift] = {
    contig: PieceWiseConstantShift(g.old_coordinate_int, g.difference)
    for contig, g in df.groupby('contig')
}

# Validate - mark incomplete CDS as pseudogenes!!!

gff = geffa.GffFile(snakemake.input.gff, fasta_file=snakemake.input.fasta)
for seqreg in gff.sequence_regions.values():
    try:
        trans = translations[seqreg.name]
    except KeyError:
        # If we have no reads in a contig, it won't show up in the translation table
        continue
    marked_for_deletion = []
    for feature in seqreg.node_registry.values():
        start, end = feature.start, feature.end
        feature.start, feature.end = trans((feature.start, feature.end))
        try:
            feature.validate()
        except Exception as e:
            print(feature)
            print(feature.sequence_region.sequence[feature.start-5:feature.end+5])
            print(f"{start}->{feature.start}, {end}->{feature.end}")
            msg = str(e)
            if ("SLAS feature needs to be of length 2" in msg) or ("__len__() should return >= 0" in msg):
                print("Deleting this feature because it's been destroyed by an indel.")
                marked_for_deletion.append(feature)
            else:
                print(list(zip(trans.coordinates, trans.shifts, trans.cumshifts)))
                raise
    for feature in marked_for_deletion:
        feature.delete()
gff.save(snakemake.output.gff, include_sequences=True)