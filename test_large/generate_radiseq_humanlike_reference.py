#!/usr/bin/env python3
"""
Generate a compact human-like reference compatible with RadiSeq mapping type 1.

RadiSeq mapping type 1 expects the SDD chromosome order:

    1, 1, 2, 2, ..., 22, 22, X, Y

This represents 46 physical chromosomes, but the input reference FASTA contains
24 sequence records:

    chr1, chr2, ..., chr22, chrX, chrY

RadiSeq duplicates each autosomal reference sequence internally and leaves X
and Y as single copies.

Outputs
-------
1. radiseq_humanlike_reference.fa
2. radiseq_humanlike_reference_sdd_mapping.tsv
3. radiseq_humanlike_reference_sdd_header.txt
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path
from typing import BinaryIO


DNA_TRANSLATION_TABLE = bytes(
    ord("ACGT"[value & 0b11]) for value in range(256)
)


def linearly_decreasing_lengths(
    count: int,
    maximum_length: int,
    minimum_length: int,
) -> list[int]:
    """Create strictly decreasing, approximately linearly spaced lengths."""
    if count < 2:
        raise ValueError("count must be at least 2")
    if maximum_length <= minimum_length:
        raise ValueError("maximum_length must exceed minimum_length")
    if minimum_length <= 0:
        raise ValueError("minimum_length must be positive")

    span = maximum_length - minimum_length
    lengths = [
        round(maximum_length - index * span / (count - 1))
        for index in range(count)
    ]

    lengths[0] = maximum_length
    lengths[-1] = minimum_length

    if any(left <= right for left, right in zip(lengths, lengths[1:])):
        raise ValueError(
            "The requested range is too narrow to produce strictly "
            "decreasing integer lengths."
        )

    return lengths


def random_dna(rng: random.Random, length: int) -> bytes:
    """Generate deterministic random DNA."""
    return rng.randbytes(length).translate(DNA_TRANSLATION_TABLE)


def write_wrapped(
    output: BinaryIO,
    sequence: bytes,
    line_width: int,
) -> None:
    for start in range(0, len(sequence), line_width):
        output.write(sequence[start : start + line_width])
        output.write(b"\n")


def write_contig(
    output: BinaryIO,
    name: str,
    length: int,
    rng: random.Random,
    line_width: int,
    chunk_size: int,
) -> None:
    output.write(f">{name} length={length}\n".encode("ascii"))

    remaining = length
    while remaining > 0:
        current_size = min(chunk_size, remaining)
        write_wrapped(output, random_dna(rng, current_size), line_width)
        remaining -= current_size


def write_reference_fasta(
    output_path: Path,
    autosome_lengths: list[int],
    x_length: int,
    y_length: int,
    seed: int,
    line_width: int,
) -> None:
    """Write chr1--chr22, chrX and chrY in that exact order."""
    rng = random.Random(seed)

    chunk_size = 1_000_000
    chunk_size -= chunk_size % line_width
    if chunk_size == 0:
        chunk_size = line_width

    with output_path.open("wb") as output:
        for chromosome, length in enumerate(autosome_lengths, start=1):
            write_contig(
                output=output,
                name=f"chr{chromosome}",
                length=length,
                rng=rng,
                line_width=line_width,
                chunk_size=chunk_size,
            )

        write_contig(
            output=output,
            name="chrX",
            length=x_length,
            rng=rng,
            line_width=line_width,
            chunk_size=chunk_size,
        )
        write_contig(
            output=output,
            name="chrY",
            length=y_length,
            rng=rng,
            line_width=line_width,
            chunk_size=chunk_size,
        )


def physical_chromosome_records(
    autosome_lengths: list[int],
    x_length: int,
    y_length: int,
) -> list[tuple[int, str, str, int]]:
    """
    Return:
        (SDD chromosome ID, biological chromosome, copy, length_bp)

    SDD IDs:
        1,2   -> chromosome 1 copies 1 and 2
        3,4   -> chromosome 2 copies 1 and 2
        ...
        43,44 -> chromosome 22 copies 1 and 2
        45    -> X
        46    -> Y
    """
    records: list[tuple[int, str, str, int]] = []
    sdd_id = 1

    for chromosome, length in enumerate(autosome_lengths, start=1):
        for copy in ("1", "2"):
            records.append((sdd_id, str(chromosome), copy, length))
            sdd_id += 1

    records.append((sdd_id, "X", "1", x_length))
    sdd_id += 1
    records.append((sdd_id, "Y", "1", y_length))

    return records


def write_mapping(
    output_path: Path,
    records: list[tuple[int, str, str, int]],
) -> None:
    with output_path.open("w", encoding="utf-8", newline="\n") as output:
        output.write(
            "sdd_chromosome_id\t"
            "biological_chromosome\t"
            "copy\t"
            "reference_contig\t"
            "length_bp\t"
            "length_mbp\n"
        )

        for sdd_id, chromosome, copy, length in records:
            output.write(
                f"{sdd_id}\t"
                f"{chromosome}\t"
                f"{copy}\t"
                f"chr{chromosome}\t"
                f"{length}\t"
                f"{length / 1_000_000:.6f}\n"
            )


def write_sdd_header(
    output_path: Path,
    records: list[tuple[int, str, str, int]],
) -> None:
    sizes = ", ".join(
        f"{length / 1_000_000:.6f}"
        for _, _, _, length in records
    )

    with output_path.open("w", encoding="utf-8", newline="\n") as output:
        output.write(f"Chromosome sizes, {len(records)}, {sizes};\n")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a compact human-like RadiSeq reference with "
            "22 autosome types plus X and Y."
        )
    )

    parser.add_argument(
        "--output-prefix",
        default="radiseq_humanlike_reference",
        help="Output prefix. Default: radiseq_humanlike_reference",
    )
    parser.add_argument(
        "--largest-autosome",
        type=int,
        default=2_000_000,
        help="Length of chr1 in bp. Default: 2000000",
    )
    parser.add_argument(
        "--smallest-autosome",
        type=int,
        default=600_000,
        help="Length of chr22 in bp. Default: 600000",
    )
    parser.add_argument(
        "--x-length",
        type=int,
        default=800_000,
        help="Length of chrX in bp. Default: 800000",
    )
    parser.add_argument(
        "--y-length",
        type=int,
        default=500_000,
        help="Length of chrY in bp. Default: 500000",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260707,
        help="Random seed. Default: 20260707",
    )
    parser.add_argument(
        "--line-width",
        type=int,
        default=80,
        help="FASTA line width. Default: 80",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    if args.x_length <= 0 or args.y_length <= 0:
        raise ValueError("X and Y lengths must be positive")
    if args.line_width <= 0:
        raise ValueError("line width must be positive")

    autosome_lengths = linearly_decreasing_lengths(
        count=22,
        maximum_length=args.largest_autosome,
        minimum_length=args.smallest_autosome,
    )

    prefix = Path(args.output_prefix)
    fasta_path = prefix.with_suffix(".fa")
    mapping_path = prefix.parent / f"{prefix.name}_sdd_mapping.tsv"
    header_path = prefix.parent / f"{prefix.name}_sdd_header.txt"

    for path in (fasta_path, mapping_path, header_path):
        path.parent.mkdir(parents=True, exist_ok=True)

    write_reference_fasta(
        output_path=fasta_path,
        autosome_lengths=autosome_lengths,
        x_length=args.x_length,
        y_length=args.y_length,
        seed=args.seed,
        line_width=args.line_width,
    )

    records = physical_chromosome_records(
        autosome_lengths=autosome_lengths,
        x_length=args.x_length,
        y_length=args.y_length,
    )

    write_mapping(mapping_path, records)
    write_sdd_header(header_path, records)

    haploid_reference_length = (
        sum(autosome_lengths) + args.x_length + args.y_length
    )
    diploid_cell_length = (
        2 * sum(autosome_lengths) + args.x_length + args.y_length
    )

    print(f"Created reference FASTA: {fasta_path}")
    print(f"Created SDD mapping:     {mapping_path}")
    print(f"Created SDD header:      {header_path}")
    print("Reference records:       24 (chr1-chr22, chrX, chrY)")
    print("Physical SDD chromosomes:46")
    print(f"Reference length:        {haploid_reference_length:,} bp")
    print(f"Diploid cell DNA length: {diploid_cell_length:,} bp")
    print(f"Seed:                    {args.seed}")


if __name__ == "__main__":
    main()
