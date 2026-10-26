UML Диаграмма отношения классов в программе
===========================================

.. mermaid:: 

   classDiagram
    	direction TB

    class Record {
        +id: str
        +__init__(id: str)
        +__repr__() str
    }

    class SequenceRecord {
        +sequence: str
        +quality: List~int~
        +__init__(id: str, sequence: str, quality: List~int~)
        +__len__() int
    }

    class AlignmentRecord {
        +chrom: str
        +start: int
        +cigar: str
        +mapq: int
        +flag: int
        +end: int
        +__init__(id: str, chrom: str, start: int, cigar: str, mapq: int)
    }

    class VariantRecord {
        +chrom: str
        +pos: int
        +var_id: str
        +ref: str
        +alt: List~str~
        +qual: float
        +filter: str
        +info: Dict
        +start: int
        +end: int
        +__init__(chrom, pos, var_id, ref, alt, qual, filt, info)
        +intersects(chrom: str, start: int, end; int) bool 
        +__repr__() str
    }

    class Reader {
        <<abstract>>
        +filepath: Path
        +file: TextIO
        +__init__(filepath: str|Path)
        +read()* Iterator~Record~
        +__enter__() Reader
        +__exit__(exc_type, exc_val, exc_tb) void
        +close() void
    }

    class SequenceReader {
        <<abstract>>
        +__init__(filepath: str|Path)
        +read()* Iterator~SequenceRecord~
        +get_sequence(id: str, seq: str)* SequenceRecord
        +validate_sequence(seq: str)* bool
        +get_seq_count()* int
        +get_mean_seq_length()* float
    }

    class GenomicReader {
        <<abstract>>
        +_header_parsed: bool
        +_chromosomes: List~str~
        +__init__(filepath: str|Path)
        +__enter__() GenomicReader
        +_parse_header()*
        +read()* Iterator~Record~
        +get_chromosomes()* List~str~
        +validate_coordinate(chrom: str, pos: int)* bool
        +get_records_in_region(chrom: str, start: int, end: int)* Iterator~Record~
        +filter_records(**filters)* Iterator~Record~
        +get_chromosome_length(chrom: str) Optional~int~
        +get_statistics() Dict
        +close() void
    }

    class FastaAnalyzer {
        -_seq_count: int
        -_total_length: int
        +__init__(filepath: str|Path)
        +read() Iterator~SequenceRecord~
        +get_sequence(id: str, seq: str) SequenceRecord
        +validate_sequence(seq: str) bool
        +get_seq_count() int
        +get_mean_seq_length() float
    }

    class FastqAnalyzer {
        -_seq_count: int
        -_total_length: int
        +__init__(filepath: str|Path)
        +read() Iterator~SequenceRecord~
        +get_sequence(id: str, seq: str) SequenceRecord
        +validate_sequence(seq: str) bool
        +get_seq_count() int
        +get_mean_seq_length() float
        +phred_to_quality(phred_char: str, offset: int) int
        +get_sequences_with_quality() Iterator~tuple~
        +per_base_sequence_quality(output_file: str) Path
        +per_base_sequence_content(output_file: str) Path
        +sequence_length_distribution(output_file: str) Path
    }

    class SamReader {
        +header: Dict
        +__init__(filepath: str|Path)
        +_parse_header()
        +_parse_header_line(line: str)
        +get_header() Dict
        +get_header_fields() Dict
        +read() Iterator~AlignmentRecord~
        +_calc_aligned_length(cigar: str) int$
        +get_chromosomes() List~str~
        +validate_coordinate(chrom: str, pos: int) bool
        +calculate_coverage(chrom: str) Dict
        +filter_aligmnents(flag: int) Iterator~AlignmentRecord~
        +count_alignments() int
        +stats_by_chromosome() DataFrame
        +filter_by_region(chrom: str, start: int, end: int) Iterator~AlignmentRecord~
        +get_records_in_region(chrom: str, start: int, end: int) List~AlignmentRecord~
        +filter_records(**filters) List~AlignmentRecord~
    }

    class VcfReader {
        +headers: List~str~
        +header_groups: Dict
        +column_names: List~str~
        +samples: List~str~
        -_variants_count: int
        +__init__(filepath: str|Path)
        +get_header() Dict
        +get_header_group(group_name: str) List~str~
        +_parse_header()
        +read() Iterator~VariantRecord~
        +_parse_variant_line(line: str) VariantRecord
        +get_records_in_region(chrom: str, start: int, end: int) Iterator~VariantRecord~
        +filter_records(**filters) Iterator~VariantRecord~
        +get_variants_count() int
        +validate_coordinate(chrom: str, pos: int) bool
        +get_statistics() Dict      
        +get_chromosomes() List~str~
        +get_region_stats() DataFrame
        +get_variant_type_stats() DataFrame
        +reset() void
        +close() void
    }

    Record <|-- SequenceRecord
    Record <|-- AlignmentRecord
    Record <|-- VariantRecord

    Reader <|-- SequenceReader
    Reader <|-- GenomicReader

    SequenceReader <|-- FastaAnalyzer
    SequenceReader <|-- FastqAnalyzer

    GenomicReader <|-- SamReader
    GenomicReader <|-- VcfReader

    FastaAnalyzer ..> SequenceRecord : creates
    FastqAnalyzer ..> SequenceRecord : creates
    SamReader ..> AlignmentRecord : creates
    VcfReader ..> VariantRecord : creates
