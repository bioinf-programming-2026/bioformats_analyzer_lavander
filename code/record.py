"""Модели записей биологических данных.

Содержит базовый класс Record и специализированные записи
для FASTA/FASTQ, SAM и VCF.

Examples:
    >>> from record import SequenceRecord
    >>> rec = SequenceRecord("seq1", "ATCG")
    >>> len(rec)
    4
"""

class Record:
    """Базовая запись биологических данных.

    Attributes:
        id: Идентификатор записи.

    Examples:
        >>> Record("rec1")
        <Record id=rec1>
    """


    def __init__(self, id: str):
        """Инициализирует базовую запись.

        Args:
            id: Идентификатор записи.

        Examples:
            >>> r = Record("x")
            >>> r.id
            'x'
        """
        self.id = id

    def __repr__(self):
        """Возвращает краткое строковое представление записи.

        Returns:
            str: Представление вида `<Record id=...>`.

        Examples:
            >>> repr(Record("abc"))
            '<Record id=abc>'
        """
        return f"<Record id={self.id}>"


class SequenceRecord(Record):
    """Запись биологической последовательности для FASTA/FASTQ.

    Attributes:
        sequence: Нуклеотидная последовательность.
        quality: Список числовых значений качества или None для FASTA.

    Examples:
        >>> rec = SequenceRecord("read1", "ATCG", quality=[40, 40, 40, 40])
        >>> rec.sequence, len(rec.quality)
        ('ATCG', 4)
    """

    def __init__(self, id: str, sequence: str, quality: list[int] | None = None):
        """Инициализирует запись последовательности.

        Args:
            id: Идентификатор последовательности.
            sequence: Нуклеотидная последовательность.
            quality: Числовые значения Phred-качества для FASTQ.

        Examples:
            >>> SequenceRecord("s1", "AT").quality is None
            True
        """
        super().__init__(id)
        self.sequence = sequence
        self.quality = quality


    def __len__(self):
        """Возвращает длину последовательности.

        Returns:
            int: Количество нуклеотидов.

        Examples:
            >>> len(SequenceRecord("s1", "ATCG"))
            4
        """
        return len(self.sequence)


class AlignmentRecord(Record):
    """Запись выравнивания из SAM-файла.

    Attributes:
        chrom: Хромосома, к которой выполнено выравнивание.
        start: Начальная позиция выравнивания.
        cigar: CIGAR-строка выравнивания.
        mapq: Качество маппинга.
        flag: SAM-флаг; может добавляться ридером динамически.
        end: Конечная позиция выравнивания; может добавляться ридером.

    Examples:
        >>> rec = AlignmentRecord("r1", "chr1", 100, "50M", 60)
        >>> rec.chrom, rec.start
        ('chr1', 100)
    """

    def __init__(self, id: str, chrom: str, start: int, cigar: str, mapq: int):
        """Инициализирует запись выравнивания.

        Args:
            id: Идентификатор прочтения.
            chrom: Название хромосомы.
            start: Начальная позиция выравнивания.
            cigar: CIGAR-строка.
            mapq: Качество маппинга.

        Examples:
            >>> AlignmentRecord("r1", "chr1", 100, "50M", 60).mapq
            60
        """
        super().__init__(id)
        self.chrom = chrom
        self.start = start
        self.cigar = cigar
        self.mapq = mapq


class VariantRecord(Record):
    """Запись геномного варианта из VCF-файла.

    Attributes:
        chrom: Название хромосомы.
        pos: Позиция варианта.
        ref: Референсный аллель.
        alt: Альтернативный аллель.
        info: Словарь дополнительных аннотаций из поля INFO.

    Examples:
        >>> rec = VariantRecord("chr1", 100, "A", "T", {"DP": "30"})
        >>> rec.id
        'chr1:100'
    """

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, info: dict):
        """Инициализирует запись варианта.

        Args:
            chrom: Название хромосомы.
            pos: Позиция варианта.
            ref: Референсный аллель.
            alt: Альтернативный аллель.
            info: Разобранное поле INFO.

        Examples:
            >>> VariantRecord("chr1", 100, "A", "T", {}).pos
            100
        """
        super().__init__(f"{chrom}:{pos}")
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info