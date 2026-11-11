"""Модели записей биологических данных.

Содержит базовый класс Record и специализированные записи
для FASTA/FASTQ, SAM и VCF.

Examples:
    >>> from record import SequenceRecord
    >>> rec = SequenceRecord("seq1", "ATCG")
    >>> len(rec)
    4
"""
from typing import List, Dict, Any, Optional

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

    Хранит полную информацию о мутации: хромосома, позиция, аллели,
    качество, фильтр и разобранный словарь INFO. Поддерживает
    интервальные проверки (`start`, `end`, `intersects`), эквивалентные
    поведению bedtools intersect.

    Attributes:
        chrom: Название хромосомы (нормализованное, stripped).
        pos: Позиция варианта (1-based).
        var_id: Поле ID из VCF ('.', если отсутствует).
        ref: Референсный аллель (в верхнем регистре).
        alt: Список альтернативных аллелей (в верхнем регистре);
            пустой список для записи `ALT="."`.
        qual: Качество варианта (float) или None для `"."`.
        filter: Поле FILTER.
        info: Словарь аннотаций из поля INFO.
        start: Свойство — начало мутации (1-based), совпадает с `pos`.
        end: Свойство — конец мутации по длине REF.

    Examples:
        >>> rec = VariantRecord(
        ...     chrom="chr1", pos=100, var_id=".",
        ...     ref="A", alt=["T"], qual=50.0, filt="PASS",
        ...     info={"DP": "30"},
        ... )
        >>> rec.id
        'chr1:100'
        >>> rec.start, rec.end
        (100, 100)
        >>> rec.intersects("chr1", 90, 110)
        True
        >>> VariantRecord(
        ...     chrom="chr1", pos=100, var_id=".",
        ...     ref="AT", alt=["A"], qual=None, filt="PASS", info={},
        ... ).end
        101
    """

    def __init__(
        self,
        chrom: str,
        pos: int,
        var_id: str,
        ref: str,
        alt: List[str],
        qual: Optional[float],
        filt: str,
        info: Dict[str, Any],
    ):
        """Инициализирует запись варианта из уже разобранных полей.

        Args:
            chrom: Название хромосомы.
            pos: Позиция варианта (1-based), уже приведённая к int.
            var_id: Поле ID из VCF.
            ref: Референсный аллель.
            alt: Список альтернативных аллелей.
            qual: Качество варианта (float) или None.
            filt: Поле FILTER.
            info: Разобранный словарь INFO.

        Examples:
            >>> VariantRecord("chr1", 100, ".", "A", ["T", "G"], 50.0,
            ...               "PASS", {"DP": "30"}).alt
            ['T', 'G']
        """
        super().__init__(f"{chrom}:{pos}")
        self.chrom = chrom
        self.pos = pos
        self.var_id = var_id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filt
        self.info = info

    @property
    def start(self) -> int:
        """Начало мутации (1-based).

        Returns:
            int: Позиция начала (совпадает с `pos`).
        """
        return self.pos

    @property
    def end(self) -> int:
        """Конец мутации.

        Вычисляется по длине референсного аллеля REF
        (SNV: len=1, INDEL: len > 1).

        Returns:
            int: Позиция конца (1-based, включительно).
        """
        return self.pos + max(len(self.ref), 1) - 1

    def intersects(self, chrom: str, start: int, end: int) -> bool:
        """Проверка пересечения с интервалом [start, end].

        Реализует логику bedtools intersect: отрезки [A, B] и [C, D]
        пересекаются, если не выполняется (B < C или A > D).

        Args:
            chrom: Хромосома запрашиваемого региона.
            start: Начало региона (1-based, включительно).
            end: Конец региона (1-based, включительно).

        Returns:
            bool: True, если вариант пересекается с регионом.

        Examples:
            >>> VariantRecord("chr1", 100, ".", "A", "T", ".", "PASS", ".")\
                    .intersects("chr1", 90, 110)
            True
            >>> VariantRecord("chr1", 100, ".", "A", "T", ".", "PASS", ".")\
                    .intersects("chr2", 90, 110)
            False
        """
        if self.chrom != str(chrom):
            return False
        return not (self.end < start or self.start > end)

    def __repr__(self) -> str:
        """Строковое представление записи варианта.

        Returns:
            str: Представление вида `<VariantRecord chrom:pos REF>ALT>`.

        Examples:
            >>> repr(VariantRecord("chr1", 100, ".", "A", "T,G", ".", "PASS", "."))
            '<VariantRecord chr1:100 A>T,G>'
        """
        alt_repr = ",".join(self.alt) if self.alt else "."
        return f"<VariantRecord {self.chrom}:{self.pos} {self.ref}>{alt_repr}>"