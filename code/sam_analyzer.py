import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, List, Any
import re
from abstract import GenomicReader
from record import AlignmentRecord

_SAM_ID_KEYS = {
    "@SQ": "SN",
    "@RG": "ID",
    "@PG": "ID",
}

"""Ридер SAM-файлов.

Разбирает заголовок и выравнивания, предоставляет статистику,
фильтрацию, покрытие и выборки по регионам.

Examples:
    >>> with SamReader("aln.sam") as r:
    ...     print(r.count_alignments())
    ...     for rec in r.filter_records(chrom="chr1", min_mapq=30):
    ...         print(rec.id, rec.start)
"""

class SamReader(GenomicReader):
    """
    Парсер SAM-файлов для работы с данными геномных выравниваний.
    
    Наследуется от GenomicReader и предоставляет методы для чтения,
    анализа и фильтрации данных выравнивания из SAM-файлов.
    
    Attributes:
        filepath (Path): Путь к SAM-файлу
        file (IO): Файловый объект для чтения данных
        header (Dict[str, List[str]]): Заголовок SAM-файла, разобранный по тегам

    Examples:
        >>> with SamReader("aln.sam") as r:
        ...     chroms = r.get_chromosomes()
        ...     stats = r.stats_by_chromosome()
        ...     print(stats.head())
    """
    
    def __init__(self, filepath: str | Path):
        """Инициализирует SAM-ридер.

        Args:
            filepath: Путь к SAM-файлу.

        Examples:
            >>> SamReader("aln.sam").header
            {}
        """
        super().__init__(filepath)
        self.header: Dict[str, List[str]] = {}

    def _parse_header(self):
        """
        Приватный метод для парсинга заголовка SAM-файла.
        Разбирает заголовок SAM-файла и оставляет позицию чтения на данных.

        Читает все строки, начинающиеся с '@', и сохраняет их в атрибуте header.
        После парсинга возвращает позицию чтения файла в начало.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     r._parse_header()
            ...     '@SQ' in r.header
            True
        """
        if not self.file:
            self.file = open(self.filepath, "r", encoding="utf-8")
        self.file.seek(0)
        while True:
            pos = self.file.tell()
            line = self.file.readline()
            if not line.startswith("@"):
                self.file.seek(pos)
                break
            self._parse_header_line(line)
        self._header_parsed = True

    def _parse_header_line(self, line: str):
        """Парсит одну строку заголовка SAM-файла.

        Args:
            line: Строка заголовка, начинающаяся с `@`.

        Examples:
            >>> r = SamReader("aln.sam")
            >>> r.header = {}
            >>> r._parse_header_line("@SQ\\tSN:chr1\\tLN:100")
            >>> r.header["@SQ"]
            ['SN:chr1\\tLN:100']
        """
        parts = line.strip().split("\t")
        tag = parts[0]
        self.header.setdefault(tag, []).append("\t".join(parts[1:]))

    def get_header(self) -> Dict[str, List[str]]:
        """
        Возвращает заголовок SAM-файла.

        Returns:
            Dict[str, List[str]]: связывает теги заголовка (например, '@SQ', '@RG')
            и значения - списками соответствующих данных

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     list(r.get_header().keys())[:2]
            ['@HD', '@SQ']
        """
        return self.header
    
    def get_header_fields(self) -> Dict[str, Dict[str, Dict[str, str]]]:
        """Возвращает полностью разобранные поля заголовка.

        Returns:
            Dict[str, Dict[str, Dict[str, str]]]: Вложенный словарь
            `{тег: {ID: {KEY: VALUE}}}` для @SQ, @RG и @PG.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     sq = r.get_header_fields().get("@SQ", {})
            ...     sq["chr1"]["LN"]
            '248956422'
        """
        result: Dict[str, Dict[str, Dict[str, str]]] = {}
        for tag, lines in self.header.items():
            id_key = _SAM_ID_KEYS.get(tag)
            if id_key is None:
                continue
            tag_result: Dict[str, Dict[str, str]] = {}
            for line in lines:
                fields: Dict[str, str] = {}
                for field in line.split("\t"):
                    if ":" in field:
                        k, v = field.split(":", 1)
                        fields[k] = v
                if id_key in fields:
                    tag_result[fields[id_key]] = fields
            if tag_result:
                result[tag] = tag_result
        return result

    def read(self) -> Iterator[AlignmentRecord]:
        """
        Основной метод для чтения записей выравнивания из SAM-файла.
        
        Пропускает строки заголовка и записи с некорректными данными.
        Генерирует объекты AlignmentRecord для каждой валидной записи.
        
        Yields:
            AlignmentRecord: Объект с данными выравнивания
        
        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     rec = next(r.read())
            ...     rec.chrom, rec.start
            ('chr1', 100)

        Notes:
            - Пропускает записи где RNAME или CIGAR равны '*'
            - Автоматически вычисляет конечную позицию выравнивания
        """
        if not self.file:
            self.file = open(self.filepath, "r", encoding="utf-8")
            self.file.seek(0)
        if not self._header_parsed:
            self._parse_header()

        for line in self.file:
            if line.startswith("@"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue

            if fields[2] == "*" or fields[5] == "*":
                continue
            
            try:
                flag = int(fields[1])
                pos = int(fields[3])
                mapq = int(fields[4])
            
            except ValueError:
                print(f"Ошибка: Пропущена битая SAM-строка: {line[:80]!r}")
                continue

            qname = fields[0] 
            rname =fields[2]
            cigar = fields[5]

            aligned_len = self._calc_aligned_length(cigar)
            end_pos = pos + aligned_len - 1 if aligned_len > 0 else int(pos)

            rec = AlignmentRecord(
                id=qname, chrom=rname, start=pos, cigar=cigar, mapq=mapq
            )
            rec.flag = flag
            rec.end = end_pos
            yield rec

    @staticmethod
    def _calc_aligned_length(cigar: str) -> int:
        """
        Вычисляет длину выравнивания на основе CIGAR-строки.

        Учитываются операции M (match), D (deletion), N (skip), = (match), X (mismatch).

        Args:
            cigar: CIGAR-строка в стандартном формате

        Returns:
            int: Суммарная длина выравнивания.

        Examples:
            >>> SamReader._calc_aligned_length("50M10D20N")
            80
            >>> SamReader._calc_aligned_length("*")
            0
        """
        if not cigar or cigar == "*":
            return 0

        total = 0
        for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
            if op in ("M", "D", "N", "=", "X"):
                total += int(length)
        return total

    def get_chromosomes(self) -> List[str]:
        """Возвращает хромосомы, встречающиеся в выравниваниях.

        Returns:
            List[str]: Отсортированный список названий хромосом.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     "chr1" in r.get_chromosomes()
            True
        """
        chromosomes = set()
        # Сохраняем текущую позицию в файле
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            if rec.chrom != "*" and rec.chrom:  # проверяем что не пустая и не *
                chromosomes.add(rec.chrom)

        # Восстанавливаем позицию в файле
        if self.file:
            self.file.seek(current_pos)
        return sorted(chromosomes)

    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """
        Проверяет корректность геномной координаты
        Args:
            chrom: Название хромосомы
            pos: Позиция на хромосоме
            
        Returns:
            True если хромосома есть в файле и позиция положительна
        
        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     r.validate_coordinate("chr1", 100)
            True
        """
        return chrom in self.get_chromosomes() and pos > 0

    def calculate_coverage(self, chrom: str) -> Dict[int, int]:
        """Вычисляет покрытие для указанной хромосомы.

        Args:
            chrom: Название хромосомы для анализап

        Returns:
            Dict[int, int]: Позиция -> количество перекрывающих оснований.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     cov = r.calculate_coverage("chr1")
            ...     max(cov.values())
            15
        """
        coverage = {}
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            if rec.chrom == chrom and rec.cigar != "*":
                aligned_len = self._calc_aligned_length(rec.cigar)
                for i in range(rec.start, rec.start + aligned_len):
                    coverage[i] = coverage.get(i, 0) + 1

        if self.file:
            self.file.seek(current_pos)
        return coverage

    def filter_alignments(self, flag: int) -> Iterator[AlignmentRecord]:
        """
        Фильтрует выравнивания по флагу SAM.

        Args:
            flag: SAM-флаг для фильтрации (битовая маска)

        Yields:
            AlignmentRecord: Записи, у которых установлен указанный флаг.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     unmapped = list(r.filter_alignments(4))
            ...     print(len(unmapped))
        """

        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            if hasattr(rec, "flag") and rec.flag & flag:
                yield rec

        if self.file:
            self.file.seek(current_pos)

    def count_alignments(self) -> int:
        """Подсчитывает количество валидных выравниваний.

        Returns:
            int: Количество записей выравнивания.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     r.count_alignments() > 0
            True
        """
        cnt = 0
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for _ in self.read():
            cnt += 1

        if self.file:
            self.file.seek(current_pos)
        return cnt

    def stats_by_chromosome(self) -> pd.DataFrame:
        """
        Возвращает статистику по количества выравниваний на хромосому.

        Returns:
            pd.DataFrame: Таблица с колонками `chrom` и `count`.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     df = r.stats_by_chromosome()
            ...     list(df.columns)
            ['chrom', 'count']
        """
        from collections import Counter

        cnt = Counter()
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            cnt[rec.chrom] += 1

        if self.file:
            self.file.seek(current_pos)

        df = pd.DataFrame(list(cnt.items()), columns=["chrom", "count"])
        return df

    def filter_by_region(
        self, chrom: str, start: int, end: int
    ) -> Iterator[AlignmentRecord]:
        """
        Фильтрует выравнивания по геномному региону.
        
        Args:
            chrom: Название хромосомы
            start: Начальная позиция региона (1-based)
            end: Конечная позиция региона (1-based)
            
        Yields:
            AlignmentRecord: Записи, пересекающие регион.
        
        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     recs = list(r.filter_by_region("chr1", 1, 1000))
            ...     all(rec.chrom == "chr1" for rec in recs)
            True
        """
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            if rec.chrom == chrom and rec.end >= start and rec.start <= end:
                yield rec

        if self.file:
            self.file.seek(current_pos)

    def get_records_in_region(
        self, chrom: str, start: int, end: int
    ) -> List[AlignmentRecord]:
        """
        Получает все записи выравнивания в указанном регионе.
        
        Args:
            chrom: Название хромосомы
            start: Начальная позиция региона
            end: Конечная позиция региона
            
        Returns:
            List[AlignmentRecord]: Список объектов в регионе.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     len(r.get_records_in_region("chr1", 1, 10000)) > 0
            True
        """
        return list(self.filter_by_region(chrom, start, end))

    def filter_records(self, **filters) -> List[AlignmentRecord]:
        """
        Универсальный метод фильтрации записей по различным критериям.

        Поддерживаются фильтры `chrom`, `min_mapq` и `flag`.
        
        Args:
            **filters: Ключевые аргументы для фильтрации:
                - chrom: фильтр по хромосоме
                - min_mapq: минимальное качество маппинга
                - flag: SAM-флаг для фильтрации

        Returns:
            List[AlignmentRecord]: Список отфильтрованных записей AlignmentRecord
        
        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     recs = r.filter_records(chrom="chr1", min_mapq=30)
            ...     all(rec.mapq >= 30 for rec in recs)
            True
        """
        results = []
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            match = True

            if "chrom" in filters and rec.chrom != filters["chrom"]:
                match = False
            if "min_mapq" in filters and rec.mapq < filters["min_mapq"]:
                match = False
            if (
                "flag" in filters
                and hasattr(rec, "flag")
                and not (rec.flag & filters["flag"])
            ):
                match = False

            if match:
                results.append(rec)

        if self.file:
            self.file.seek(current_pos)
        return results