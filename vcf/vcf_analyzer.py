from typing import Iterator, List, Dict
import pandas as pd
from abstract import GenomicDataReader, VcfParsingError
from record import VariantRecord


class VcfAnalyzer(GenomicDataReader):
    """Потоковый парсер и анализатор файлов формата VCF."""

    def __init__(self, filename: str):
        super().__init__(filename)
        self.raw_headers: List[str] = []
        self.header_groups: Dict[str, List[str]] = {}
        self.column_names: List[str] = []
        self._parse_headers()

    def _parse_headers(self):
        """
        Парсит метаданные VCF (строки ##) и группирует их.
        Проверяет обязательную спецификацию файла.
        """
        self.open()
        self._file.seek(0)
        self.raw_headers = []
        self.header_groups = {}

        first_line = self._file.readline()
        if not first_line.startswith("##fileformat=VCF"):
            raise VcfParsingError(
                f"Файл '{self.filename}' не соответствует формату VCF! "
                "Первая строка должна начинаться с '##fileformat=VCF'."
            )
        self.raw_headers.append(first_line.strip())

        for line in self._file:
            line_str = line.strip()
            if line_str.startswith("##"):
                self.raw_headers.append(line_str)
                # Извлекаем ключ группы: ##INFO=<...> -> 'INFO'
                group_key = line_str[2:].split("=", 1)[0]
                self.header_groups.setdefault(group_key, []).append(line_str)
            elif line_str.startswith("#CHROM"):
                self.raw_headers.append(line_str)
                self.column_names = line_str[1:].split("\t")
                break
            else:
                raise VcfParsingError("Отсутствует строка заголовка столбцов '#CHROM'.")

    def get_header(self) -> List[str]:
        """Требование ТЗ: Получение полного заголовка файла."""
        return self.raw_headers

    def get_header_group(self, group_name: str) -> List[str]:
        """
        Требование ТЗ: Получение информации по отдельным группам заголовков
        (например, 'INFO', 'FILTER', 'FORMAT', 'contig').
        """
        return self.header_groups.get(group_name, [])

    def read(self) -> Iterator[VariantRecord]:
        """
        Ленивый генератор: построчно считывает варианты (мутации).
        Не хранит все варианты в оперативной памяти.
        """
        self.open()
        self._file.seek(0)

        line_num = 0
        for line in self._file:
            line_num += 1
            if line.startswith("#") or not line.strip():
                continue

            fields = line.strip().split("\t")
            if len(fields) < 8:
                raise VcfParsingError(
                    f"Ошибка в строке {line_num}: запись VCF должна содержать "
                    f"минимум 8 обязательных колонок, получено: {len(fields)}"
                )

            try:
                yield VariantRecord(
                    chrom=fields[0],
                    pos=int(fields[1]),
                    var_id=fields[2],
                    ref=fields[3],
                    alt=fields[4],
                    qual=fields[5],
                    filt=fields[6],
                    info_str=fields[7]
                )
            except ValueError as e:
                raise VcfParsingError(f"Ошибка парсинга позиции/качества в строке {line_num}: {e}")

    def get_variants_count(self) -> int:
        """Требование ТЗ: Получение общего количества вариантов."""
        count = 0
        for _ in self.read():
            count += 1
        return count

    def get_region_statistics(self) -> pd.DataFrame:
        """
        Требование ТЗ: Получение статистики 'количество вариантов - регион'
        (с использованием pandas).
        """
        
        chrom_stream = (record.chrom for record in self.read())
        df = pd.DataFrame(chrom_stream, columns=["CHROM"])

        if df.empty:
            return pd.DataFrame(columns=["CHROM", "variant_count"])

        
        stats = (
            df["CHROM"]
            .value_counts()
            .reset_index()
        )
        stats.columns = ["CHROM", "variant_count"]
        return stats

    def intersect(self, chrom: str, start: int, end: int) -> Iterator[VariantRecord]:
        """
        Требование ТЗ: Получение вариантов, лежащих в определенном геномном отрезке
        (аналог bedtools intersect).
        Работает как генератор, фильтруя поток мутаций.
        """
        for record in self.read():
            if record.intersects(chrom, start, end):
                yield record