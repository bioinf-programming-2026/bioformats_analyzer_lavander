import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, List, Any
from abstract import GenomicReader
from record import VariantRecord

"""Ридер VCF-файлов.

Разбирает заголовок и варианты, предоставляет статистику,
фильтрацию, выборки по регионам и классификацию вариантов.

Examples:
    >>> with VcfReader("variants.vcf") as r:
    ...     print(r.get_statistics()["total_variants"])
    ...     for v in r.get_records_in_region("chr1", 1000, 2000):
    ...         print(v.id, v.ref, "->", v.alt)
"""

class VcfReader(GenomicReader):
    """
    Реализация ридера для файлов в формате VCF (Variant Call Format).

    VCF - стандартный формат для хранения геномных вариаций, содержащий
    информацию о хромосоме, позиции, референсном и альтернативном аллелях,
    качестве вызова и дополнительных аннотациях.

    Наследует от абстрактного класса `GenomicReader` и реализует
    специфичную для VCF логику парсинга и анализа данных.

    Attributes:
        headers (List[str]): Список строк заголовка VCF-файла в исходном порядке.
        header_groups (Dict[str, List[str]]): Словарь {тег: [строки]} — сгруппированный заголовок.
            Ключи содержат префикс: '##fileformat', '##INFO', '##FILTER',
            '##FORMAT', '##contig', '#CHROM'.
        column_names (List[str]): Список имён колонок из строки #CHROM.
        samples (List[str]): Список образцов (индивидуумов) из строки #CHROM.
        _variants_count: Счётчик вариантов в текущем проходе `read()`
            (для внутреннего использования).

    Examples:
        >>> with VcfReader("variants.vcf") as r:
        ...     stats = r.get_region_stats()
        ...     types = r.get_variant_type_stats()
        ...     print(types)
    """

    def __init__(self, filepath: str | Path):
        """
        Инициализирует VCF-ридер.

        Заголовок **не** разбирается здесь: файл ещё не открыт.
        Разбор произойдёт в `GenomicReader.__enter__` (или в `reset()`).

        Args:
            filepath: Путь к VCF-файлу.

        Examples:
            >>> VcfReader("variants.vcf").samples
            []
        """
        super().__init__(filepath)
        self.headers: List[str] = []
        self.header_groups: Dict[str, List[str]] = {}
        self.column_names: List[str] = []
        self.samples: List[str] = []
        self._variants_count = 0

    # ------------------------------------------------------------------
    #  Заголовок
    # ------------------------------------------------------------------

    def get_header(self) -> Dict[str, List[str]]:
        """
        Возвращает заголовок VCF, сгруппированный по тегам.

        Returns:
            Dict[str, List[str]]: ключи — теги с префиксом ('##fileformat',
            '##INFO', '#CHROM' и т.п.), значения — списки строк.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     "##fileformat" in r.get_header()
            True
        """
        return self.header_groups

    def get_header_group(self, group_name: str) -> List[str]:
        """
        Возвращает строки заголовка указанной группы.

        Имя группы можно передавать как с префиксом, так и без него:
        `'INFO'`, `'##INFO'` — эквивалентны; `'#CHROM'` — как есть.

        Args:
            group_name: Имя группы.

        Returns:
            List[str]: Список строк указанной группы или пустой список.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     for line in r.get_header_group("INFO")[:2]:
            ...         print(line)
        """
        key = group_name if group_name.startswith("#") else f"##{group_name}"
        return self.header_groups.get(key, [])

    def _parse_header(self):
        """
        Разбирает заголовок VCF, валидирует формат и извлекает группы.

        - Проверяет, что первая строка начинается с '##fileformat=VCF'.
        - Проверяет наличие строки '#CHROM'.
        - Заполняет `self.headers`, `self.header_groups`, `self.column_names`,
          `self.samples`.
        - Устанавливает `self._header_parsed = True`.

        Raises:
            ValueError: Если файл не соответствует формату VCF
                (нет '##fileformat=VCF' или отсутствует '#CHROM').

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     r._parse_header()
            ...     len(r.headers) > 0
            True
        """
        self.headers = []
        self.header_groups = {}
        self.column_names = []
        self.samples = []

        if not self.file:
            return

        # Сохраняем начальную позицию для восстановления после парсинга
        current_pos = self.file.tell()
        self.file.seek(0)

        try:
            first_line = self.file.readline()
            if not first_line.startswith("##fileformat=VCF"):
                raise ValueError(
                    f"Файл '{self.filepath}' не соответствует формату VCF! "
                    "Первая строка должна начинаться с '##fileformat=VCF'."
                )

            self.headers.append(first_line.strip())
            group_key = first_line.split("=", 1)[0]          # '##fileformat'
            self.header_groups.setdefault(group_key, []).append(first_line.strip())

            line = self.file.readline()
            while line:
                line_str = line.strip()
                if line_str.startswith("##"):
                    self.headers.append(line_str)
                    group_key = line_str.split("=", 1)[0]    # '##INFO', '##FILTER'
                    self.header_groups.setdefault(group_key, []).append(line_str)
                elif line_str.startswith("#CHROM"):
                    # Парсим строку с названиями колонок
                    # Формат: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SAMPLE1 ...]
                    self.headers.append(line_str)
                    parts = line_str.split("\t")
                    self.column_names = [p.lstrip("#") for p in parts]
                    if len(parts) >= 10:  # Минимум: обязательные поля + samples
                        self.samples = parts[9:]  # Образцы начинаются с 9-й колонки
                    self.header_groups.setdefault("#CHROM", []).append(line_str)
                    self._header_parsed = True
                    return
                else:
                    raise ValueError(
                        f"Файл '{self.filepath}': отсутствует строка "
                        f"заголовка столбцов '#CHROM'."
                    )
                line = self.file.readline()

            raise ValueError(
                f"Файл '{self.filepath}': отсутствует строка "
                f"заголовка столбцов '#CHROM'."
            )
        finally:
            # Возвращаемся к начальной позиции чтения данных
            self.file.seek(current_pos)

    # ------------------------------------------------------------------
    #  Чтение
    # ------------------------------------------------------------------

    def read(self) -> Iterator[VariantRecord]:
        """
        Основной метод чтения вариантов из VCF-файла.

        Генератор, который последовательно читает и парсит варианты из файла,
        возвращая объекты VariantRecord для каждой валидной строки с вариантом.
        Некорректные строки (меньше 8 колонок, невалидный POS/QUAL) пропускаются
        с выводом предупреждения — поток не прерывается.

        Yields:
            VariantRecord: Объект, содержащий информацию о геномном варианте

        Raises:
            ValueError: Если файл не был открыт перед чтением.

        Процесс:
        1. Пропускает все заголовочные строки
        2. Читает файл построчно
        3. Для каждой не-пустой и не-заголовочной строки вызывает _parse_variant_line()
        4. Увеличивает счетчик вариантов для каждого успешно распарсенного варианта

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     first = next(r.read())
            ...     first.chrom
            'chr1'
        """
        if not self.file:
            raise ValueError(
                "Файл не открыт. Используйте VcfReader как контекстный менеджер."
            )
        
        # Переходим к началу данных (после заголовков)
        self.file.seek(0)
        line = self.file.readline()
        while line and line.startswith('#'):
            line = self.file.readline()

        # Читаем и парсим данные о вариантах
        self._variants_count = 0
        while line:
            line = line.strip()
            if line and not line.startswith('#'):
                variant = self._parse_variant_line(line)
                if variant:
                    self._variants_count += 1
                    yield variant
            line = self.file.readline()

    def _parse_variant_line(self, line: str) -> VariantRecord | None:
        """
        Парсит одну строку с вариантом из VCF-файла.

        VCF-строки вариантов содержат минимум 8 обязательных полей:
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO. Здесь выполняется
        вся нормализация и приведение типов:

        - CHROM, ID, FILTER — `.strip()`;
        - POS — `int()`;
        - REF — `.strip().upper()`;
        - ALT — `.strip().upper()`, split по `,`, пустой список при `"."`;
        - QUAL — `float()`, `None` при `"."`;
        - INFO — разбор `k1=v1;k2;k3=v3` в `dict`.

        Готовые значения передаются в конструктор `VariantRecord`.

        Args:
            line (str): Строка из VCF-файла, содержащая описание варианта.

        Returns:
            VariantRecord | None: Объект варианта или None при ошибке парсинга
            (короткая строка, невалидный POS/QUAL и т.п.).

        Examples:
            >>> r = VcfReader("variants.vcf")
            >>> rec = r._parse_variant_line(
            ...     "chr1\\t100\\t.\\tA\\tT\\t50\\tPASS\\tDP=30")
            >>> rec.pos, rec.ref, rec.alt, rec.qual, rec.info["DP"]
            (100, 'A', ['T'], 50.0, '30')
            >>> rec = r._parse_variant_line(
            ...     "chr1\\t100\\t.\\tA\\t.\\t.\\tPASS\\tDB;DP=5")
            >>> rec.alt, rec.qual, rec.info
            ([], None, {'DB': True, 'DP': '5'})
        """
        parts = line.split('\t')
        if len(parts) < 8:
            return None

        try:
            chrom = parts[0].strip()
            pos = int(parts[1])
            var_id = parts[2].strip()
            ref = parts[3].strip().upper()

            # ALT может содержать несколько аллелей через запятую
            alt_raw = parts[4]
            alt: List[str] = (
                [a.strip().upper() for a in alt_raw.split(",")]
                if alt_raw != "." else []
            )

            # QUAL может быть числом или точкой '.'
            qual_raw = parts[5]
            qual = float(qual_raw) if qual_raw != "." else None

            filt = parts[6].strip()

            # Парсим INFO поле в словарь
            info_str = parts[7]
            info: Dict[str, Any] = {}
            if info_str and info_str != ".":
                for item in info_str.strip().split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        info[k] = v
                    elif item:
                        info[item] = True

            return VariantRecord(
                chrom=chrom,
                pos=pos,
                var_id=var_id,
                ref=ref,
                alt=alt,
                qual=qual,
                filt=filt,
                info=info,
            )
        except (ValueError, IndexError) as e:
            print(f"Ошибка парсинга варианта: {e}")
            return None

    # ------------------------------------------------------------------
    #  Регион / пересечение
    # ------------------------------------------------------------------

    def get_records_in_region(
        self, chrom: str, start: int, end: int
    ) -> Iterator[VariantRecord]:
        """
        Возвращает варианты, пересекающие указанный регион [start, end].

        Использует `VariantRecord.intersects` (учёт длины REF для INDEL).
        Реализует логику bedtools intersect. Работает как генератор,
        фильтруя поток мутаций.

        Args:
            chrom: Название хромосомы.
            start: Начало региона (1-based, включительно).
            end: Конец региона (1-based, включительно).

        Yields:
            VariantRecord: Варианты, пересекающие регион.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     recs = list(r.get_records_in_region("chr1", 100, 500))
            ...     all(v.chrom == "chr1" for v in recs)
            True
        """
        for record in self.read():
            if record.intersects(chrom, start, end):
                yield record

    # ------------------------------------------------------------------
    #  Фильтрация и подсчёт
    # ------------------------------------------------------------------

    def filter_records(self, **filters) -> Iterator[VariantRecord]:
        """Фильтрует варианты по атрибутам VariantRecord.

        Args:
            **filters: Имена атрибутов и ожидаемые значения.

        Yields:
            VariantRecord: Варианты, удовлетворяющие фильтрам.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     snvs = list(r.filter_records(chrom="chr1"))
            ...     print(len(snvs))
        """
        for record in self.read():
            if all(getattr(record, k, None) == v for k, v in filters.items()):
                yield record

    def get_variants_count(self) -> int:
        """
        Потоковый подсчёт общего количества вариантов в файле.

        Итерируется по `read()`, не загружая данные в память целиком.
        Считаются только валидные варианты (те же, что отдаёт `read()`).

        Returns:
            int: Общее число валидных вариантов.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     r.get_variants_count() >= 0
            True
        """
        return sum(1 for _ in self.read())

    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """Проверяет корректность геномной координаты.

        Args:
            chrom: Название хромосомы.
            pos: Позиция (1-based).

        Returns:
            bool: True, если хромосома есть в файле и позиция положительна.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     r.validate_coordinate("chr1", 100)
            True
        """
        return chrom in self.get_chromosomes() and pos > 0

    # ------------------------------------------------------------------
    #  Статистика
    # ------------------------------------------------------------------

    def get_statistics(self) -> Dict[str, Any]:
        """
        Возвращает базовую статистику по VCF-файлу.

        Returns:
            Dict[str, Any]: Словарь со статистикой:
                - total_variants: общее количество вариантов
                - samples_count: количество образцов
                - header_lines: количество строк в заголовке

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     stats = r.get_statistics()
            ...     "total_variants" in stats
            True
        """
        return {
            'total_variants': self.get_variants_count(),
            'samples_count': len(self.samples),
            'header_lines': len(self.headers),
        }

    def get_chromosomes(self) -> List[str]:
        """
        Возвращает список уникальных хромосом, представленных в файле.

        Returns:
            List[str]: Отсортированный список названий хромосом.

        Примечание:
            Проходит через `read()`, поэтому учитываются только валидные
            варианты (битые строки пропускаются, как и при чтении).

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     "chr1" in r.get_chromosomes()
            True
        """
        return sorted({record.chrom for record in self.read()})

    def get_region_stats(self) -> pd.DataFrame:
        """
        Возвращает распределение вариантов по хромосомам.

        Returns:
            pd.DataFrame: Колонки ['chromosome', 'variant_count'].

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     df = r.get_region_stats()
            ...     list(df.columns)
            ['chromosome', 'variant_count']
        """
        region_counts: Dict[str, int] = {}

        for record in self.read():
            region_counts[record.chrom] = region_counts.get(record.chrom, 0) + 1

        # Создаем DataFrame из собранных данных
        return pd.DataFrame(
            [{'chromosome': c, 'variant_count': n}
             for c, n in region_counts.items()],
            columns=['chromosome', 'variant_count'],
        )

    def get_variant_type_stats(self) -> pd.DataFrame:
        """
        Анализирует и возвращает статистику по типам геномных вариантов.

        Классификация: SNV, Insertion, Deletion, Complex,
        MultiAllelic, Symbolic.

        Returns:
            pd.DataFrame: Колонки ['variant_type', 'count'].

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     df = r.get_variant_type_stats()
            ...     set(df["variant_type"]).issubset(
            ...         {"SNV", "Insertion", "Deletion", "Complex",
            ...          "MultiAllelic", "Symbolic"})
            True
        """
        type_counts = {
            'SNV': 0,
            'Insertion': 0,
            'Deletion': 0,
            'Complex': 0,
            'MultiAllelic': 0,
            'Symbolic': 0,
        }

        for record in self.read():
            alts = record.alt
            if len(alts) > 1:
                type_counts['MultiAllelic'] += 1
            elif not alts:
                # ALT == "." в исходной строке
                type_counts['Symbolic'] += 1
            else:
                a = alts[0]
                if a in ('.', '*') or a.startswith('<'):
                    type_counts['Symbolic'] += 1
                elif len(record.ref) == 1 and len(a) == 1:
                    type_counts['SNV'] += 1
                elif len(record.ref) < len(a):
                    type_counts['Insertion'] += 1
                elif len(record.ref) > len(a):
                    type_counts['Deletion'] += 1
                else:
                    type_counts['Complex'] += 1

        # Создаем DataFrame только для ненулевых типов
        return pd.DataFrame(
            [{'variant_type': t, 'count': c}
             for t, c in type_counts.items() if c > 0],
            columns=['variant_type', 'count'],
        )

    # ------------------------------------------------------------------
    #  Управление состоянием
    # ------------------------------------------------------------------

    def reset(self) -> None:
        """
        Сбрасывает состояние ридера к начальному.

        Действия:
            - Перемещает позицию чтения файла в начало
            - Перепарсивает заголовок файла
            - Сбрасывает внутренние счетчики и состояния

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     list(r.read())
            ...     r.reset()
            ...     len(list(r.read())) >= 0
            True
        """
        if self.file:
            self.file.seek(0)
            self._parse_header()

    def close(self) -> None:
        """Закрывает файл и сбрасывает VCF-специфичное состояние.

        Examples:
            >>> r = VcfReader("variants.vcf")
            >>> r.__enter__()
            >>> r.close()
            >>> r.headers, r.samples
            ([], [])
        """
        super().close()
        self.headers = []
        self.header_groups = {}
        self.column_names = []
        self.samples = []
        self._variants_count = 0