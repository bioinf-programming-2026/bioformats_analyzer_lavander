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
    
    Наследует от абстрактного класса GenomicDataReader и реализует
    специфичную для VCF логику парсинга и анализа данных.

    Attributes:
        headers (List[str]): Список строк заголовка VCF файла
        samples (List[str]): Список образцов (индивидуумов) из строки #CHROM VCF файла
        _variants_count (int): Счетчик вариантов (для внутреннего использования)

    Examples:
        >>> with VcfReader("variants.vcf") as r:
        ...     stats = r.get_region_stats()
        ...     types = r.get_variant_type_stats()
        ...     print(types)
    """     

    def __init__(self, filepath: str | Path):
        """
        Инициализирует VCF-ридер.

        Args:
            filepath: Путь к VCF-файлу.

        Examples:
            >>> VcfReader("variants.vcf").samples
            []
        """
        super().__init__(filepath)
        self.headers = []
        self.samples = []
        self._variants_count = 0

    def get_header(self) -> Dict[str, List[str]]:
        """
        Возвращает заголовок VCF, сгруппированный по тегам.

        Returns:
            Dict[str, List[str]]: ключи — теги (включая '##' или '#'),
            значения — списки строк с этим тегом (без изменений).

        Пример:
            {
                '##fileformat': ['##fileformat=VCFv4.2'],
                '##fileDate':   ['##fileDate=20250101'],
                '##reference':  ['##reference=GRCh38'],
                '##source':     ['##source=myPipeline'],
                '##INFO': [
                    '##INFO=<ID=DP,...>',
                    '##INFO=<ID=AF,...>',
                    '##INFO=<ID=DB,...>',
                ],
                '#CHROM': ['#CHROM\tPOS\tID\t...'],
            }
        
        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     "##fileformat" in r.get_header()
            True
        """
        grouped: Dict[str, List[str]] = {}
        for line in self.headers:
            if line.startswith("##"):
                tag = line.split("=", 1)[0]          # '##INFO', '##fileformat'
            elif line.startswith("#"):
                tag = line.split("\t", 1)[0]         # '#CHROM'
            else:
                continue
            grouped.setdefault(tag, []).append(line)
        return grouped

    def _parse_header(self):
        """
        Разбирает заголовок VCF и извлекает имена образцов.
        
        VCF заголовок состоит из строк, начинающихся с '#' и содержит:
        - Мета-информацию (##-строки)
        - Строку с названиями колонок (#CHROM...)
        - Информацию о формате и образцах
        
        Алгоритм:
        1. Сохраняет текущую позицию в файле
        2. Читает файл с начала до первой не-заголовочной строки
        3. Извлекает имена образцов из строки #CHROM
        4. Восстанавливает исходную позицию чтения
        
        Устанавливает:
            self.headers: все строки заголовка
            self.samples: список образцов (если есть)
            self._header_parsed: True после завершения парсинга

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     r._parse_header()
            ...     len(r.headers) > 0
            True
        """
        self.headers = []
        self.samples = []
        
        if not self.file:
            return

        # Сохраняем начальную позицию для восстановления после парсинга
        current_pos = self.file.tell()
        self.file.seek(0)
        
        line = self.file.readline()
        while line and line.startswith('#'):
            self.headers.append(line.strip())
            
            if line.startswith('#CHROM'):
                # Парсим строку с названиями колонок
                # Формат: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SAMPLE1 ...]
                parts = line.strip().split('\t')
                if len(parts) >= 10:  # Минимум: обязательные поля + samples
                    self.samples = parts[9:]  # Образцы начинаются с 9-й колонки
            
            line = self.file.readline()
        
        # Возвращаемся к начальной позиции чтения данных
        self.file.seek(current_pos)
        self._header_parsed = True

    def validate_coordinate(self, chrom, pos):
        """Проверяет корректность геномной координаты.

        Args:
            chrom: Название хромосомы.
            pos: Позиция.

        Returns:
            bool: True, если хромосома есть в файле и позиция положительна.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     r.validate_coordinate("chr1", 100)
            True
        """
        return chrom in self.get_chromosomes() and pos > 0
    
    def get_records_in_region(
        self, chrom: str, start: int, end: int
        ) -> Iterator[VariantRecord]:
        """Возвращает варианты, пересекающие указанный регион.

        Интервал варианта считается как `[pos, pos + len(ref) - 1]`.

        Args:
            chrom: Название хромосомы.
            start: Начало региона.
            end: Конец региона.

        Yields:
            VariantRecord: Варианты, пересекающие регион.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     recs = list(r.get_records_in_region("chr1", 100, 500))
            ...     all(v.chrom == "chr1" for v in recs)
            True
        """
        for record in self.read():
            if record.chrom != chrom:
                continue
            rec_start = record.pos
            rec_end = record.pos + len(record.ref) - 1
            if rec_start <= end and rec_end >= start:
                yield record
    
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


    def read(self) -> Iterator[VariantRecord]:
        """
        Основной метод чтения вариантов из VCF файла.
        
        Генератор, который последовательно читает и парсит варианты из файла,
        возвращая объекты VariantRecord для каждой валидной строки с вариантом.
        
        Yields:
            VariantRecord: Объект, содержащий информацию о геномном варианте
            
        Raises:
            ValueError: Если файл не был открыт перед чтением
            
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
            raise ValueError("Файл не открыт")
        
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
        Парсит одну строку с вариантом из VCF файла.
        
        VCF строки вариантов содержат минимум 8 обязательных полей:
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
        
        Args:
            line (str): Строка из VCF файла, содержащая описание варианта
            
        Returns:
            VariantRecord | None: Объект варианта или None при ошибке парсинга
            
        Обрабатываемые поля:
            - CHROM: хромосома (str)
            - POS: позиция (int)
            - REF: референсный аллель (str)
            - ALT: альтернативный аллель (str)
            - INFO: поле с дополнительной информацией (парсится в dict)
        
        Examples:
            >>> r = VcfReader("variants.vcf")
            >>> rec = r._parse_variant_line("chr1\\t100\\t.\\tA\\tT\\t50\\tPASS\\tDP=30")
            >>> rec.ref, rec.alt, rec.info["DP"]
            ('A', 'T', '30')
        """
        parts = line.split('\t')
        if len(parts) < 8:
            return None
        
        try:
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            
            # Парсим INFO поле в словарь
            info_dict = {}
            if parts[7] != '.':
                for info_item in parts[7].split(';'):
                    if '=' in info_item:
                        key, value = info_item.split('=', 1)
                        info_dict[key] = value
                    else:
                        info_dict[info_item] = True
            
            return VariantRecord(chrom, pos, ref, alt, info_dict)
            
        except (ValueError, IndexError) as e:
            print(f"Ошибка парсинга варианта: {e}")
            return None

    def get_statistics(self) -> Dict[str, Any]:
        """
        Возвращает базовую статистику по VCF файлу.
        
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
        # Пересчитываем количество вариантов для актуальности
        count = self._count_variants()
        return {
            'total_variants': count,
            'samples_count': len(self.samples),
            'header_lines': len(self.headers)
        }

    def _count_variants(self) -> int:
        """
        Подсчитывает количество валидных вариантов в файле.
        
        Returns:
            int: Количество строк с вариантами (не-заголовочных строк)
            
        Примечание:
            Метод независимо подсчитывает варианты, временно перемещая
            позицию чтения файла и восстанавливая её после подсчета.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     r._count_variants() >= 0
            True
        """
        if not self.file:
            return 0
        
        current_pos = self.file.tell()
        self.file.seek(0)
        
        # Пропускаем заголовки
        line = self.file.readline()
        while line and line.startswith('#'):
            line = self.file.readline()
        
        # Считаем варианты (все не-пустые не-заголовочные строки)
        count = 0
        while line:
            line = line.strip()
            if line and not line.startswith('#'): 
                parts = line.split('\t')
                if len(parts) >= 8:
                    try:
                        int(parts[1])
                        count += 1
                    except ValueError:
                        pass
            line = self.file.readline()
        
        # Восстанавливаем исходную позицию чтения
        self.file.seek(current_pos)
        return count

    def get_chromosomes(self) -> List[str]:
        """
        Возвращает список уникальных хромосом, представленных в файле.
        
        Returns:
            List[str]: Отсортированный список названий хромосом
            
        Примечание:
            Для больших файлов может быть затратной операцией,
            так требует полного прохода по файлу.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     "chr1" in r.get_chromosomes()
            True
        """
        chromosomes = set()
        
        if not self.file:
            return list(chromosomes)
        
        current_pos = self.file.tell()
        self.file.seek(0)
        
        # Пропускаем заголовки
        line = self.file.readline()
        while line and line.startswith('#'):
            line = self.file.readline()
        
        # Собираем уникальные хромосомы
        while line:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 8:
                    chromosomes.add(parts[0])  # CHROM - первая колонка
            line = self.file.readline()
        
        # Восстанавливаем позицию чтения
        self.file.seek(current_pos)
        
        return sorted(list(chromosomes))

    def get_region_stats(self) -> pd.DataFrame:
        """
        Возвращает распределение вариантов по хромосомам.
        
        Returns:
            pd.DataFrame: DataFrame с колонками:
                - chromosome: название хромосомы
                - variant_count: количество вариантов на хромосоме
        
        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     df = r.get_region_stats()
            ...     list(df.columns)
            ...     print(stats_df)  # Показывает таблицу с распределением
            ['chromosome', 'variant_count']
                    
        """
        region_counts = {}
        
        if not self.file:
            return pd.DataFrame(columns=['chromosome', 'variant_count'])
        
        current_pos = self.file.tell()
        self.file.seek(0)
        
        # Пропускаем заголовки
        line = self.file.readline()
        while line and line.startswith('#'):
            line = self.file.readline()
        
        # Считаем варианты по хромосомам
        while line:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split('\t')
                if parts:
                    chrom = parts[0]
                    region_counts[chrom] = region_counts.get(chrom, 0) + 1
            line = self.file.readline()
        
        # Восстанавливаем позицию чтения
        self.file.seek(current_pos)
        
        # Создаем DataFrame из собранных данных
        data = []
        for chrom, count in region_counts.items():
            data.append({'chromosome': chrom, 'variant_count': count})
        
        return pd.DataFrame(data)

    def get_variant_type_stats(self) -> pd.DataFrame:
        """
        Анализирует и возвращает статистику по типам геномных вариантов.
        
        Классификация вариантов:
            - SNV (Single Nucleotide Variant): одинарные нуклеотидные замены
            - Insertion: вставки (REF короче ALT)
            - Deletion: делеции (REF длиннее ALT)  
            - Complex: сложные варианты (равная длина, но разные последовательности)
            - MultiAllelic
            - Symbolic
            
        Returns:
            pd.DataFrame: DataFrame с колонками:
                - variant_type: тип варианта
                - count: количество вариантов этого типа
        
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
        
        if not self.file:
            return pd.DataFrame(columns=['variant_type', 'count'])
        
        current_pos = self.file.tell()
        self.file.seek(0)
        
        # Пропускаем заголовки
        line = self.file.readline()
        while line and line.startswith('#'):
            line = self.file.readline()
        
        # Анализируем типы вариантов по длине REF и ALT
        while line:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 5:
                    ref = parts[3]
                    alt = parts[4]
                    alts = alt.split(',')

                    # Определяем тип варианта по длине последовательностей
                    if len(alts) > 1:
                        type_counts['MultiAllelic'] += 1
                    else:
                        a = alts[0]
                        if a in ('.', '*') or a.startswith('<'):
                            type_counts['Symbolic'] += 1
                        elif len(ref) == 1 and len(alt) == 1:
                            type_counts['SNV'] += 1
                        elif len(ref) < len(alt):
                            type_counts['Insertion'] += 1
                        elif len(ref) > len(alt):
                            type_counts['Deletion'] += 1
                        else:
                            type_counts['Complex'] += 1
            
            line = self.file.readline()
        
        # Восстанавливаем позицию чтения
        self.file.seek(current_pos)
        
        # Создаем DataFrame только для ненулевых типов
        data = []
        for var_type, count in type_counts.items():
            if count > 0:
                data.append({'variant_type': var_type, 'count': count})
        
        return pd.DataFrame(data)

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
        self.samples = []
        self._variants_count = 0