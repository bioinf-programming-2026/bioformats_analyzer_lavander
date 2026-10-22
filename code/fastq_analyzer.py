from pathlib import Path
from typing import Iterator
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from abstract import SequenceReader
from record import SequenceRecord

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "test" / "plots"

"""Анализатор FASTQ-файлов.

Читает последовательности и Phred-качество, считает базовую статистику
и строит диагностические графики.

Examples:
    >>> with FastqAnalyzer("reads.fastq") as a:
    ...     list(a.read())
    ...     a.get_seq_count()
    1000
"""

class FastqAnalyzer(SequenceReader):
    """
    Анализатор FASTQ файлов для биоинформатического анализа.
    
    Наследуется от абстрактного класса SequenceReader и предоставляет функциональность
    для чтения и анализа данных из файлов в формате FASTQ.

    Поддерживает чтение FASTQ-записей, подсчёт количества и средней длины,
    конвертацию Phred-качества и построение графиков.
    
    FASTQ формат содержит последовательности ДНК/РНК с оценками качества:
    - Строка 1: Заголовок (начинается с @)
    - Строка 2: Нуклеотидная последовательность
    - Строка 3: Разделитель (обычно +)
    - Строка 4: Строка качества в формате Phred
    
    Attributes:
        filepath (Path): Путь к анализируемому FASTQ файлу
        _seq_count (int): Количество прочитанных последовательностей
        _total_length (int): Общая длина всех последовательностей

    Examples:
        >>> with FastqAnalyzer("reads.fastq") as a:
        ...     a.per_base_sequence_quality("q.png")
        ...     a.per_base_sequence_content("c.png")
        ...     a.sequence_length_distribution("l.png")
    """

    def __init__(self, filepath: str | Path):
        """
        Инициализирует анализатор FASTQ файлов.
        Args:
            filepath: Путь к FASTQ-файлу в виде строки или объекта Path.

        Examples:
            >>> FastqAnalyzer("reads.fastq")._seq_count
            0
        """
        super().__init__(filepath)
        self._seq_count = 0
        self._total_length = 0

    def read(self) -> Iterator[SequenceRecord]:
        """
        Читает последовательности из FASTQ файла и возвращает записи последовательностей.

        Сбрасывает внутреннюю статистику, читает файл блоками по четыре строки
        и создаёт SequenceRecord с последовательностью и качеством.

        Генератор, который читает файл построчно, группирует строки по 4 (полная запись FASTQ)
        и создает объекты SequenceRecord для каждой последовательности.

        Yields:
            SequenceRecord: Очередная FASTQ-запись с информацией о последовательности (ID и последовательность)

        Examples:
            >>> with FastqAnalyzer("reads.fastq") as a:
            ...     rec = next(a.read())
            ...     rec.id
            'read1'
            
        Пример FASTQ записи:
            @seq1
            ATGCGATCGATCG
            +
            IIIIIIIIIIIII
        """
        self._seq_count = 0
        self._total_length = 0

        with open(self.filepath, "r", encoding="utf-8") as file:
            lines = [file.readline().strip() for _ in range(4)]
            while lines[0]:
                header = lines[0][1:]  # Убираем символ @ из заголовка
                sequence = lines[1]
                quality = [self.phred_to_quality(c) for c in lines[3]] # Строка качества
                record = SequenceRecord(id=header, sequence=sequence, quality=quality)
                    
                # Обновляем статистику
                self._seq_count += 1
                self._total_length += len(sequence)
                    
                yield record
                lines = [file.readline().strip() for _ in range(4)]  # Сбрасываем буфер для следующей записи

    def get_sequence(self, id: str, seq: str) -> SequenceRecord:
        """Создаёт объект SequenceRecord.

        Args:
            id: Идентификатор последовательности.
            seq: Нуклеотидная последовательность.

        Returns:
            SequenceRecord: Созданная запись.

        Examples:
            >>> FastqAnalyzer("x.fastq").get_sequence("r", "AC").sequence
            'AC'
        """
        return SequenceRecord(id=id, sequence=seq)
    
    def validate_sequence(self, seq: str) -> bool:
        """Проверяет, что последовательность непустая и состоит из букв.

        Args:
            seq: Проверяемая последовательность.

        Returns:
            bool: Результат проверки.

        Examples:
            >>> FastqAnalyzer("x.fastq").validate_sequence("ACGT")
            True
        """
        return bool(seq) and seq.isalpha()
    
    def get_seq_count(self) -> int:
        """Возвращает количество прочитанных последовательностей.

        Returns:
            int: Количество последовательностей.

        Examples:
            >>> with FastqAnalyzer("reads.fastq") as a:
            ...     list(a.read())
            ...     a.get_seq_count() > 0
            True
        """
        return self._seq_count

    def get_mean_seq_length(self) -> float:
        """Возвращает среднюю длину последовательностей.

        Returns:
            float: Средняя длина последовательности; 0.0, если последовательностей нет.

        Examples:
            >>> with FastqAnalyzer("reads.fastq") as a:
            ...     list(a.read())
            ...     round(a.get_mean_seq_length(), 1)
            5.0
        """
        if self._seq_count == 0:
            return 0.0
        return self._total_length / self._seq_count

    def phred_to_quality(self, phred_char: str, offset: int = 33) -> int:
        """
        Конвертирует символ Phred в числовое значение качества.
        
        Phred качество кодируется одним символом ASCII. По умолчанию используется
        кодировка Sanger (offset=33), где качество = ASCII код - 33.
        
        Args:
            phred_char: Один символ из строки качества
            offset: Смещение для декодирования (33 для Sanger, 64 для Illumina 1.3+)
            
        Returns:
            int: Численное значение качества от 0 до 40+

        Raises:
            ValueError: Если символ качества выходит за допустимый диапазон

        Пример:
            'I' (ASCII 73) -> 73 - 33 = 40

        Examples:
            >>> FastqAnalyzer("x.fastq").phred_to_quality('I')
            40
            >>> FastqAnalyzer("x.fastq").phred_to_quality('h', offset=64)
            40
        """
        q = ord(phred_char) - offset
        if not 0 <= q <= 93:
            raise ValueError(f"Недопустимый символ качества: {phred_char!r}")
        return q

    def get_sequences_with_quality(self) -> Iterator[tuple[str, str, list[int]]]:
        """
        Генератор для чтения последовательностей вместе с информацией о качестве.
        
        Читает полные записи FASTQ и возвращает кортежи с заголовком, 
        последовательностью и строкой качества.
        
        Yields:
            tuple[str, str, list[int]]: Идентификатор, последовательность и список качеств для каждой записи FASTQ
            
        Используется для анализа качества последовательностей.

        Examples:
            >>> with FastqAnalyzer("reads.fastq") as a:
            ...     rid, seq, q = next(a.get_sequences_with_quality())
            ...     len(seq) == len(q)
            True
        """
        for record in self.read():
            yield record.id, record.sequence, record.quality or []

    def per_base_sequence_quality(self, output_file: str = "fastq_quality_plot.png"):
        """
        Строит график качества оснований по позициям в ридах.
        
        Анализирует качество каждого основания в каждой позиции рида и
        визуализирует среднее качество по всем последовательностям.
        
        Args:
            output_file: Имя файла для сохранения графика
            
        Returns:
            Path | str: Путь к сохранённому графику или имя файла при отсутствии данных.
        
        Examples:
            >>> with FastqAnalyzer("reads.fastq") as a:
            ...     path = a.per_base_sequence_quality("quality.png")
            ...     path.exists()
            True
        """
        # Словарь для хранения качества по позициям: {позиция: [качества]}
        quality_by_position = defaultdict(list)
        
        # Сбор данных о качестве для каждой позиции
        for _, sequence, quality in self.get_sequences_with_quality():
            for pos, q in enumerate(quality):
                quality_by_position[pos].append(q)
        
        # Проверка наличия данных
        if not quality_by_position:
            print("Нет данных для графика качества")
            return output_file
            
        # Подготовка данных для графика
        positions = sorted(quality_by_position.keys())
        mean_qualities = [np.mean(quality_by_position[pos]) for pos in positions]
        
        # Построение графика
        plt.figure(figsize=(12, 6))
        plt.plot(positions, mean_qualities, linewidth=2)
        plt.xlabel('Position in read (bp)')
        plt.ylabel('Quality score')
        plt.title('Per Base Sequence Quality')
        plt.grid(True, alpha=0.3)
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        output_path = OUTPUT_DIR / output_file
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"График качества создан: {output_path}")
        return output_path

    def per_base_sequence_content(self, output_file: str = "fastq_content_plot.png"):
        """
        Анализирует и визуализирует содержание нуклеотидов по позициям.
        
        Вычисляет процентное содержание A, T, G, C оснований для каждой позиции
        во всех последовательностях и строит линейный график.

        Args:
            output_file: Имя файла для сохранения графика.

        Returns:
            Path: Путь к сохраненному графику.

        Examples:
            >>> with FastqAnalyzer("reads.fastq") as a:
            ...     a.per_base_sequence_content("content.png").exists()
            True
        """
        base_content = {'A': [], 'T': [], 'G': [], 'C': []}
        max_length = 0
        
        # Первый проход: сбор всех последовательностей и определение максимальной длины
        sequences_data = []
        for _, sequence, _ in self.get_sequences_with_quality():
            sequences_data.append(sequence)
            max_length = max(max_length, len(sequence))
        
        if max_length == 0:
            print("Нет данных для графика содержания")
            return output_path
        
        # Инициализация счетчиков для каждой позиции
        for base in base_content:
            base_content[base] = [0] * max_length
        total_bases = [0] * max_length  # Общее количество оснований на каждой позиции
        
        # Второй проход: подсчет оснований по позициям
        for sequence in sequences_data:
            for pos, base in enumerate(sequence):
                if base in base_content and pos < max_length:
                    base_content[base][pos] += 1
                    total_bases[pos] += 1
        
        # Расчет процентного содержания для каждого основания
        percentages = {}
        for base in base_content:
            percentages[base] = [100 * count / total if total > 0 else 0 
                               for count, total in zip(base_content[base], total_bases)]
        
        # Построение графика
        plt.figure(figsize=(12, 6))
        positions = list(range(max_length))
        
        # Цвета для разных оснований
        for base, color in zip(['A', 'T', 'G', 'C'], ['green', 'red', 'black', 'blue']):
            plt.plot(positions, percentages[base], label=base, color=color, linewidth=2)
        
        plt.xlabel('Position in read (bp)')
        plt.ylabel('Percentage (%)')
        plt.title('Per Base Sequence Content')
        plt.legend()
        plt.grid(True, alpha=0.3)
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        output_path = OUTPUT_DIR / output_file
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"График содержания создан: {output_path}")
        return output_path

    def sequence_length_distribution(self, output_file: str = "fastq_length_plot.png"):
        """
        Анализирует и визуализирует распределение длин последовательностей.
        
        Строит гистограмму, показывающую частоту встречаемости ридов разной длины.
        
        Args:
            output_file: Имя файла для сохранения гистограммы
            
        Returns:
            Path: Путь к сохраненному графику

        Examples:
            >>> with FastqAnalyzer("reads.fastq") as a:
            ...     a.sequence_length_distribution("lengths.png").exists()
            True
        """
        # Сбор данных о длинах последовательностей
        lengths = []
        for record in self.read():
            lengths.append(len(record.sequence))
        
        if not lengths:
            print("Нет данных для графика распределения длин")
            return output_path
            
        # Построение гистограммы
        plt.figure(figsize=(10, 6))
        plt.hist(lengths, bins=20, alpha=0.7, edgecolor='black')
        plt.xlabel('Sequence Length (bp)')
        plt.ylabel('Frequency')
        plt.title('Sequence Length Distribution')
        plt.grid(True, alpha=0.3)
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        output_path = OUTPUT_DIR / output_file
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"График распределения длин создан: {output_path}")
        return output_path


def demo_fastq_analysis(filepath: str):
    """
    Демонстрационная функция для комплексного анализа FASTQ файла.
    
    Выполняет полный анализ FASTQ файла включая:
    - Проверку читаемости файла
    - Базовую статистику
    - Создание диагностических графиков
    
    Args:
        filepath: Путь к FASTQ файлу для анализа
    
    Examples:
        >>> demo_fastq_analysis("reads.fastq")  # doctest: +SKIP
        ========== Анализ FASTQ ==========
        1. Количество последовательностей: 1000
        2. Средняя длина последовательностей: 150.00
    """
    print("========== Анализ FASTQ ==========")
    print(f"Выбранный файл: {filepath}")
    
    # Предварительная проверка файла
    with open(filepath, "r", encoding="utf-8") as f:
        print("Файл читается:")
        first_lines = [line.strip() for line in f.readlines()[:4]]
        if len(first_lines) < 4:
            print(f"❌ Файл содержит только {len(first_lines)} строк - он не является FASTQ файлом."
                  f"(Нужно минимум 4 строки на одну запись)")
            return
        if not first_lines[0].startswith("@"):
            print (f"❌ Строка 1 должна начинаться с '@', а начинается с "
                  f"{first_lines[0][:1]!r} - файл не является FASTQ файлом.")
            return
        if not first_lines[2].startswith("+"):
            print (f"❌ Строка 3 должна начинаться с '+', а начинается с "
                  f"{first_lines[2][:1]!r} - файл не является FASTQ файлом.")
            return
        
    # Создание анализатора и выполнение анализа
    analyzer = FastqAnalyzer(filepath)
    
    # Базовая статистика
    sequences = list(analyzer.read())
    print(f"1. Количество последовательностей: {analyzer.get_seq_count()}")
    print(f"2. Средняя длина последовательностей: {analyzer.get_mean_seq_length():.2f}")
    
    if analyzer.get_seq_count() == 0:
        print("Нет последовательностей для анализа!")
        return
    
    # Создание диагностических графиков
    print("3. Создание графиков...")
    quality_plot = analyzer.per_base_sequence_quality()
    content_plot = analyzer.per_base_sequence_content()
    length_plot = analyzer.sequence_length_distribution()
    
    print(f"4. График качества: {quality_plot}")
    print(f"5. График содержания: {content_plot}")
    print(f"6. Распределение длин: {length_plot}")    




    