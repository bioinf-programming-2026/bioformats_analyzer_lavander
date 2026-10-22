from pathlib import Path
from typing import Iterator
from abstract import SequenceReader
from record import SequenceRecord

"""Анализатор FASTA-файлов.

Читает последовательности, считает их количество и среднюю длину,
поддерживает работу через контекстный менеджер.

Examples:
    >>> with FastaAnalyzer("seqs.fasta") as a:
    ...     list(a.read())
    ...     (a.get_seq_count(), a.get_mean_seq_length())
    (2, 4.0)
"""

class FastaAnalyzer(SequenceReader):
    """Анализатор FASTA-файлов.

    Анализатор FASTA файлов для чтения и анализа биологических последовательностей.
    
    Реализует:
    1. Подсчет количества последовательностей в файле
    2. Расчет средней длины последовательностей
    3. Итеративное чтение записей с поддержкой контекстного менеджера
    
    Attributes:
        filepath (Path): Путь к анализируемому FASTA файлу
        file (file object): Открытый файловый объект
        _seq_count (int): Внутренний счетчик последовательностей
        _total_length (int): Суммарная длина всех последовательностей

    Examples:
        >>> with FastaAnalyzer("seqs.fasta") as a:
        ...     for rec in a.read():
        ...         print(rec.id, rec.sequence)
        seq1 ATCG
        seq2 GGCC
    """

    def __init__(self, filepath: str | Path):
        """Инициализирует анализатор FASTA.

        Args:
            filepath: Путь к FASTA-файлу.

        Examples:
            >>> FastaAnalyzer("seqs.fasta")._seq_count
            0
        """
        super().__init__(filepath)
        self._seq_count = 0  # Счетчик последовательностей
        self._total_length = 0  # Накопитель общей длины последовательностей

    def read(self) -> Iterator[SequenceRecord]:
        """
        Основной метод чтения FASTA файла с подсчетом статистики.
        
        Генератор, который последовательно возвращает записи из FASTA файла.
        Одновременно обновляет внутреннюю статистику для последующих запросов.
        
        Логика обработки:
        - Пропускает пустые строки
        - Строки, начинающиеся с '>' - заголовки последовательностей
        - Последующие строки до следующего заголовка - данные последовательности
        
        Yields:
            SequenceRecord: Объект записи с идентификатором и последовательностью
        
        Raises:
            RuntimeError: Если файл не был открыт через контекстный менеджер.

        Обновляет:
            _seq_count: увеличивает на 1 для каждой найденной последовательности
            _total_length: добавляет длину каждой последовательности

        Examples:
            >>> with FastaAnalyzer("seqs.fasta") as a:
            ...     [rec.id for rec in a.read()]
            ['seq1', 'seq2']
        """
        # Автоматическое открытие файла если он не был открыт через контекстный менеджер
        if self.file is None:
            raise RuntimeError(
                "FastaAnalyzer должен использоваться через контекстный менеджер: "
                "`with FastaAnalyzer(path) as a:`"
            )

        current_header = None  # Текущий обрабатываемый заголовок
        current_sequence = []  # Накопитель для строк последовательности
        
        for line in self.file:
            line = line.strip()
            
            if not line:  # Пропуск пустых строк
                continue
                
            if line.startswith(">"):  # Обнаружен новый заголовок
                # Если уже обрабатывали предыдущую последовательность - отдаем ее
                if current_header is not None:
                    sequence = "".join(current_sequence)

                    
                    # ОБНОВЛЕНИЕ СТАТИСТИКИ
                    self._seq_count += 1
                    self._total_length += len(sequence)
                    
                    yield SequenceRecord(id=current_header, sequence=sequence)
                
                # Начинаем новую последовательность
                current_header = line[1:].strip()  # Убираем '>' и лишние пробелы
                current_sequence = []
            else:
                # Накопление данных последовательности
                current_sequence.append(line)
        
        # Обработка последней последовательности в файле
        if current_header is not None:
            sequence = "".join(current_sequence)
            
            self._seq_count += 1
            self._total_length += len(sequence)
            
            yield SequenceRecord(id=current_header, sequence=sequence)

    def get_sequence(self, id: str, seq: str) -> SequenceRecord:
        """Создаёт объект SequenceRecord.

        Args:
            id: Идентификатор последовательности.
            seq: Нуклеотидная последовательность.

        Returns:
            SequenceRecord: Созданная запись.

        Examples:
            >>> FastaAnalyzer("x.fasta").get_sequence("s1", "AC").sequence
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
            >>> FastaAnalyzer("x.fasta").validate_sequence("ACGT")
            True
        """
        return bool(seq) and seq.isalpha()

    def get_seq_count(self) -> int:
        """
        Возвращает количество прочитанных последовательностей
        
        Returns:
            int: Общее количество последовательностей в файле
            
        Note:
            Для корректного подсчета требуется предварительный вызов read()
        
        Examples:
            >>> with FastaAnalyzer("seqs.fasta") as a:
            ...     list(a.read())
            ...     a.get_seq_count()
            2
        """
        return self._seq_count

    def get_mean_seq_length(self) -> float:
        """
        Возвращает среднюю длину последовательностей.

        Returns:
            float: Средняя длина последовательности

        Examples:
            >>> with FastaAnalyzer("seqs.fasta") as a:
            ...     list(a.read())
            ...     a.get_mean_seq_length()
            4.0
            
        Note:
            Возвращает 0.0 если файл не содержит последовательностей
        """
        if self._seq_count == 0:
            return 0.0
        return self._total_length / self._seq_count


def demo_fasta_analysis(filepath: str):
    """
    Демонстрационная функция анализа FASTA файла.

    Показывает использование основных возможностей класса FastaAnalyzer:
    1. Чтение всех последовательностей
    2. Подсчет общего количества
    3. Расчет средней длины

    Args:
        filepath: Путь к FASTA-файлу для анализа.

    Examples:
        >>> demo_fasta_analysis("seqs.fasta")  # doctest: +SKIP
        ========== Анализ FASTA ==========
        1. Количество последовательностей: 2
        2. Средняя длина последовательностей: 4.00
    """
    print("========== Анализ FASTA ==========")
    print(f"Выбранный файл: {filepath}")
    
    # Использование контекстного менеджера для гарантированного закрытия файла
    with FastaAnalyzer(filepath) as analyzer:
        # Чтение всех последовательностей (одновременно обновляет статистику)
        sequences = list(analyzer.read())
        
        # ВЫВОД РЕЗУЛЬТАТОВ СОГЛАСНО ТРЕБОВАНИЯМ ЗАДАНИЯ:
        print(f"1. Количество последовательностей: {analyzer.get_seq_count()}")
        print(f"2. Средняя длина последовательностей: {analyzer.get_mean_seq_length():.2f}")