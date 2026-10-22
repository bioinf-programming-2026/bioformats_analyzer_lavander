from abc import ABC, abstractmethod
from typing import Iterator, List, Dict, Any, Optional
from record import Record, SequenceRecord
from pathlib import Path

"""Абстрактные базовые классы для чтения биологических файлов.

    Модуль задаёт общий интерфейс ридеров FASTA, FASTQ, SAM и VCF:
    открытие файла через контекстный менеджер, чтение записей,
    закрытие ресурсов и базовые операции анализа.

    Examples:
    Ридеры нельзя создавать напрямую — используются наследники:

    >>> from fasta_analyzer import FastaAnalyzer
    >>> with FastaAnalyzer("demo.fasta") as reader:
    ...     for rec in reader.read():
    ...         print(rec.id, len(rec.sequence))
    """

class Reader(ABC):
    """Базовый абстрактный класс для всех ридеров биологических файлов.

    Описывает общий жизненный цикл ридера: открытие файла через
    контекстный менеджер, чтение записей и освобождение ресурсов.

    Attributes:
        filepath: Путь к читаемому файлу.
        file: Открытый файловый объект или None.

    Examples:
        >>> class MyReader(Reader):
        ...     def read(self):
        ...         yield from self.file
        >>> with MyReader("data.txt") as r:
        ...     list(r.read())
    """

    def __init__(self, filepath: str | Path):
        """Инициализирует ридер.

        Args:
            filepath: Путь к файлу в виде строки или объекта Path.

        Examples:
            >>> r = FastaAnalyzer("seqs.fasta")
            >>> r.filepath
            PosixPath('seqs.fasta')
        """
        self.filepath = Path(filepath)
        self.file = None

    @abstractmethod
    def read(self) -> Iterator[Record]:
        """Возвращает итератор по записям файла.

        Yields:
            Record: Очередная запись из файла.

        Raises:
            NotImplementedError: Если метод не реализован в наследнике.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     first = next(r.read())
            ...     first.chrom
            'chr1'
        """
        pass

    def __enter__(self):
        """Открывает файл при входе в блок `with`.

        Returns:
            Reader: Текущий экземпляр ридера.

        Examples:
            >>> with FastaAnalyzer("seqs.fasta") as r:
            ...     r.file.closed
            False
        """
        self.file = open(self.filepath, "r", encoding="utf-8")
        return self
    
    def close(self):
        """Закрывает файл, если он открыт, и сбрасывает ссылку на него.

        Examples:
            >>> r = FastaAnalyzer("seqs.fasta")
            >>> r.__enter__()
            >>> r.close()
            >>> r.file is None
            True
        """
        if self.file and not self.file.closed:
            self.file.close()
            self.file = None

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Закрывает файл при выходе из блока `with`.

        Args:
            exc_type: Тип исключения, если оно возникло.
            exc_val: Значение исключения.
            exc_tb: Traceback исключения.

        Returns:
            None: Исключения не подавляются.

        Examples:
            >>> with FastaAnalyzer("seqs.fasta") as r:
            ...     pass
            >>> r.file is None
            True
            
        Действия:
            - Вызывает родительский __exit__ для закрытия файла
            - Сбрасывает все внутренние состояния ридера
        """
        self.close()


class SequenceReader(Reader):
    """Абстрактный ридер для форматов последовательностей FASTA и FASTQ.

    Добавляет операции, характерные для последовательностей:
    создание записи, валидацию, подсчёт количества и средней длины.

    Examples:
        >>> with FastqAnalyzer("reads.fastq") as reader:
        ...     list(reader.read())
        ...     reader.get_seq_count()
        1000
    """

    def __init__(self, filepath: str | Path):
        """Инициализирует ридер последовательностей.

        Args:
            filepath: Путь к FASTA/FASTQ-файлу.

        Examples:
            >>> FastaAnalyzer("seqs.fasta").filepath
            PosixPath('seqs.fasta')
        """
        super().__init__(filepath)

    @abstractmethod
    def read(self) -> Iterator[SequenceRecord]:
        """Возвращает итератор по объектам SequenceRecord.

        Yields:
            SequenceRecord: Очередная биологическая последовательность.

        Examples:
            >>> with FastaAnalyzer("seqs.fasta") as r:
            ...     next(r.read()).id
            'seq1'
        """
        pass

    @abstractmethod
    def get_sequence(self, id: str, seq: str) -> SequenceRecord:
        """Создаёт объект SequenceRecord из идентификатора и последовательности.

        Args:
            id: Идентификатор последовательности.
            seq: Нуклеотидная последовательность.

        Returns:
            SequenceRecord: Созданная запись последовательности.

        Examples:
            >>> FastaAnalyzer("x.fasta").get_sequence("seq1", "ATCG").sequence
            'ATCG'
        """

        return SequenceRecord(id=id, sequence=seq)

    @abstractmethod
    def validate_sequence(self, seq: str) -> bool:
        """Проверяет корректность последовательности.

        Args:
            seq: Проверяемая последовательность.

        Returns:
            bool: True, если последовательность непустая и состоит из букв.

        Examples:
            >>> FastaAnalyzer("x.fasta").validate_sequence("ATCG")
            True
            >>> FastaAnalyzer("x.fasta").validate_sequence("")
            False
        """
        return bool(seq) and seq.isalpha()

    @abstractmethod
    def get_seq_count(self) -> int:
        """Возвращает количество прочитанных последовательностей.

        Returns:
            int: Количество последовательностей.

        Examples:
            >>> with FastaAnalyzer("seqs.fasta") as r:
            ...     list(r.read())
            ...     r.get_seq_count()
            2
        """
        pass
    
    @abstractmethod
    def get_mean_seq_length(self) -> float:
        """Возвращает среднюю длину последовательностей.

        Returns:
            float: Средняя длина; 0.0, если последовательностей нет.

        Examples:
            >>> with FastaAnalyzer("seqs.fasta") as r:
            ...     list(r.read())
            ...     r.get_mean_seq_length()
            4.0
        """
        pass


class GenomicReader(Reader):
    """Абстрактный ридер для геномных форматов SAM и VCF.

    Инкапсулирует общую логику работы с заголовком файла,
    хромосомами, координатами и фильтрацией записей.

    Attributes:
        _header_parsed: Флаг, указывающий, был ли разобран заголовок.
        _chromosomes: Внутренний список хромосом.

    Examples:
        >>> with VcfReader("variants.vcf") as reader:
        ...     reader.get_chromosomes()
        ['chr1', 'chr2']
    """ 

    def __init__(self, filepath: str | Path):
        """Инициализирует геномный ридер.

        Args:
            filepath: Путь к SAM/VCF-файлу.

        Examples:
            >>> SamReader("aln.sam").filepath
            PosixPath('aln.sam')
        """
        super().__init__(filepath)
        self._header_parsed = False
        self._chromosomes = []

    def __enter__(self):
        """Открывает файл и разбирает заголовок при входе в блок `with`.

        Returns:
            GenomicReader: Текущий экземпляр ридера.

        Raises:
            RuntimeError: Если файл не удалось открыть или разобрать заголовок.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     r._header_parsed
            True
        """
        try:
            self.file = open(self.filepath, "r", encoding="utf-8")
            self._parse_header()
            return self
        except Exception as e:
            # если ошибка — аккуратно закрываем файл, чтобы не остался висеть
            if self.file and not self.file.closed:
                self.file.close()
            raise RuntimeError(f"Ошибка при открытии или парсинге файла {self.filepath}") from e

    @abstractmethod
    def _parse_header(self):
        """Разбирает заголовок геномного файла.

        Raises:
            NotImplementedError: Если метод не реализован в наследнике.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     r._parse_header()
            ...     r.samples
            ['sample1', 'sample2']
        """
        pass

    @abstractmethod
    def read(self) -> Iterator[Record]:
        """Возвращает итератор по геномным записям, пропуская заголовок.

        Yields:
            Record: Очередная геномная запись.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     next(r.read()).chrom
            'chr1'
        """
        pass

    @abstractmethod
    def get_chromosomes(self) -> List[str]:
        """Возвращает список хромосом, найденных в файле.

        Returns:
            List[str]: Отсортированный список названий хромосом.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     r.get_chromosomes()[:2]
            ['chr1', 'chr2']
        """
        pass

    @abstractmethod
    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """Проверяет корректность геномной координаты.

        Args:
            chrom: Название хромосомы.
            pos: Позиция на хромосоме.

        Returns:
            bool: True, если координата корректна.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     r.validate_coordinate("chr1", 100)
            True
            >>> with SamReader("aln.sam") as r:
            ...     r.validate_coordinate("chrX", 100)
            False
        """
        pass

    @abstractmethod
    def get_records_in_region(self, chrom: str, start: int, end: int) -> Iterator[Record]:
        """Возвращает записи, попадающие в указанный регион.

        Args:
            chrom: Название хромосомы.
            start: Начало региона.
            end: Конец региона.

        Yields:
            Record: Записи, пересекающие регион.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     list(r.get_records_in_region("chr1", 1000, 2000))
            [<Record id=chr1:1500>]
        """
        pass

    @abstractmethod
    def filter_records(self, **filters) -> Iterator[Record]:
        """Фильтрует записи по заданным критериям.

        Args:
            **filters: Критерии фильтрации, зависящие от реализации.

        Yields:
            Record: Записи, удовлетворяющие фильтрам.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     list(r.filter_records(chrom="chr1", min_mapq=30))[:1]
            [<Record id=read1>]
        """
        pass

    def get_chromosome_length(self, chrom: str) -> Optional[int]:
        """Возвращает длину хромосомы, если она известна из заголовка.

        Args:
            chrom: Название хромосомы.

        Returns:
            Optional[int]: Длина хромосомы или None, если данных нет.

        Examples:
            >>> with SamReader("aln.sam") as r:
            ...     r.get_chromosome_length("chr1")
            248956422
        """
        # Базовая реализация, может быть переопределена в дочерних классах
        return None

    def get_statistics(self) -> Dict[str, Any]:
        """Возвращает базовую статистику по файлу.

        Returns:
            Dict[str, Any]: Путь, размер файла, список хромосом и их количество.

        Examples:
            >>> with VcfReader("variants.vcf") as r:
            ...     r.get_statistics()["chromosome_count"]
            2
        """
        return {
            "file_path": str(self.filepath),
            "file_size": self.filepath.stat().st_size if self.filepath.exists() else 0,
            "chromosomes": self.get_chromosomes(),
            "chromosome_count": len(self.get_chromosomes()),
        }

    def close(self):
        """Закрывает файл и сбрасывает флаг разбора заголовка.

        Examples:
            >>> r = SamReader("aln.sam")
            >>> r.__enter__()
            >>> r.close()
            >>> r._header_parsed
            False
        """
        super().close()
        self._header_parsed = False