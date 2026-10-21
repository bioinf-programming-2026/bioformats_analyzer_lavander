from abc import ABC, abstractmethod
import os
from typing import Iterator, Any


class BioParserError(Exception):
    """Базовое исключение для ошибок парсинга био-форматов."""
    pass

class EmptyFileError(BioParserError):
    """Вызывается, если переданный файл пуст."""
    pass

class FastaParsingError(BioParserError):
    """Вызывается при нарушении спецификации формата FASTA."""
    pass


class Reader(ABC):
    """Абстрактный класс ридера с безопасным управлением ресурсами."""

    def __init__(self, filename: str):
        if not isinstance(filename, str):
            raise TypeError(f"Имя файла должно быть строкой, получено: {type(filename)}")
        self.filename = filename
        self._file = None

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def open(self):
        """Открытие файла с проверками на существование и доступность."""
        if not os.path.exists(self.filename):
            raise FileNotFoundError(f"Файл '{self.filename}' не найден.")
        if os.path.isdir(self.filename):
            raise IsADirectoryError(f"Указанный путь является директорией: '{self.filename}'")
        if os.path.getsize(self.filename) == 0:
            raise EmptyFileError(f"Файл '{self.filename}' пуст.")

        if self._file is None or self._file.closed:
            self._file = open(self.filename, "r", encoding="utf-8", errors="replace")

    def close(self):
        """Безопасное закрытие файла."""
        if self._file and not self._file.closed:
            self._file.close()

    @abstractmethod
    def read(self) -> Iterator[Any]:
        """Генератор чтения записей."""
        pass

    def __iter__(self):
        return self.read()

class SequenceReader(Reader):
    """Абстрактный класс для последовательностей."""

    @abstractmethod
    def validate_sequence(self, sequence: str) -> bool:
        """Проверка последовательности на валидность символов."""
        pass

