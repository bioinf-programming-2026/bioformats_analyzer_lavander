from abc import ABC, abstractmethod
import os
from typing import Iterator, Any

class BioParserError(Exception):
    """Базовая ошибка парсера."""
    pass

class EmptyFileError(BioParserError):
    """Файл пуст."""
    pass

class VcfParsingError(BioParserError):
    """Ошибка валидации или формата VCF."""
    pass

class Reader(ABC):
    """Абстрактный класс ридера файлов с контекстным менеджером."""

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
        """Безопасное открытие файла."""
        if not os.path.exists(self.filename):
            raise FileNotFoundError(f"Файл '{self.filename}' не найден.")
        if os.path.isdir(self.filename):
            raise IsADirectoryError(f"Путь является папкой: '{self.filename}'")
        if os.path.getsize(self.filename) == 0:
            raise EmptyFileError(f"Файл '{self.filename}' пуст.")

        if self._file is None or self._file.closed:
            self._file = open(self.filename, "r", encoding="utf-8", errors="replace")

    def close(self):
        """Закрытие дескриптора файла."""
        if self._file and not self._file.closed:
            self._file.close()

    @abstractmethod
    def read(self) -> Iterator[Any]:
        """Генератор записей."""
        pass

    def __iter__(self):
        return self.read()

class GenomicDataReader(Reader):
    """Абстрактный класс для позиционных геномных форматов (VCF/SAM)."""

    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """Проверка валидности геномной координаты (1-based)."""
        return isinstance(pos, int) and pos > 0