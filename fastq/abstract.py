"""Базовые классы записи и анализатора."""


class AbstractRecord:
    """Задаёт общий метод преобразования записи в текст."""

    def to_fastq(self):
        """Вернуть запись в формате FASTQ."""
        raise NotImplementedError


class AbstractAnalyzer:
    """Задаёт общий интерфейс анализатора."""

    def __init__(self, file_path):
        """Сохранить путь к исходному файлу."""
        self.file_path = file_path

    def read_records(self):
        """Прочитать записи из файла."""
        raise NotImplementedError

    def get_statistics(self):
        """Рассчитать статистику файла."""
        raise NotImplementedError