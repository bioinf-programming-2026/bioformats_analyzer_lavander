"""Класс одной FASTQ-записи."""

from abstract import AbstractRecord


class FastqRecord(AbstractRecord):
    """Хранит заголовок, последовательность и строку качества."""

    def __init__(self, header, sequence, quality):
        """Создать запись и проверить её данные."""
        if not header.strip():
            raise ValueError("Пустой заголовок.")

        if not sequence:
            raise ValueError("Пустая последовательность.")

        if len(sequence) != len(quality):
            raise ValueError(
                "Длины последовательности и строки качества не совпадают."
            )

        for symbol in quality:
            if not 33 <= ord(symbol) <= 126:
                raise ValueError("Недопустимый символ качества.")

        self.header = header
        self.sequence = sequence
        self.quality = quality

    def __len__(self):
        """Вернуть длину последовательности."""
        return len(self.sequence)

    def mean_quality(self):
        """Вернуть среднее качество Phred+33."""
        scores = [ord(symbol) - 33 for symbol in self.quality]
        return sum(scores) / len(self)

    def gc_content(self):
        """Вернуть процент G и C от полной длины последовательности."""
        sequence = self.sequence.upper()
        gc_count = sequence.count("G") + sequence.count("C")
        return gc_count / len(self) * 100

    def to_fastq(self):
        """Вернуть запись в четырёхстрочном формате FASTQ."""
        return (
            f"@{self.header}\n"
            f"{self.sequence}\n"
            f"+\n"
            f"{self.quality}\n"
        )