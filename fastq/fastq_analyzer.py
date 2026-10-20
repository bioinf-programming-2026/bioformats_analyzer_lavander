"""Чтение, анализ и фильтрация FASTQ."""

from abstract import AbstractAnalyzer
from record import FastqRecord


class FastqAnalyzer(AbstractAnalyzer):
    """Анализирует четырёхстрочный FASTQ с качеством Phred+33."""

    def __init__(self, file_path):
        """Сохранить путь к файлу и создать список записей."""
        super().__init__(file_path)
        self.records = []

    def read_records(self):
        """Прочитать файл и вернуть список записей."""
        records = []

        with open(self.file_path, "r", encoding="ascii") as file:
            while True:
                header = file.readline()

                if header == "":
                    break

                sequence = file.readline()
                separator = file.readline()
                quality = file.readline()

                if sequence == "" or separator == "" or quality == "":
                    raise ValueError("Неполная FASTQ-запись.")

                header = header.rstrip("\r\n")
                sequence = sequence.rstrip("\r\n")
                separator = separator.rstrip("\r\n")
                quality = quality.rstrip("\r\n")

                if not header.startswith("@"):
                    raise ValueError("Заголовок должен начинаться с @.")

                if not separator.startswith("+"):
                    raise ValueError("Третья строка должна начинаться с +.")

                record = FastqRecord(
                    header[1:],
                    sequence,
                    quality,
                )

                records.append(record)

        self.records = records
        return self.records

    def get_statistics(self):
        """Рассчитать статистику загруженных записей."""
        if not self.records:
            return {
                "reads_count": 0,
                "total_bases": 0,
                "mean_length": None,
                "gc_percent": None,
                "mean_quality": None,
            }

        total_bases = 0
        total_gc = 0
        total_quality = 0

        for record in self.records:
            sequence = record.sequence.upper()

            total_bases += len(record)
            total_gc += sequence.count("G") + sequence.count("C")

            for symbol in record.quality:
                total_quality += ord(symbol) - 33

        return {
            "reads_count": len(self.records),
            "total_bases": total_bases,
            "mean_length": total_bases / len(self.records),
            "gc_percent": total_gc / total_bases * 100,
            "mean_quality": total_quality / total_bases,
        }

    def filter_records(self, min_quality=20, min_length=1):
        """Вернуть записи с подходящими качеством и длиной."""
        filtered = []

        for record in self.records:
            if (
                record.mean_quality() >= min_quality
                and len(record) >= min_length
            ):
                filtered.append(record)

        return filtered

    def save_records(self, records, output_path):
        """Сохранить записи в новый файл без перезаписи существующего."""
        with open(output_path, "x", encoding="ascii") as file:
            for record in records:
                file.write(record.to_fastq())