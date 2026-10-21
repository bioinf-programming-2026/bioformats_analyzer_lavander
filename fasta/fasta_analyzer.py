from typing import Iterator
from abstract import SequenceReader, FastaParsingError
from record import FastaRecord


class FastaAnalyzer(SequenceReader):
    """Ридер и анализатор файлов FASTA."""

    def validate_sequence(self, sequence: str) -> bool:
        """
        Проверяет допустимость символов (ДНК, РНК, белки по IUPAC и стоп-кодоны).
        Возвращает False, если есть посторонние символы (цифры, знаки препинания).
        """
        if not sequence:
            return False
        valid_chars = set("ACGTURYKMSWBDHVNacgturykmswbdhvn*-")
        return set(sequence).issubset(valid_chars)

    def read(self) -> Iterator[FastaRecord]:
        """
        Ленивый генератор с валидацией формата FASTA на лету.
        Отлавливает:
        - отсутствие символа '>' в начале файла
        - пустые ID (> без имени)
        - пустые последовательности между заголовками
        """
        self.open()
        self._file.seek(0)

        current_id = None
        current_desc = ""
        current_seq_parts = []
        line_num = 0
        first_entry_found = False

        for line in self._file:
            line_num += 1
            line = line.strip()
           
            if not line:
                continue

            if line.startswith(">"):
                first_entry_found = True

                # Сохраняем предыдущую запись
                if current_id is not None:
                    full_seq = "".join(current_seq_parts)
                    if not full_seq:
                        raise FastaParsingError(
                            f"Ошибка в строке {line_num}: Запись '>{current_id}' не содержит последовательности!"
                        )
                    yield FastaRecord(current_id, full_seq, current_desc)

                # Разбираем заголовок
                header_content = line[1:].strip()
                if not header_content:
                    raise FastaParsingError(f"Ошибка в строке {line_num}: Пустой заголовок '>' без ID!")

                parts = header_content.split(maxsplit=1)
                current_id = parts[0]
                current_desc = parts[1] if len(parts) > 1 else ""
                current_seq_parts = []
            else:
                if not first_entry_found:
                    raise FastaParsingError(
                        f"Ошибка в строке {line_num}: Файл не начинается с символа '>'. Это не FASTA!"
                    )
                current_seq_parts.append(line)

        if current_id is not None:
            full_seq = "".join(current_seq_parts)
            if not full_seq:
                raise FastaParsingError(
                    f"Ошибка в конце файла: Запись '>{current_id}' не содержит последовательности!"
                )
            yield FastaRecord(current_id, full_seq, current_desc)

    def get_sequences_count(self) -> int:
        """Подсчет количества последовательностей (требование 14.7.2.1)."""
        count = 0
        for _ in self.read():
            count += 1
        return count

    def get_average_length(self) -> float:
        """Расчет средней длины последовательностей (требование 14.7.2.1)."""
        total_len = 0
        count = 0
        for record in self.read():
            total_len += len(record)
            count += 1
        return total_len / count if count > 0 else 0.0