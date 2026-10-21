class FastaRecord:
    """Представление одной FASTA-записи (заголовок + последовательность)."""

    def __init__(self, record_id, sequence, description):
        if not record_id or not isinstance(record_id, str):
            raise ValueError("ID записи не может быть пустым и должен быть строкой.")
        if not isinstance(sequence, str):
            raise TypeError("Последовательность должна быть строкой.")

        self.id = record_id.strip()
        self.description = description.strip()
        self.sequence = sequence.strip()

    def __len__(self) -> int:
        """Возвращает длину последовательности."""
        return len(self.sequence)

    def __repr__(self) -> str:
        return f"<FastaRecord id='{self.id}' len={len(self)}>"