from typing import List, Dict, Optional


class VariantRecord:
    """Хранит данные одной мутации VCF и реализует интервальные проверки."""

    def __init__(
        self,
        chrom: str,
        pos: int,
        var_id: str,
        ref: str,
        alt: str,
        qual: str,
        filt: str,
        info_str: str
    ):
        self.chrom = str(chrom).strip()
        self.pos = int(pos)
        self.var_id = var_id.strip()
        self.ref = ref.strip().upper()
        # ALT может содержать несколько аллелей через запятую
        self.alt: List[str] = [a.strip().upper() for a in alt.split(",")] if alt != "." else []
       
        # QUAL может быть числом или точкой '.'
        self.qual: Optional[float] = float(qual) if qual != "." else None
        self.filter = filt.strip()

        self.info: Dict[str, Any] = {}
        if info_str and info_str != ".":
            for item in info_str.strip().split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    self.info[k] = v
                elif item:
                    self.info[item] = True

    @property
    def start(self) -> int:
        """Начало мутации (1-based)."""
        return self.pos

    @property
    def end(self) -> int:
        """
        Конец мутации.
        Вычисляется по длине референсного аллеля REF (SNV: len=1, INDEL: len > 1).
        """
        return self.pos + max(len(self.ref), 1) - 1

    def intersects(self, chrom: str, start: int, end: int) -> bool:
        """
        Проверка пересечения с интервалом [start, end] (логика bedtools intersect).
        Отрезки [A, B] и [C, D] пересекаются, если не выполняется (B < C или A > D).
        """
        if self.chrom != str(chrom):
            return False
        return not (self.end < start or self.start > end)

    def __repr__(self) -> str:
        alt_repr = ",".join(self.alt) if self.alt else "."
        return f"<VariantRecord {self.chrom}:{self.pos} {self.ref}>{alt_repr}>"