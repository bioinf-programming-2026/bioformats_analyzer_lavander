import os
import sys
import tempfile

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fastq_analyzer import FastqAnalyzer, demo_fastq_analysis
from fasta_analyzer import FastaAnalyzer, demo_fasta_analysis
from sam_analyzer import SamReader
from vcf_analyzer import VcfReader

"""Демонстрация и тесты анализаторов FASTA, FASTQ, SAM и VCF.

Модуль запускает небольшие проверки и демонстрирует использование
соответствующих ридеров на пользовательских файлах.

Examples:
    Запуск из командной строки:

    $ python demo.py
    Введите формат файла: fasta
    Введите путь к файлу: seqs.fasta
"""

def print_header(title: str):
    """Оформление заголовка блока.
        
    Args:
        title: Заголовок блока.
        
    Examples:
        >>> print_header("FASTA: ТЕСТЫ")
        <BLANKLINE>
        ==========================================================================
        >>> FASTA: ТЕСТЫ
        ==========================================================================
    """
    print("\n" + "=" * 70)
    print(f" >>> {title}")
    print("=" * 70)



# ====== FASTA ======

def tests_fasta_analysis():
    """Проверяет базовый анализ FASTA на временном файле.

    Проверяет количество последовательностей и среднюю длину.
    Временный файл удаляется после завершения теста.

    Raises:
        AssertionError: Если результаты анализа не совпадают с ожидаемыми.

    Examples:
        >>> tests_fasta_analysis()
        Запуск тестов FASTA...
        Тесты FASTA пройдены!
    """
    print("\nЗапуск тестов FASTA...")
    fasta_content = """>seq1
ATCG
>seq2
GGCC
"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
        f.write(fasta_content)
        temp_file = f.name
    
    try:
        with FastaAnalyzer(temp_file) as analyzer:
            sequences = list(analyzer.read())
            assert analyzer.get_seq_count() == 2
            assert analyzer.get_mean_seq_length() == 4.0
            print("Тесты FASTA пройдены!", end="")
    finally:
        os.unlink(temp_file)


# ====== FASTQ ======

def test_fastq_analyzer():
    """Проверяет подсчёт последовательностей, среднюю длину и Phred-качество FASTQ.

    Создаёт временные FASTQ-файлы и проверяет основные методы FastqAnalyzer.

    Raises:
        AssertionError: Если результаты анализа не совпадают с ожидаемыми.

    Examples:
        >>> test_fastq_analyzer()  # doctest: +SKIP
    """
    
    def create_test_fastq(content: str) -> str:
        """Создаёт временный FASTQ-файл с заданным содержимым.

        Args:
            content: Содержимое FASTQ-файла.

        Returns:
            str: Путь к созданному временному файлу.

        Examples:
            >>> path = create_test_fastq("@r\\nAT\\n+\\nII\\n")
            >>> open(path).readline().strip()
            '@r'
        """
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fastq') as f:
            f.write(content)
        return f.name
    

    # Тест подсчета последовательностей в FASTQ.
    fastq_content = """@read1
ATCG
+
IIII
@read2
GGCC
+
IIII
"""
    test_file = create_test_fastq(fastq_content)
        
    try:
        with FastqAnalyzer(test_file) as analyzer:
            sequences = list(analyzer.read())
            assert analyzer.get_seq_count() == 2
    finally:
        os.unlink(test_file)
    
    # Тест расчета средней длины в FASTQ
    fastq_content = """@read1
ATCG
+
IIII
@read2
GGCCAA
+
IIIIII
"""
    test_file = create_test_fastq(fastq_content)
        
    try:
        with FastqAnalyzer(test_file) as analyzer:
            sequences = list(analyzer.read())
            assert analyzer.get_mean_seq_length() == (4 + 6) / 2
    finally:
        os.unlink(test_file)
    
    # Тест конвертации Phred качества
    analyzer = FastqAnalyzer("dummy.fastq")
    # Символ 'I' имеет ASCII код 73, при offset 33: 73-33=40
    assert analyzer.phred_to_quality('I') == 40


def tests_fastq_analysis():
    """Запускает тесты FASTQ-анализатора.

    Создаёт временный FASTQ-файл, вызывает test_fastq_analyzer()
    и удаляет файл после завершения.

    Examples:
        >>> tests_fastq_analysis()
        Запуск тестов FASTQ...
        Тесты FASTQ пройдены!
    """
    print("\nЗапуск тестов FASTQ...")

    # Создаем тестовый FASTQ файл
    fastq_content = """@test_read
ATCG
+
IIII
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write(fastq_content)
        temp_file = f.name

    
    try:
        # Проверяем, что функция выполняется без ошибок
        test_fastq_analyzer()
        print("Тесты FASTQ пройдены!", end="")
    finally:
        os.unlink(temp_file)


# ====== SAM ======

def demo_sam(file_path, start=1, end=10000):
    """Демонстрирует работу SamReader на SAM-файле.

    Выводит заголовок, количество выравниваний, хромосомы,
    статистику по хромосомам, записи в регионе и первую запись.

    Args:
        file_path: Путь к SAM-файлу.
        start: Начало геномного отрезка для bedtools-intersect.
        end: Конец геномного отрезка.

    Examples:
        >>> demo_sam("aln.sam")  # doctest: +SKIP
        Тестирование SAM ридера...
        ========== Анализ SAM ==========
    """
    if not os.path.exists(file_path):
        print(f"Файл {file_path} не найден!")
        print("Убедитесь, что файл находится в той же директории")
        return
    
    try:
        print("\nТестирование SAM ридера...")
        print("\n============================= Анализ SAM =============================\n")
        print(f"Выбранный файл: {file_path}")
        print("=" * 70)


        with SamReader(file_path) as reader:
            # --------------------------------------------------------------
            # 1. Заголовок и группы метаданных (@SQ, @RG, @PG)
            # --------------------------------------------------------------
            print_header("1. ЗАГОЛОВОК И ГРУППЫ МЕТАДАННЫХ")

            header = reader.get_header()                # группы (уже есть)
            fields = reader.get_header_fields()         # разбор @SQ/@RG/@PG

            print(f"Всего групп: {len(header)}")
            print(f"Группы в файле: {list(header.keys())}")
            print()

            # @SQ — референсы
            sq = fields.get("@SQ", {})
            if sq:
                print(f" • Группа '@SQ' (референсов: {len(sq)}):")
                for name, info in list(sq.items())[:2]:
                    print(f"     {name}: length={info.get('LN', '?')}")
                if len(sq) > 2:
                    print(f"     ... и ещё {len(sq) - 2} референсов")
                print()

            # @RG — read groups
            rg = fields.get("@RG", {})
            if rg:
                print(f" • Группа '@RG' (read groups: {len(rg)}):")
                for rg_id in list(rg.keys())[:2]:        
                    print(f"     ID={rg_id}")                 # read groups по ID
                if len(rg) > 2:
                    print(f"     ... и ещё {len(rg) - 2}")
                print()

            # @PG — программы
            pg = fields.get("@PG", {})
            if pg:
                print(f" • Группа '@PG' (программ: {len(pg)}):")
                for pg_id in list(pg.keys())[:2]:
                    print(f"     ID={pg_id}")
                if len(pg) > 2:
                    print(f"     ... и ещё {len(pg) - 2}")
                print()

            # --------------------------------------------------------------
            # 2. Подсчёт количества выравниваний (генератор)
            # --------------------------------------------------------------
            print_header("2. ПОДСЧЁТ КОЛИЧЕСТВА ВЫРАВНИВАНИЙ (ГЕНЕРАТОР)")

            count = reader.count_alignments()
            print(f"Общее число выравниваний: {count}")

            # --------------------------------------------------------------
            # 3. Статистика «количество выравниваний — регион» (pandas)
            # --------------------------------------------------------------
            print_header("3. СТАТИСТИКА «КОЛИЧЕСТВО ВЫРАВНИВАНИЙ — РЕГИОН» (PANDAS)")

            stats_df = reader.stats_by_chromosome()
            if not stats_df.empty:
                for _, row in stats_df.iterrows():
                    print(f"  {row.iloc[0]}: {row.iloc[1]} выравниваний")
            else:
                print("  Данные по регионам отсутствуют")

            # --------------------------------------------------------------
            # 4. Пересечение с геномным отрезком (bedtools intersect)
            # --------------------------------------------------------------
            chroms = reader.get_chromosomes()
            if chroms:
                target_chrom = chroms[0]
                start_pos, end_pos = start, end

                print_header(
                    f"4. ПЕРЕСЕЧЕНИЕ: {target_chrom}:{start_pos}-{end_pos} "
                    f"(BEDTOOLS INTERSECT)"
                )

                matched = list(reader.get_records_in_region(
                    target_chrom, start_pos, end_pos))
                count_in_matched = len(matched)
                print(f"Найдено выравниваний в интервале: {count_in_matched}")
                print()

                if count_in_matched > 0:
                    print("Примеры найденных выравниваний:")
                    for i, rec in enumerate(matched[:5], 1):
                        print(f"  [{i}] {rec.chrom}:{rec.start}  mapq={rec.mapq}")
                    if count_in_matched > 5:
                        print(f"  ... и ещё {len(matched) - 5} выравниваний")

    except (FileNotFoundError, IsADirectoryError, ValueError) as e:
        print(f"\n[ОШИБКА ФОРМАТА SAM]: {e}")
    except Exception as e:
        print(f"\n[НЕПРЕДВИДЕННАЯ ОШИБКА]: {e}")
        import traceback
        traceback.print_exc()


# ====== VCF ======

def demo_vcf(vcf_file, start=1, end=10000):
    """Демонстрирует работу VcfReader на VCF-файле.

    Показывает статистику файла, группы метаданных заголовка, количество
    вариантов, первые варианты, распределение по хромосомам, типы
    вариантов, пересечение с регионом и потоковое чтение.

    Args:
        vcf_file: Путь к VCF-файлу.
        start: Начало геномного отрезка для bedtools-intersect.
        end: Конец геномного отрезка.

    Examples:
        >>> demo_vcf("variants.vcf")  # doctest: +SKIP
        ================ Анализ VCF ==================
    """
    if not os.path.exists(vcf_file):
        print(f"Файл {vcf_file} не найден!")
        print("Убедитесь, что файл находится в той же директории")
        return

    try:
        print("\n============================= Анализ VCF =============================\n")
        print(f"Выбранный файл: {vcf_file}")
        print("=" * 70)

        with VcfReader(vcf_file) as reader:
            # ---------------------------------------------------------------------
            # 1. Заголовок и ВСЕ группы метаданных (стандартные + пользовательские)
            # ---------------------------------------------------------------------
            print_header("1. ЗАГОЛОВОК И ГРУППЫ МЕТАДАННЫХ")

            header = reader.get_header()

            print(f"   Всего групп: {len(header)}")
            print(f"   Всего строк заголовка: {len(reader.headers)}")
            print(f"   Группы в файле: {list(header.keys())}")
            print()

            # Стандартные группы показываем первыми (если они есть в файле),
            # затем — все остальные (включая кастомные) в порядке появления
            STANDARD_GROUPS = [
            "##fileformat", "##fileDate", "##source", "##reference",
            "##INFO", "##FILTER", "##FORMAT", "##contig",
            "##ALT", "##SAMPLE", "##PEDIGREE",
            ]
            ordered_keys = [k for k in STANDARD_GROUPS if k in header]
            ordered_keys += [k for k in header if k not in STANDARD_GROUPS]

            for group_key in ordered_keys:
                group_lines = header[group_key]
                print(f"   • Группа '{group_key}' (записей: {len(group_lines)}):")
                for line in group_lines[:3]:
                    print(f"       {line}")
                if len(group_lines) > 3:
                    print(f"       ... и ещё {len(group_lines) - 3} строк")
                print()

            # ---------------------------------------------------------------------
            # 2. Подсчёт количества вариантов (ТЗ 14.7.2.4)
            # ---------------------------------------------------------------------
            print_header("2. ПОДСЧЁТ КОЛИЧЕСТВА ВАРИАНТОВ (ГЕНЕРАТОР)")

            total = reader.get_variants_count()
            print(f"   Общее число вариантов: {total}")

            # ---------------------------------------------------------------------
            # 3. Распределение по хромосомам
            # ---------------------------------------------------------------------
            print_header("3. СТАТИСТИКА «КОЛИЧЕСТВО ВАРИАНТОВ - РЕГИОН» (PANDAS)")

            region_stats = reader.get_region_stats()
            if not region_stats.empty:
                for _, row in region_stats.iterrows():
                    print(f"   {row['chromosome']}: "
                          f"{row['variant_count']:,} вариантов")
            else:
                print("   Данные по регионам отсутствуют")

            # ---------------------------------------------------------------------
            # 4. Пересечение с регионом (bedtools intersect)
            # ---------------------------------------------------------------------
            if not region_stats.empty:
                target_chrom = region_stats["chromosome"].iloc[0]
                start_pos, end_pos = start, end

                print_header(
                    f"4. ПЕРЕСЕЧЕНИЕ: {target_chrom}:{start_pos}-{end_pos} "
                    f"(BEDTOOLS INTERSECT)"
                )

                matched = list(reader.get_records_in_region(
                    chrom=target_chrom, start=start_pos, end=end_pos))
                count_in_matched = len(matched)
                print(f"   Найдено вариантов в интервале: {count_in_matched}")
                print()
                if count_in_matched > 0:
                    print("   Примеры найденных мутаций в отрезке:")
                    for i, var in enumerate(matched[:5], 1):
                        alt_str = ",".join(var.alt) if var.alt else "."
                        print(f"   [{i}] {var.chrom}:{var.pos} {var.ref}>{alt_str}")
                        print(f"       ID: {var.var_id} | "
                              f"QUAL: {var.qual} | FILTER: {var.filter}")
                        print(f"       start={var.start}, end={var.end}")
                        print(f"       INFO: {format_info_preview(var.info)}")
                        print("   " + "-" * 50)

    except (FileNotFoundError, IsADirectoryError, ValueError) as e:
        print(f"\n[ОШИБКА ФОРМАТА VCF]: {e}")
    except Exception as e:
        print(f"\n[НЕПРЕДВИДЕННАЯ ОШИБКА]: {e}")
        import traceback
        traceback.print_exc()

def format_info_preview(info_dict):
    """Форматирует словарь INFO для краткого вывода.

    Args:
        info_dict: Словарь аннотаций VCF.

    Returns:
        str: Строка с первыми тремя полями INFO или сообщение об отсутствии данных.

    Examples:
        >>> format_info_preview({"DP": "30", "AF": "0.5"})
        'DP=30; AF=0.5'
        >>> format_info_preview({})
        'нет информации'
    """
    if not info_dict:
        return "нет информации"
    
    preview_items = []
    for key, value in list(info_dict.items())[:3]:  # Показываем первые 3 поля
        if value is True:
            preview_items.append(key)
        else:
            preview_items.append(f"{key}={value}")
    
    preview = "; ".join(preview_items)
    if len(info_dict) > 3:
        preview += f" ... (еще {len(info_dict) - 3} полей)"
    
    return preview

if __name__ == "__main__":
    supported_formats = ("fasta", "fastq", "sam", "vcf")
    print(f"Поддерживаемые форматы: {', '.join(supported_formats)}")

    fmt = input("Введите формат файла: ").strip().lower()

    if fmt not in supported_formats:
        print(f"[ОШИБКА] Неизвестный формат: '{fmt}'. "
              f"Допустимые: {', '.join(supported_formats)}")
        sys.exit(1)

    file_path = input("Введите путь к файлу: ").strip().strip('"').strip("'")

    if not os.path.exists(file_path):
        print(f"[ОШИБКА] Файл не найден: {file_path}")
        sys.exit(1)
    if os.path.isdir(file_path):
        print(f"[ОШИБКА] Указан путь к каталогу, а не файлу: {file_path}")
        sys.exit(1)

    if fmt == "fasta":
        tests_fasta_analysis()
        print()
        demo_fasta_analysis(file_path)
    elif fmt == "fastq":
        tests_fastq_analysis()
        print()
        demo_fastq_analysis(file_path)
    elif fmt == "sam":
        demo_sam(file_path)
    elif fmt == "vcf":
        demo_vcf(file_path)