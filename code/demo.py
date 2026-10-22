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
        Тест пройден!
    """

    print("Запуск тестов FASTA...")
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
            print("Тест пройден!")
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
    print("Запуск тестов FASTQ...")

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
        print("Тесты FASTQ пройдены!")
    finally:
        os.unlink(temp_file)


# ====== SAM ======

def demo_sam(file_path):
    """Демонстрирует работу SamReader на SAM-файле.

    Выводит заголовок, количество выравниваний, хромосомы,
    статистику по хромосомам, записи в регионе и первую запись.

    Args:
        file_path: Путь к SAM-файлу.

    Examples:
        >>> demo_sam("aln.sam")  # doctest: +SKIP
        Тестирование SAM ридера...
        ========== Анализ SAM ==========
    """
    print("Тестирование SAM ридера...")
    print("========== Анализ SAM ==========")
    print(f"Выбранный файл: {file_path}")
    try:
        with SamReader(file_path) as reader:
            # Базовые методы
            header = reader.get_header()
            count = reader.count_alignments()
            chroms = reader.get_chromosomes()
            
            print(f"1. Заголовок: {len(header)} строк")

            header_groups = reader.get_header()                # группы (уже есть)
            fields = reader.get_header_fields()                # разбор @SQ/@RG/@PG

            print(f"Группы заголовка: {list(header_groups.keys())}")
            # → ['@SQ', '@RG', '@PG', '@HD']

            sq = fields.get('@SQ', {})
            print(f"@SQ: {len(sq)} хромосом")
            for name, info in list(sq.items())[:3]:
                
                print(f"  {name}: length={info.get('LN', '?')}")

            rg = fields.get('@RG', {})
            print(f"@RG: {list(rg.keys())}")                    # read groups по ID

            print(f"2. Выравнивания: {count}")
            print(f"3. Хромосомы: {chroms[:3]}")  # первые 3
            print("\n4. Статистика по хромосомам:")
            print(reader.stats_by_chromosome().to_string(index=False))

            if chroms:
                c = chroms[0]
                region_records = reader.get_records_in_region(c, 1, 10_000)
                print(f"\n5. Выравниваний в {c}:1-1000: {len(region_records)}")

            # Первая запись
            for rec in reader.read():
                print(f"6. Первая запись: {rec.chrom}:{rec.start}")
                break
                
    except Exception as e:
        print(f"Ошибка SAM: {e}")


# ====== VCF ======

def demo_vcf(vcf_file):
    """Демонстрирует работу VcfReader на VCF-файле.

    Выводит статистику файла, заголовок, первые варианты,
    распределение по хромосомам и типы вариантов.

    Args:
        vcf_file: Путь к VCF-файлу.

    Examples:
        >>> demo_vcf("variants.vcf")  # doctest: +SKIP
        ========== Анализ VCF ==========
    """
    # vcf_file = input("Введите путь к VCF файлу: ") - убираем input и используем переданный аргумент
    
    if not os.path.exists(vcf_file):
        print(f"Файл {vcf_file} не найден!")
        print("Убедитесь, что файл находится в той же директории")
        return
    
    try:
        print(f"========== Анализ VCF ==========")
        print(f"Выбранный файл: {vcf_file}")
        print("=" * 60)
        
        with VcfReader(vcf_file) as reader:
            # Получаем всю статистику
            print("1. Статистика файла:")
            print("-" * 30)
            
            stats = reader.get_statistics()
            chromosomes = reader.get_chromosomes()
            region_stats = reader.get_region_stats()
            type_stats = reader.get_variant_type_stats()
            
            print(f"   Всего вариантов: {stats.get('total_variants', 0):,}")
            print(f"   Образцы: {stats.get('samples_count', 0)}")
            print(f"   Строк заголовка: {stats.get('header_lines', 0)}")
            print(f"   Хромосомы: {chromosomes}")
            print()

            header = reader.get_header()
            print(f"Группы заголовка: {list(header.keys())}")
            print(f"  ##fileformat: {len(header.get('##fileformat', []))} записей")
            print(f"  ##fileDate {len(header.get('##fileDate', []))} записей")
            print(f"  ##source:   {len(header.get('##source', []))} записей")
            print(f"  ##reference: {len(header.get('##reference', []))} записей")
            print(f"  ##INFO: {len(header.get('##INFO', []))} записей")
            print(f"  #CHROM: {len(header.get('#CHROM', []))} записей")

            # Показываем первые несколько вариантов
            print("2. Первые 5 вариантов:")
            print("-" * 30)
            
            variant_count = 0
            for variant in reader.read():
                if variant_count >= 5:
                    break
                
                # Форматируем информацию для красивого вывода
                info_preview = format_info_preview(variant.info)
                print(f"   {variant.id}")
                print(f"     REF: {variant.ref} -> ALT: {variant.alt}")
                print(f"     INFO: {info_preview}")
                print()
                variant_count += 1
            
            if variant_count == 0:
                print("   Варианты не найдены")
            print()

            # Статистика по регионам
            print("3. Распределение по хромосомам:")
            print("-" * 30)
            
            if not region_stats.empty:
                for _, row in region_stats.iterrows():
                    print(f"   {row['chromosome']}: {row['variant_count']:,} вариантов")
            else:
                print("   Данные по регионам отсутствуют")
            print()

            # Статистика по типам вариантов
            print("4. Типы вариантов:")
            print("-" * 30)
            
            if not type_stats.empty:
                total_variants = type_stats['count'].sum()
                for _, row in type_stats.iterrows():
                    percentage = (row['count'] / total_variants * 100) if total_variants > 0 else 0
                    print(f"   {row['variant_type']}: {row['count']:,} ({percentage:.1f}%)")
            else:
                print("   Данные по типам вариантов отсутствуют")

    except Exception as e:
        print(f"Ошибка при чтении VCF файла: {e}")
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
    print("Поддерживаемые форматы: fasta, fastq, sam, vcf")

    fmt = input("Введите формат файла: ").strip().lower()
    file_path = input("Введите путь к файлу: ").strip().strip('"').strip("'")

    if not os.path.exists(file_path):
        print(f"Файл не найден: {file_path}")
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
    else:
        print(f"Неизвестный формат: {fmt}")