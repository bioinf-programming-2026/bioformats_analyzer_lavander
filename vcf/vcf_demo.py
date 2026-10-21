import os
import sys
from vcf_analyzer import VcfAnalyzer
from abstract import BioParserError

VCF_FILE = "homo_sapiens-chrMT.vcf"


def print_header(title: str):
    """Оформление заголовка блока."""
    print("\n" + "=" * 70)
    print(f" >>> {title}")
    print("=" * 70)

def main():
    print(f"Запуск анализа VCF файла: {VCF_FILE}")

    # 1. Проверка существования файла
    if not os.path.exists(VCF_FILE):
        print(f"\n[ОШИБКА]: Файл '{VCF_FILE}' не найден в текущей папке!")
        print("Пожалуйста, положите файл 'homo_sapiens-chrMT.vcf' рядом с этим скриптом.")
        sys.exit(1)

    try:
        with VcfAnalyzer(VCF_FILE) as analyzer:
            print_header("1. АНАЛИЗ ЗАГОЛОВКА И ГРУПП МЕТАДАННЫХ")
           
            all_headers = analyzer.get_header()
            print(f"Всего строк в заголовке файла: {len(all_headers)}")

            # Извлечение групп: INFO, FILTER, FORMAT
            for group in ["INFO", "FILTER", "FORMAT"]:
                group_lines = analyzer.get_header_group(group)
                print(f"\n• Группа '##{group}' (найдено записей: {len(group_lines)}):")
                # Покажем первые 2 записи для примера
                for line in group_lines[:2]:
                    print(f"   {line}")
                if len(group_lines) > 2:
                    print(f"   ... и ещё {len(group_lines) - 2} строк в группе {group}")

            print_header("2. ПОДСЧЕТ КОЛИЧЕСТВА ВАРИАНТОВ (ГЕНЕРАТОР)")
           
            total_variants = analyzer.get_variants_count()
            print(f"✔ [ТЗ] Общее число мутаций/вариантов в файле: {total_variants}")
            print_header("3. СТАТИСТИКА 'КОЛИЧЕСТВО ВАРИАНТОВ - РЕГИОН' (PANDAS)")
           
            stats_df = analyzer.get_region_statistics()
            print("Таблица распределения вариантов по хромосомам/регионам:")
            print("-" * 40)
            print(stats_df.to_string(index=False))
            print("-" * 40)

            # Определяем точное имя митохондриальной хромосомы в файле (MT, chrM или chrMT)
            target_chrom = stats_df["CHROM"].iloc[0] if not stats_df.empty else "MT"

            start_pos = 1
            end_pos = 1000
           
            print_header(f"4. ПЕРЕСЕЧЕНИЕ РЕГИОНА: {target_chrom}:{start_pos}-{end_pos} (BEDTOOLS INTERSECT)")

            matched_variants = []
            for var in analyzer.intersect(chrom=target_chrom, start=start_pos, end=end_pos):
                matched_variants.append(var)

            print(f"✔ Найдено вариантов в заданном интервале: {len(matched_variants)}")
           
            # Вывод первых 5 найденных мутаций в регионе
            print("\nПримеры найденных мутаций в отрезке:")
            for i, var in enumerate(matched_variants[:5], 1):
                print(f"  [{i}] Позиция: {var.chrom}:{var.pos}")
                print(f"      ID: {var.var_id} | Мутация: {var.ref} -> {','.join(var.alt)}")
                print(f"      QUAL: {var.qual} | FILTER: {var.filter}")
                # Вывод разобранных параметров из поля INFO
                info_preview = ", ".join(f"{k}={v}" for k, v in list(var.info.items())[:3])
                print(f"      INFO: {info_preview} ...")
                print("  " + "-" * 40)

            print_header("5. ПРОВЕРКА ПОТОКОВОГО ЧТЕНИЯ (read)")
            print("Первая запись из генератора без загрузки файла в ОЗУ:")
            first_record = next(iter(analyzer))
            print(f"  Объект:  {first_record}")
            print(f"  Старт:   {first_record.start}")
            print(f"  Конец:   {first_record.end}")
            print(f"  INFO dp: {first_record.info.get('DP', 'нет параметра DP')}")

        print("\n" + "=" * 70)
        print(" Все тесты и требования ТЗ 14.7.2.4 успешно выполнены!")
        print("=" * 70)

    except BioParserError as e:
        print(f"\n[ОШИБКА ФОРМАТА VCF]: {e}")
    except Exception as e:
        print(f"\n[НЕПРЕДВИДЕННАЯ ОШИБКА]: {e}")

if __name__ == "__main__":
    main()