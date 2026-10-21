import os
import sys
from fasta_analyzer import FastaAnalyzer
from abstract import BioParserError

# Список целевых файлов для демонстрации
TARGET_FILES = [
    "covid-genome.fasta",
    "flu_hemagglutinin.fasta",
    "botulinum_neurotoxin.fasta"
]

def print_banner(text: str):
    """Красивый заголовок блока."""
    print("\n" + "=" * 65)
    print(f" {text}")
    print("=" * 65)

def analyze_file(filename: str):
    """Полный анализ одного FASTA-файла по всем требованиям ТЗ."""
    print_banner(f"АНАЛИЗ ФАЙЛА: {filename}")

    # Проверка наличия файла перед началом
    if not os.path.exists(filename):
        print(f"[ОШИБКА]: Файл '{filename}' не найден в текущей папке!")
        print("Пожалуйста, убедитесь, что файл лежит рядом со скриптом.")
        return None

    try:
        with FastaAnalyzer(filename) as analyzer:
            count = analyzer.get_sequences_count()
            print(f"✔ [ТЗ] Количество последовательностей: {count}")

            avg_len = analyzer.get_average_length()
            print(f"✔ [ТЗ] Средняя длина последовательностей: {avg_len:.2f} символов")

            print("\n--- Потоковое чтение записей (через генератор read()) ---")
           
            shown_limit = 3
            for idx, record in enumerate(analyzer, start=1):
                if idx <= shown_limit:
                    is_valid = analyzer.validate_sequence(record.sequence)
                    status_str = "[OK] Валидна (IUPAC)" if is_valid else "[ВНИМАНИЕ] Нестандартные символы"

                    print(f"  Запись #{idx}:")
                    print(f"    • ID:          {record.id}")
                    print(f"    • Описание:    {record.description if record.description else '(нет)'}")
                    print(f"    • Длина:       {len(record)} симв.")
                    print(f"    • Валидность:  {status_str}")
                    # Показываем первые 40 символов
                    preview = record.sequence[:40] + ("..." if len(record) > 40 else "")
                    print(f"    • Начало seq:  {preview}")
                    print("  " + "-" * 40)
                else:
                    if idx == shown_limit + 1:
                        print(f"  ... и ещё {count - shown_limit} записей успешно обработано генератором.")
                    break

            return {
                "file": filename,
                "count": count,
                "avg_len": avg_len
            }

    except BioParserError as err:
        print(f"[КРИТИЧЕСКАЯ ОШИБКА ФОРМАТА]: {err}")
        return None
    except Exception as err:
        print(f"[НЕПРЕДВИДЕННАЯ ОШИБКА]: {err}")
        return None

def main():
    print("Запуск демонстрационной программы для модуля FASTA...")
    print(f"Рабочая директория: {os.getcwd()}")

    summary_stats = []

    for fasta_file in TARGET_FILES:
        res = analyze_file(fasta_file)
        if res:
            summary_stats.append(res)

    if summary_stats:
        print_banner("ИТОГОВАЯ СВОДНАЯ ТАБЛИЦА")
        print(f"{'Файл':<30} | {'Кол-во seq':<12} | {'Ср. длина (bp/aa)':<18}")
        print("-" * 65)
        for stat in summary_stats:
            print(f"{stat['file']:<30} | {stat['count']:<12} | {stat['avg_len']:<18.2f}")
        print("-" * 65)
        print("Анализ всех биологических образцов успешно завершен!")

if __name__ == "__main__":
    main()

