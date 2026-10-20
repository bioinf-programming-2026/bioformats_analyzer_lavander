"""Пример использования FASTQ-анализатора."""

from fastq_analyzer import FastqAnalyzer


def main():
    
    """Прочитать файл, спросив название у юзера, вывести статистику и отфильтровать записи."""
    filename = "test/files/" + input("Введите название файла>")
    analyzer = FastqAnalyzer(filename)
    analyzer.read_records()

    print("Статистика:")

    for name, value in analyzer.get_statistics().items():
        print(name, ":", value)

    print("\nЗаписи:")

    for record in analyzer.records:
        print(
            record.header,
            "Длина:", len(record),
            "GC:", round(record.gc_content(), 2),
            "Качество:", round(record.mean_quality(), 2),
        )

    filtered = analyzer.filter_records(
        min_quality=20,
        min_length=8,
    )

    print("\nПосле фильтрации:", len(filtered))

    # Раскомментируй для сохранения в новый файл:
    # analyzer.save_records(filtered, "test/files/filtered.fastq")


if __name__ == "__main__":
    main()