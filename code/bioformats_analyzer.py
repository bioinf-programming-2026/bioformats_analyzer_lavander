#!/usr/bin/env python3
"""Интерфейс командной строки для анализатора биологических файлов.

Поддерживаемые форматы: FASTA, FASTQ, SAM, VCF.
Формат определяется автоматически по расширению файла
(или по содержимому, если расширение неизвестно).

Примеры:

    # Автоопределение формата
    python bioformats_analyzer.py  samples.fasta

    # Явное указание формата + запуск тестов
    python bioformats_analyzer.py reads.fq --format fastq --test

    # Сохранить вывод в файл
    python bioformats_analyzer.py variants.vcf --output

    # Указать свой геномный регион для SAM/VCF
    python bioformats_analyzer.py reads.sam --start 1000 --end 5000
"""
import os
import sys
from pathlib import Path

import click

# Настройка sys.path
_current_dir = os.path.dirname(os.path.abspath(__file__))
for _sub in ('', 'fasta', 'fastq', 'sam', 'vcf'):
    _p = os.path.join(_current_dir, _sub) if _sub else _current_dir
    if _p not in sys.path:
        sys.path.insert(0, _p)


# Карта расширений -> формат
_EXT_MAP = {
    '.fasta': 'fasta', '.fa': 'fasta', '.fna': 'fasta',
    '.ffn': 'fasta', '.faa': 'fasta', '.frn': 'fasta',
    '.fastq': 'fastq', '.fq': 'fastq',
    '.sam': 'sam',
    '.vcf': 'vcf',
}


def detect_format(filepath: str) -> str | None:
    """Определяет формат файла по расширению, а при неудаче — по содержимому.

    Args:
        filepath: Путь к файлу.

    Returns:
        Одно из ``'fasta'``, ``'fastq'``, ``'sam'``, ``'vcf'`` или ``None``.

    Examples:
        >>> detect_format("data.fasta")
        'fasta'
        >>> detect_format("data.unknown")
    """
    ext = Path(filepath).suffix.lower()
    if ext in _EXT_MAP:
        return _EXT_MAP[ext]

    # Анализ файла по первой строке
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            first = f.readline().strip()
    except OSError:
        return None

    if first.startswith('##fileformat=VCF'):
        return 'vcf'
    if first.startswith(('@HD', '@SQ', '@RG', '@PG')):
        return 'sam'
    if first.startswith('>'):
        return 'fasta'
    if first.startswith('@'):
        return 'fastq'
    return None

@click.group(
        help='Biological file analysis: FASTA, FASTQ, SAM or VCF.',
        context_settings={'help_option_names': ['-h', '--help']},
        no_args_is_help=True,
)
def main():
    """Анализирует биологический файл (FASTA, FASTQ, SAM или VCF).

    \b
    Примеры:
      python bioformats_analyzer.py samples.fasta
      python bioformats_analyzer.py reads.fq -f fastq --test
      python bioformats_analyzer.py variants.vcf --output
      python bioformats_analyzer.py aln.sam --start 1000 --end 5000
    """
    pass

@main.command(
        name='analyze',
        help='analyze biological file (format is auto-detected).',
)

@click.argument('file', type=click.Path())
@click.option(
    '--format', '-f', 'fmt',
    type=click.Choice(['fasta', 'fastq', 'sam', 'vcf']),
    help='File format (auto-detected if not specified).',
)
@click.option(
    '--test', '-t', is_flag=True,
    help='Run tests before analyzsis (FASTA/FASTQ only).',
)
@click.option(
    '--output', '-o', is_flag=True,
    help='Write output to file (output_<format>.txt).',
)
@click.option(
    '--start', '-s', type=click.IntRange(min=1), default=None,
    help='Genomic region start (SAM/VCF, default 1).',
)
@click.option(
    '--end', '-e', type=click.IntRange(min=1), default=None,
    help='Genomic region end (SAM/VCF, default 10000).',
)
def analyze(file, fmt, test, output, start, end):
    """Analyze biological file (FASTA, FASTQ, SAM или VCF).

    \b
    Examples:
      python bioformats_analyzer.py samples.fasta
      python bioformats_analyzer.py reads.fq -f fastq --test
      python bioformats_analyzer.py variants.vcf --output
      python bioformats_analyzer.py aln.sam --start 1000 --end 5000
    """
    # Проверка файла
    if not os.path.exists(file):
        raise click.ClickException(f"Файл '{file}' не существует")
    if os.path.isdir(file):
        raise click.ClickException(f"'{file}' - папка, а не файл")

    # Проверка интервала пересечения
    

    # Определение формата
    if fmt is None:
        fmt = detect_format(file)
        if fmt is None:
            raise click.UsageError(
                f"Не удалось определить формат файла '{file}'. "
                f"Укажите --format: fasta, fastq, sam или vcf."
            )

    start_given = start is not None
    end_given = end is not None
    test_given = test is True

    # --test применим только к fasta/fastq
    if test_given and fmt in ('sam', 'vcf'):
        click.echo(
            f"Предупреждение: --test не поддерживается для {fmt.upper()} "
            f"(проигнорирован).",
            err=True,
        )
        test = False

    # --start / --end применимы только к sam/vcf
    if fmt in ('fasta', 'fastq'):
        if start_given or end_given:
            click.echo(
                f"Предупреждение: --start/--end не поддерживаются для {fmt.upper()} "
                f"(проигнорированы).",
                err=True,
            )
    else:
        start = 1 if start is None else start
        end = 10000 if end is None else end
        if end < start:
            raise click.ClickException(f"--end ({end}) должен быть >= --start ({start}).")
    
    # Импорт демонстрационных функций
    try:
        from demo import (
            demo_fasta_analysis, tests_fasta_analysis,
            demo_fastq_analysis, tests_fastq_analysis,
            demo_sam, demo_vcf,
        )
    except ImportError as exc:
        raise click.ClickException(f"Не удалось импортировать demo.py: {exc}")

    def run() -> None:
        """Внутренняя функция — сам анализ файла."""
        if fmt == 'fasta':
            if test:
                tests_fasta_analysis()
                click.echo()
            demo_fasta_analysis(file)
        elif fmt == 'fastq':
            if test:
                tests_fastq_analysis()
                click.echo()
            demo_fastq_analysis(file)
        elif fmt == 'sam':
            demo_sam(file, start=start, end=end)
        elif fmt == 'vcf':
            demo_vcf(file, start=start, end=end)

    # Выполнение с перенаправлением stdout при необходимости
    if output:
        output_file = f"output_{fmt}.txt"
        with open(output_file, 'w', encoding='utf-8') as f:
            original_stdout = sys.stdout
            sys.stdout = f
            try:
                run()
            finally:
                sys.stdout = original_stdout
        click.echo(f"Результаты сохранены в {output_file}")
        if fmt == 'fastq':
            click.echo("Графики сохранены в директории test/plots/")
    else:
        run()


if __name__ == '__main__':
    main()