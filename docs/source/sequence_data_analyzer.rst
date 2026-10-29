Анализатор последовательностей
==============================

.. autoclass:: abstract.SequenceReader
   :members:
   :undoc-members:
   :show-inheritance:
   
.. autoclass:: record.SequenceRecord
   :members:
   :undoc-members:
   :show-inheritance:
   
Анализатор файлов в формате FASTA
---------------------------------
Данный анализатор предназначен для работы с последовательностями в формате FASTA. Он анализирует FASTA файл и выводит следующие данные о последовательности: 

- количество последовательностей

- средняя длина последовательности

.. automodule:: fasta_analyzer
   :members:
   :undoc-members:
   :show-inheritance:
   
.. autofunction:: demo.tests_fasta_analysis

   
Анализатор файлов в формате FASTQ
---------------------------------
Данный анализатор предназначен для работы с последовательностями в формате FASTQ. Он анализирует FASTQ файл и выводит следующие данные о последовательности: 

- количество последовательностей

- средняя длина последовательности

А также строит следующие графики:

- график качества последовательности по нуклеотидам

- график содержания нуклеотидов в последовательностях

- график распределения длин последовательностей

.. automodule:: fastq_analyzer
   :members:
   :undoc-members:
   :show-inheritance:
   
.. autofunction:: demo.tests_fastq_analysis




