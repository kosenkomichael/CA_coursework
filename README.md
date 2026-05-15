# ECDLP-python

Реализация решения задачи ECDLP на эллиптических кривых над конечным полем
с использованием rho-алгоритма Полларда.

## Возможности

- проверка, что точки лежат на кривой;
- сложение и умножение точек;
- поиск порядка точки;
- решение ECDLP через Pollard rho;
- генерация и прогон тестовых наборов;
- опциональная сверка через SageMath.

<details>
<summary>Установка и запуск</summary>

## Установка sage-core jupyter notebook

```bash
docker run -p 8888:8888 sagemath/sagemath:10.5 sage-jupyter --JupyterApp.token='' --JupyterApp.password=''
```

## Запуск приложения

```bash
python -m venv .venv
.\.venv\Scripts\activate
pip install -e .[dev]
```

## Запуск CLI

```bash
ecdlp
```

## Запуск Tests

```bash
python -c "from ecdlp.test_runner import run_test_file; run_test_file()"
```
alternative:
```bash
tests
```
</details>

## Структура проекта

```text
ecdlp-python/
├─ pyproject.toml - Файл конфигурации проекта
├─ README.md - Описание проекта, инструкция по установке и запуску.
├─ src/
│  └─ ecdlp/
│     ├─ __init__.py - пакет импортируемый
│     ├─ __main__.py - точка входа при ecdlp
│     ├─ cli.py - интерактивный ввод параметров
│     ├─ models.py - dataclass модели точки и тд
│     ├─ curve.py - реализация класса точки с необходимыми операциями
│     ├─ utils.py - взятие по модулю, решение сравнения и тд
│     ├─ step_calculator.py - итерации алгоритма rho
│     ├─ solver.py - основной модуль для алгоритма, запускает итерации
│     ├─ sage_provider.py - потенциальный запуск через sage
│     └─ test_runner.py - запуск тестов
└─ data/
   └─ ro_alg_ecdlp_data.txt - тестовые данные
```