# ecdlp-python

Реализация решения задачи ECDLP на эллиптических кривых над конечным полем
с использованием rho-алгоритма Полларда.

## Возможности

- проверка, что точки лежат на кривой;
- сложение и умножение точек;
- поиск порядка точки;
- решение ECDLP через Pollard rho;
- генерация и прогон тестовых наборов;
- опциональная сверка через SageMath.

## Установка

```bash
python -m venv .venv
.\.venv\Scripts\activate
pip install -e .[dev]
```

## Запуск CLI

```bash
ecdlp
tests
```

## Запуск Tests

```bash
python -c "from ecdlp.test_runner import run_test_file; run_test_file()"
```