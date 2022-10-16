# x-spline
Реализация x-splines по [статье от 15 сентября 1995 года](https://dl.acm.org/doi/10.1145/218380.218488).
Реализован только двухмерный режим.
## Зависимости
- [openmp](https://www.openmp.org/)
- [matplotlib](https://pypi.org/project/matplotlib/)
## Компиляция
### Release
```bash
cmake -DCMAKE_BUILD_TYPE=Release .
make
```
### Debug
```bash
cmake -DCMAKE_BUILD_TYPE=Debug .
make
```
## Использование
### Запуск
```bash
./x-spline
```
### Визуализация
```python
python visualiser.py
```
