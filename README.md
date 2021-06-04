# x-spline

## Зависимости
- [openmp](https://www.openmp.org/)
- [matplotlib](https://pypi.org/project/matplotlib/)
## Компиляция
### Release
```bash
mkdir cmake-build-release
cd cmake-build-release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
### Debug
```bash
mkdir cmake-build-debug
cd cmake-build-debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```
## Визуализация
```python
python visualiser.py
```
