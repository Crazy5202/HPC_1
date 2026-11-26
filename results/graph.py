import matplotlib.pyplot as plt

# Данные
x = [10, 100, 1000, 2500, 5000, 7500, 10000]
y_usual = [0.000101, 0.000961, 0.099199, 0.644105, 2.662025, 6.086427, 11.054989]

y_mpi_1 = [0.000124, 0.000956, 0.105401, 0.673295, 2.731849, 6.269190, 11.458012]
y_mpi_2 = [0.000122, 0.000706, 0.061207, 0.375164, 1.520324, 3.419562, 6.230278]
y_mpi_4 = [0.000444, 0.000657, 0.038297, 0.229380, 0.907373, 2.045004, 3.640418]

# Создаем график
#plt.plot(x, y_usual, marker='o', linestyle='--', color='b', label = 'Обычная реализация')
plt.plot(x, y_mpi_1, marker='o', linestyle='-', color='r', label = '1 поток')
plt.plot(x, y_mpi_2, marker='o', linestyle='-', color='b', label = '2 потока')
plt.plot(x, y_mpi_4, marker='o', linestyle='-', color='g', label = '4 потока')

# Настройки осей и заголовок
plt.xlabel('Размер матрицы и вектора m')
plt.ylabel('Время выполнения, с')

plt.grid(True, which="both", ls="-", alpha=0.2)

plt.xticks([100, 1000, 2500, 5000, 7500, 10000])

plt.legend()

# Отображаем график
plt.show()