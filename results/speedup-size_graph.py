import matplotlib.pyplot as plt

# Данные
x = [100, 1000, 10000]

y = [0.000956/0.000657, 0.105401/0.038297, 11.458012/3.640418]

y_ideal = [4, 4, 4]

# Создаем график
plt.plot(x, y_ideal, marker=None, linestyle='--', color='black', label = 'Идеальное ускорение')
plt.plot(x, y, marker='o', linestyle='-', color='black', label = 'Реальное ускорение')

# Настройки осей и заголовок
plt.xlabel('Размерность задачи m')
plt.ylabel('Ускорение')

plt.xticks(x)
plt.yticks(y + [4])

plt.grid(True, which="both", ls="-", alpha=0.2)

plt.title('Ускорение для 4 потоков в зависимости от размерности')

plt.legend()

plt.title

# Отображаем график
plt.show()